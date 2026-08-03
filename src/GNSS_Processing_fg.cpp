/*
 * BSD 3-Clause License

 *  Copyright (c) 2025, Dongjiao He
 *  All rights reserved.
 *
 *  Author: Dongjiao HE <hdj65822@connect.hku.hk>
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the Universitaet Bremen nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#include "GNSS_Processing_fg.h"
#include "LambdaAmbiguityResolver.h"
#include "RtkSignalUtils.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <set>
#include <sstream>

namespace
{
const char *rtkBandName(RtkSignalBand band)
{
  return band == RtkSignalBand::Primary ? "L1/primary" :
         band == RtkSignalBand::Secondary ? "L2/E5b/B2" :
         band == RtkSignalBand::L5 ? "L5/E5a/B2a" :
         band == RtkSignalBand::Extra ? "E6/B3/G3" :
         band == RtkSignalBand::Wide ? "E5ab/B2ab" :
         band == RtkSignalBand::WideLane ? "wide-lane" :
         band == RtkSignalBand::NarrowLane ? "narrow-lane" : "unsupported";
}

const char *rtkSystemName(uint32_t system)
{
  return system == SYS_GPS ? "GPS" :
         system == SYS_GLO ? "GLONASS" :
         system == SYS_GAL ? "Galileo" :
         system == SYS_BDS ? "BeiDou" : "UNKNOWN";
}

bool isRtkDiagnosticEpoch(size_t epoch, int interval)
{
  return epoch <= 3 || interval <= 1 ||
         epoch % static_cast<size_t>(interval) == 0;
}

EphemBasePtr closestVisualizationEphemeris(const GNSSAssignment &assignment,
                                           const ObsPtr &observation)
{
  if (!observation) return EphemBasePtr();
  const auto &ephemerides = assignment.ephem_from_rinex
      ? assignment.sat2ephem_rnx : assignment.sat2ephem;
  const auto &time_indices = assignment.ephem_from_rinex
      ? assignment.sat2time_index_rnx : assignment.sat2time_index;
  const auto ephemeris_it = ephemerides.find(observation->sat);
  const auto time_it = time_indices.find(observation->sat);
  if (ephemeris_it == ephemerides.end() || time_it == time_indices.end())
    return EphemBasePtr();

  const double observation_time = time2sec(observation->time);
  double closest_delta = EPH_VALID_SECONDS;
  size_t closest_index = 0;
  bool found = false;
  for (const auto &entry : time_it->second)
  {
    const double delta = std::abs(entry.first - observation_time);
    if (delta < closest_delta)
    {
      closest_delta = delta;
      closest_index = entry.second;
      found = true;
    }
  }
  if (!found || closest_index >= ephemeris_it->second.size())
    return EphemBasePtr();
  return ephemeris_it->second[closest_index];
}

Eigen::Vector3d visualizationSatelliteEcef(const ObsPtr &observation,
                                           const EphemBasePtr &ephemeris)
{
  if (!observation || !ephemeris) return Eigen::Vector3d::Zero();
  if (satsys(observation->sat, nullptr) == SYS_GLO)
  {
    const GloEphemPtr glo_ephemeris =
        std::dynamic_pointer_cast<GloEphem>(ephemeris);
    return glo_ephemeris
        ? geph2pos(observation->time, glo_ephemeris, nullptr)
        : Eigen::Vector3d::Zero();
  }
  const EphemPtr broadcast_ephemeris =
      std::dynamic_pointer_cast<Ephem>(ephemeris);
  return broadcast_ephemeris
      ? eph2pos(observation->time, broadcast_ephemeris, nullptr)
      : Eigen::Vector3d::Zero();
}
}  // namespace

// Processing order:
//   processGNSS -> GNSSLIAlign (until initialized)
//   Evaluate -> initializeFrameEstimate -> AddFactor -> runISAM2opt
//            -> attemptIntegerAmbiguityResolution -> updateStateFromEstimate

// --------------------------------------------------------------------------
// Lifecycle and input-buffer management
// --------------------------------------------------------------------------

GNSSProcess::GNSSProcess()
    : diff_t_gnss_local(0.0)
{
  Reset();
  // initNoises();
}

GNSSProcess::~GNSSProcess() {}

void GNSSProcess::clearGnssBuffer(size_t index)
{
  gnss_meas_buf[index].clear();
  gnss_ephem_buf[index].clear();
}

void GNSSProcess::shiftInitializationWindow(size_t first_index)
{
  // Retain only the continuous suffix after a bad initialization epoch/gap.
  const size_t retained = wind_size + 1 - first_index;
  for (size_t destination = 0; destination < retained; ++destination)
  {
    const size_t source = destination + first_index;
    gnss_meas_buf[destination] = std::move(gnss_meas_buf[source]);
    gnss_ephem_buf[destination] = std::move(gnss_ephem_buf[source]);
    rot_window[destination] = rot_window[source];
    pos_window[destination] = pos_window[source];
    vel_window[destination] = vel_window[source];
    pos_ecef_window[destination] = pos_ecef_window[source];
  }
  frame_count = static_cast<int>(retained);
  for (size_t index = retained; index < static_cast<size_t>(wind_size + 1); ++index)
  {
    clearGnssBuffer(index);
  }
}

int GNSSProcess::firstValidClockBias(const Eigen::Matrix<double, 7, 1> &rough_xyzt) const
{
  for (int index = 0; index < 4; ++index)
  {
    if (std::abs(rough_xyzt(3 + index)) > 0.0)
    {
      return index;
    }
  }
  return -1;
}

void GNSSProcess::Reset()
{
  ROS_WARN("Reset GNSSProcess");
  // Clear both the legacy time-differenced carrier history and GNSS window.
  std::map<sat_first, std::map<uint32_t, double[6]>> empty_map_c;
  sat2cp.swap(empty_map_c);
  // sat2time_index.swap(empty_map_i);
  // sat2ephem.swap(empty_map_e);
  for (size_t i = 0; i < WINDOW_SIZE+1; i++)
  {
    std::vector<ObsPtr> empty_vec_o;
    std::vector<EphemBasePtr> empty_vec_e;
    gnss_meas_buf[i].swap(empty_vec_o);
    gnss_ephem_buf[i].swap(empty_vec_e);
  }
  p_assign->change_ext = 1;
  std::map<uint32_t, uint32_t> empty_map_t;
  std::map<uint32_t, double> empty_map_st;
  p_assign->sat_track_status.swap(empty_map_t);
  p_assign->sat_track_time.swap(empty_map_st);
  p_assign->sat_track_last_time.swap(empty_map_st);
  p_assign->hatch_filter_meas.swap(empty_map_st);
  p_assign->last_cp_meas.swap(empty_map_st);
  p_assign->gtSAMgraph.resize(0);
  p_assign->initialEstimate.clear();
  p_assign->isamCurrentEstimate.clear();
  p_assign->sum_d = 0;
  p_assign->sum_d2 = 0;
  // p_assign->hatch_filter_meas = 0;
  // p_assign->last_cp = 0;
  // index_delete = 0;
  frame_delete = 0;
  p_assign->factor_id_frame.clear();
  id_accumulate = 0;
  ambiguities_.clear();
  reference_satellites_.clear();
  raw_satellite_directions_.clear();
  valid_satellite_directions_.clear();
  rtk_satellite_directions_.clear();
  next_ambiguity_id_ = 0;
  rtk_dd_added_this_epoch_ = false;
  integer_solution_available = false;
  last_lambda_ratio = 0.0;
  rtk_synchronized_epochs = 0;
  rtk_base_epoch_misses = 0;
  rtk_zero_factor_epochs = 0;
  rtk_double_difference_factors = 0;
  rtk_double_difference_pseudorange_factors = 0;
  rtk_secondary_factors = 0;
  rtk_l5_factors = 0;
  rtk_fix_attempts = 0;
  rtk_fix_successes = 0;
  rtk_float_ambiguity_count = 0;
  rtk_fixed_ambiguity_count = 0;
  rtk_fix_status = "NOT_ATTEMPTED";
  frame_num = 0;
  last_gnss_time = 0.0;
  first_gnss_time = 0.0;
  frame_count = 0;
  invalid_lidar = false;
  Rot_gnss_init.setIdentity();
  p_assign->process_feat_num = 0;
  gnss_ready = false;
  {
    // Reinitialize the IMU delta with zero biases; optimized biases will be
    // injected after the first no-LiDAR iSAM2 update.
    pre_integration->repropagate(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());
  }

  // iSAM2 is recreated because Reset() starts an entirely new graph, rather
  // than continuing to relinearize variables from the previous trajectory.
  gtsam::ISAM2Params parameters;
  parameters.relinearizeThreshold = 0.1;
  parameters.relinearizeSkip = 5; // may matter? improtant!
  p_assign->isam = gtsam::ISAM2(parameters);
}

void GNSSProcess::inputIonoParams(double ts, const std::vector<double> &iono_params)
{
  // The broadcast Klobuchar-style model is represented by exactly 8 values.
  if (iono_params.size() != 8)    return;

  // update ionosphere parameters
  std::vector<double> empty_vec_d;
  p_assign->latest_gnss_iono_params.swap(empty_vec_d);
  std::copy(iono_params.begin(), iono_params.end(), std::back_inserter(p_assign->latest_gnss_iono_params));
}

void GNSSProcess::inputpvt(double ts, double lat, double lon, double alt, int float_sol, int diff_sol)
{
  Eigen::Vector3d lla;
  lla << lat, lon, alt;
  if (pvt_time.empty())
  {
    first_lla_pvt = lla;
    first_xyz_ecef_pvt = geo2ecef(lla);
    // first_lla_pvt << 22.609671,114.017229,98.401000;
    // first_xyz_ecef_pvt << -2397684.725162,5380932.949883,2436910.600325;
    printf("first ecef xyz 3:%f,%f,%f\n",first_xyz_ecef_pvt(0),first_xyz_ecef_pvt(1),first_xyz_ecef_pvt(2));
    printf("first lla:%f,%f,%f\n",first_lla_pvt(0),first_lla_pvt(1),first_lla_pvt(2));
  }
  Eigen::Vector3d xyz_ecef = geo2ecef(lla);
  Eigen::Vector3d xyz_enu = ecef2enu(first_lla_pvt, xyz_ecef - first_xyz_ecef_pvt);
  pvt_time.push_back(ts);
  pvt_holder.push_back(xyz_enu);
  diff_holder.push_back(diff_sol);
  float_holder.push_back(float_sol);
}

void GNSSProcess::inputlla(double ts, double lat, double lon, double alt) //
{
  Eigen::Vector3d lla;
  lla << lat, lon, alt;
  if (lla_time.empty())
  {
    first_lla_lla = lla;
    first_xyz_ecef_lla = geo2ecef(lla);
  }
  Eigen::Vector3d xyz_ecef = geo2ecef(lla);
  Eigen::Vector3d xyz_enu = ecef2enu(first_lla_lla, xyz_ecef - first_xyz_ecef_lla);
  lla_time.push_back(ts);
  lla_holder.push_back(xyz_enu);
}

Eigen::Vector3d GNSSProcess::local2enu(Eigen::Matrix3d R_enu_local_, Eigen::Vector3d anc, Eigen::Vector3d &pos)
{
  Eigen::Vector3d enu_pos;
  if (!nolidar)
  {
    enu_pos = R_enu_local_ * pos; // - anc_local); //

    // Eigen::Matrix3d R_ecef_enu_ = ecef2rotation(anc);
    // Eigen::Vector3d ecef_pos_ = anc + R_ecef_enu_ * enu_pos;
    Eigen::Vector3d ecef_pos_ = anc + enu_pos;
    // Eigen::Vector3d lla_pos = ecef2geo(first_xyz_enu_pvt);
    enu_pos = ecef2enu(first_lla_pvt, ecef_pos_ - first_xyz_ecef_pvt);
  }
  else
  {
    Eigen::Vector3d pos_r = pos;
    // Eigen::Vector3d lla_pos = ecef2geo(first_xyz_enu_pvt);
    enu_pos = ecef2enu(first_lla_pvt, pos_r - first_xyz_ecef_pvt);
  }
  return enu_pos;
}

void GNSSProcess::inputGNSSTimeDiff(const double t_diff) //
{
    diff_t_gnss_local = t_diff;
}

bool GNSSProcess::loadBaseStationFile(const std::string &path)
{
  std::string error;
  use_double_differences =
      base_station_data_.load(path, error, base_carrier_std_cycles);
  if (!use_double_differences)
  {
    ROS_ERROR_STREAM("Base-station GNSS file rejected: " << error);
    return false;
  }
  ROS_INFO_STREAM(
      "[RTK-TIME] loaded base observations path=" << path
      << " rinex_time_system=" << base_station_data_.timeSystem()
      << " epochs=" << base_station_data_.epochCount()
      << " gpst_range=[" << std::fixed << std::setprecision(3)
      << base_station_data_.firstTimestamp() << ','
      << base_station_data_.lastTimestamp() << ']'
      << " base_ecef_m=[" << base_station_data_.ecefPosition().transpose()
      << "]");
  return true;
}

void GNSSProcess::processGNSS(const std::vector<ObsPtr> &gnss_meas, state_output &state)
{
  static size_t received_epoch_count = 0;
  ++received_epoch_count;
  // processGNSSBase performs observation quality filtering, ephemeris
  // association, satellite tracking, and Hatch pseudorange smoothing.
  std::vector<double>().swap(psr_meas_hatch_filter);
  std::vector<ObsPtr> valid_meas;
  std::vector<EphemBasePtr> valid_ephems;
  valid_ephems.clear();
  valid_meas.clear();
  if (gnss_meas.empty())
  {
    raw_satellite_directions_.clear();
    valid_satellite_directions_.clear();
    if (gnss_ready)
    {
      clearGnssBuffer(0);
    }
    return;
  }

  if (gnss_ready)
  {
    // Convert the IMU origin to the GNSS antenna point before updating the
    // approximate ECEF position used by measurement preprocessing.
    Eigen::Vector3d pos_gnss = state.pos + state.rot * Tex_imu_r; // .normalized()
    updateGNSSStatistics(pos_gnss);
  }
  p_assign->processGNSSBase(gnss_meas, psr_meas_hatch_filter, valid_meas, valid_ephems, gnss_ready, ecef_pos, last_gnss_time);
  ROS_INFO_STREAM_THROTTLE(
      5.0, "GNSS/RTK runtime status: received epochs=" << received_epoch_count
           << ", input satellites=" << gnss_meas.size()
           << ", valid satellites=" << valid_meas.size()
           << ", graph initialized=" << (gnss_ready ? "yes" : "no")
           << ", synchronized base epochs=" << rtk_synchronized_epochs
           << ", DD factors=" << rtk_double_difference_factors);

  if (!gnss_ready)
  {
    // At least four satellites are required for coarse position/clock
    // initialization. Store synchronized local motion for frame alignment.
    if (valid_meas.size() < 4) return; // right or not? valid_meas.empty() ||
    {
      rot_window[frame_count] = state.rot; //.normalized().toRotationMatrix();
      pos_window[frame_count] = state.pos + state.rot * Tex_imu_r; // .normalized()
      Eigen::Matrix3d omg_skew;
      omg_skew << SKEW_SYM_MATRX(state.omg);
      vel_window[frame_count] = state.vel + state.rot * omg_skew * Tex_imu_r; // .normalized().toRotationMatrix()
      // vel_window[frame_count] = state.vel;
    }
    gnss_meas_buf[frame_count] = valid_meas;
    gnss_ephem_buf[frame_count] = valid_ephems;
    frame_count ++;
    gnss_ready = GNSSLIAlign();
    if (gnss_ready)
    {
      ROS_INFO("GNSS Initialization is done");
      state_const_ = state;
      // state_const_last = state;
    }
  }
  else
  {
    // After initialization only slot zero is used as the current epoch.
    gnss_meas_buf[0] = valid_meas;
    gnss_ephem_buf[0] = valid_ephems;
  }

  if (gnss_ready && !nolidar)
  {
    updateRoverSatelliteDirections(gnss_meas, valid_meas, valid_ephems,
                                   state);
  }
}

// --------------------------------------------------------------------------
// Incremental optimization and fixed-lag cleanup
// --------------------------------------------------------------------------

void GNSSProcess::runISAM2opt(void) //
{
  gtsam::FactorIndices delete_factor;
  gtsam::FactorIndices().swap(delete_factor);

  if (gnss_ready)
  {
    bool delete_happen = false;
    if (frame_num - frame_delete > delete_thred) // (graph_whole1.size() - index_delete > 4000)
    {
      // Collect all factors owned by frames leaving the active graph window.
      delete_happen = true;
      while (frame_num - frame_delete > delete_thred) // (graph_whole1.size() - index_delete > 3000)
      {
        if (!p_assign->factor_id_frame.empty())
        {
          // if (frame_delete > 0)
          {
          for (size_t i = 0; i < p_assign->factor_id_frame[0].size(); i++)
          {
            {
              delete_factor.push_back(p_assign->factor_id_frame[0][i]);
            }
          }
          // index_delete += p_assign->factor_id_frame[0].size();
          }

          p_assign->factor_id_frame.pop_front();
          frame_delete ++;
        }
        if (p_assign->factor_id_frame.empty()) break;
      }
    }

    if (delete_happen)
    {
      // GNSSAssignment remaps factor identifiers while marginalizing/removing
      // old variables, so this path replaces the ordinary update below.
      p_assign->delete_variables(nolidar, frame_delete, frame_num, id_accumulate, delete_factor);
    }
    else
    {
      p_assign->isam.update(p_assign->gtSAMgraph, p_assign->initialEstimate);
      p_assign->gtSAMgraph.resize(0); // will the initialEstimate change?
      p_assign->initialEstimate.clear();
      p_assign->isam.update();
    }
  }
  else
  {
    p_assign->isam.update(p_assign->gtSAMgraph, p_assign->initialEstimate);
    p_assign->gtSAMgraph.resize(0); // will the initialEstimate change?
    p_assign->initialEstimate.clear();
    p_assign->isam.update();
  }
  p_assign->isamCurrentEstimate = p_assign->isam.calculateEstimate();

  if (nolidar) // || invalid_lidar)
  {
    // Continue IMU integration from the newly optimized accelerometer and
    // gyroscope biases.
    pre_integration->repropagate(p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(6),
                                p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(9));
  }
  pruneCarrierPhaseHistoryByFrame();
}

void GNSSProcess::pruneCarrierPhaseHistoryByFrame()
{
  while (!sat2cp.empty() && sat2cp.begin()->first.frame_num < frame_delete)
  {
    sat2cp.erase(sat2cp.begin());
  }
}

void GNSSProcess::pruneCarrierPhaseHistoryByTime(double time_current)
{
  while (!sat2cp.empty() &&
         time_current - sat2cp.begin()->first.timecur > gnss_cp_time_threshold)
  {
    sat2cp.erase(sat2cp.begin());
  }
}

// --------------------------------------------------------------------------
// GNSS/local-frame initialization
// --------------------------------------------------------------------------

bool GNSSProcess::GNSSLIAlign()
{
  // A full motion/GNSS window is needed before estimating a stable global
  // anchor and local-to-ECEF rotation.
  if (frame_count < wind_size + 1) return false;

  for (uint32_t i = 0; i < wind_size; i++)
  {
    if (time2sec(gnss_meas_buf[i+1][0]->time) - time2sec(gnss_meas_buf[i][0]->time) > 15 * gnss_sample_period) // need IMU to prop
    {
      // Discard the discontinuous prefix but keep later usable epochs.
      // if (frame_count == wind_size + 1)
      // {
        shiftInitializationWindow(i + 1);
      // }
      return false;
    }
  }

  std::vector<std::vector<ObsPtr>> curr_gnss_meas_buf;
  std::vector<std::vector<EphemBasePtr>> curr_gnss_ephem_buf;
  if ((pos_window[wind_size] - pos_window[0]).norm() < 10.0) // && pos_window[0].norm() < 10.0)
  {
    // Low-excitation initialization: use coarse SPP for the ECEF anchor and
    // the conventional ENU rotation because trajectory alignment is weak.
    for (uint32_t i = 0; i < (wind_size+1); ++i)
    {
        curr_gnss_meas_buf.push_back(gnss_meas_buf[i]);
        curr_gnss_ephem_buf.push_back(gnss_ephem_buf[i]);
    }

    GNSSLIInitializer gnss_li_initializer(curr_gnss_meas_buf, curr_gnss_ephem_buf, p_assign->latest_gnss_iono_params);

    // 1. get a rough global location
    Eigen::Matrix<double, 7, 1> rough_xyzt;
    // Eigen::Matrix<double, 3, 1> rough_xyz;
    rough_xyzt.setZero();
    // rough_xyz.setZero();
    if (!gnss_li_initializer.coarse_localization(rough_xyzt))
    {
        std::cerr << "Fail to obtain a coarse location.\n";
        for (uint32_t i = 0; i < (wind_size); ++i)
        {
          gnss_meas_buf[i] = gnss_meas_buf[i+1];
          gnss_ephem_buf[i] = gnss_ephem_buf[i+1];
          rot_window[i] = rot_window[i+1];
          pos_window[i] = pos_window[i+1];
          vel_window[i] = vel_window[i+1];
        }
        frame_count = wind_size;
        std::vector<ObsPtr> empty_vec_o;
        std::vector<EphemBasePtr> empty_vec_e;
        gnss_meas_buf[frame_count].swap(empty_vec_o);
        gnss_ephem_buf[frame_count].swap(empty_vec_e);
        return false;
    }
    {
      const int dt_idx = firstValidClockBias(rough_xyzt);
      if (dt_idx < 0)
      {
        std::cerr << "Fail to quick init anchor point.\n";
        for (uint32_t i = 0; i < (wind_size); ++i)
        {
          gnss_meas_buf[i] = gnss_meas_buf[i+1]; // change the strategy
          gnss_ephem_buf[i] = gnss_ephem_buf[i+1];

          rot_window[i] = rot_window[i+1];
          pos_window[i] = pos_window[i+1];
          vel_window[i] = vel_window[i+1];
        }
        frame_count = wind_size;
        std::vector<ObsPtr> empty_vec_o;
        std::vector<EphemBasePtr> empty_vec_e;
        gnss_meas_buf[frame_count].swap(empty_vec_o);
        gnss_ephem_buf[frame_count].swap(empty_vec_e);
        return false;
      }
      anc_local = pos_window[0] - rot_window[0] * Tex_imu_r; // [wind_size]; // ?
      yaw_enu_local = 0.0; // -2418165.665753, 5385967.410215, 2405315.115443; //
      para_rcv_ddt[0] = 0.0; // 128.0;
      // rough_xyz = rough_xyzt.head<3>();
      // if (anc_local.norm() > 100)
      // {
      //   std::vector<Eigen::Vector3d> local_vs;
      //   for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
      //       local_vs.push_back(vel_window[i]); // values at gnss measurement
      //   if (gnss_li_initializer.yaw_alignment(local_vs, rough_xyz, yaw_enu_local, para_rcv_ddt[0]))
      //   {
      //     printf("yaw_enu_local:%f\n",yaw_enu_local);
      //   }
      //   else
      //   {
      //     yaw_enu_local = 0.0;
      //     para_rcv_ddt[0] = 0.0;
      //     rough_xyz = rough_xyzt.head<3>();
      //   }
      // }
      anc_ecef = rough_xyzt.head<3>(); // - anc_local; << -2418181.50, 5385962.29, 2405305.18;
      R_ecef_enu = ecef2rotation(anc_ecef); // * Eigen::AngleAxisd(yaw_enu_local, Eigen::Vector3d::UnitZ()).matrix(); // * yawAngle; // * pitchAngle * rollAngle; //<< 0.772234, 0.501306, -0.390316,
                      // 0.047633, 0.566933, 0.822386,
                      // 0.633550, -0.653666, 0.413926; //
      anc_ecef -= R_ecef_enu * anc_local; // anc_local too large: need initialize yaw
      // R_ecef_enu = ecef2rotation(anc_ecef); // * Eigen::AngleAxisd(yaw_enu_local, Eigen::Vector3d::UnitZ()).matrix(); // * yawAngle; // * pitchAngle * rollAngle; //<< 0.772234, 0.501306, -0.390316,
      para_rcv_dt[4*wind_size] = rough_xyzt(3+dt_idx);
    }
  }
  else
  {
    // Excited-motion initialization: obtain an ECEF point for every local
    // point and align the two trajectories to estimate the global transform.
    for (uint32_t i = 0; i < (wind_size+1); ++i)
    {
      // for (uint32_t j = 0; j < gnss_meas_buf[i].size(); j++)
      // {
      //   int freq_idx_ = -1;
      //   double freq = L1_freq(gnss_meas_buf[i][j], &freq_idx_); // L1_freq NEEDED
      //   // if (freq_idx_ < 0)   continue;
      // }
        curr_gnss_meas_buf.push_back(gnss_meas_buf[i]);
        curr_gnss_ephem_buf.push_back(gnss_ephem_buf[i]);
      // }

      GNSSLIInitializer gnss_li_initializer(curr_gnss_meas_buf, curr_gnss_ephem_buf, p_assign->latest_gnss_iono_params);

      // 1. get a rough global location
      Eigen::Matrix<double, 7, 1> rough_xyzt;
      Eigen::Matrix<double, 3, 1> rough_xyz;
      rough_xyzt.setZero();
      rough_xyz.setZero();
      if (!gnss_li_initializer.coarse_localization(rough_xyzt))
      {
        for (uint32_t j = i; j < wind_size; ++j)
        {
          gnss_meas_buf[j] = gnss_meas_buf[j+1];
          gnss_ephem_buf[j] = gnss_ephem_buf[j+1];
          rot_window[j] = rot_window[j+1];
          pos_window[j] = pos_window[j+1];
          vel_window[j] = vel_window[j+1];
        }
        frame_count -= 1;
        for (uint32_t j = frame_count; j < wind_size+1; ++j) // wind_size-i
        {
          std::vector<ObsPtr> empty_vec_o;
          std::vector<EphemBasePtr> empty_vec_e;
          gnss_meas_buf[j].swap(empty_vec_o);
          gnss_ephem_buf[j].swap(empty_vec_e);
        }
        return false;
      }
      const int dt_idx = firstValidClockBias(rough_xyzt);
      // if (num_dt_fail == 4)
      if (dt_idx < 0)
      {
        std::cerr << "Fail to quick init anchor point.\n";
        for (uint32_t j = i; j < wind_size; ++j)
        {
          gnss_meas_buf[j] = gnss_meas_buf[j+1]; // change the strategy
          gnss_ephem_buf[j] = gnss_ephem_buf[j+1];

          rot_window[j] = rot_window[j+1];
          pos_window[j] = pos_window[j+1];
          vel_window[j] = vel_window[j+1];
        }
        frame_count = wind_size;
        std::vector<ObsPtr> empty_vec_o;
        std::vector<EphemBasePtr> empty_vec_e;
        gnss_meas_buf[frame_count].swap(empty_vec_o);
        gnss_ephem_buf[frame_count].swap(empty_vec_e);
        return false;
      }
      pos_ecef_window[i] = rough_xyzt.head<3>();
      para_rcv_dt[4*i] = rough_xyzt(3+dt_idx);
    }

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_in (new pcl::PointCloud<pcl::PointXYZ>(wind_size+1, 1));
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_out (new pcl::PointCloud<pcl::PointXYZ>(wind_size+1, 1));

    // Fill in the CloudIn data
    // for (auto& point : *cloud_in)
    for (size_t i = 0; i < wind_size + 1; i++)
    {
      cloud_in->points[i].x = pos_window[i](0);
      cloud_in->points[i].y = pos_window[i](1);
      cloud_in->points[i].z = pos_window[i](2);
      cloud_out->points[i].x = pos_ecef_window[i](0);
      cloud_out->points[i].y = pos_ecef_window[i](1);
      cloud_out->points[i].z = pos_ecef_window[i](2);
    }
    pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;
    icp.setInputSource(cloud_in);
    icp.setInputTarget(cloud_out);

    pcl::PointCloud<pcl::PointXYZ> Final;
    icp.align(Final);

    std::cout << "ICP has " << (icp.hasConverged()?"converged":"not converged") << ", score: " << icp.getFitnessScore() << std::endl;
    // std::cout << icp.getFinalTransformation() << std::endl;
    // Eigen::Vector3d pos_trans;
    // Eigen::Matrix3d rot_trans;
    // TrajAlign(local_traj, enu_traj, pos_trans, rot_trans);
    if (!icp.hasConverged())
    {
        std::cerr << "Fail to obtain a coarse location.\n";
        for (uint32_t i = 0; i < (wind_size); ++i)
        {
          gnss_meas_buf[i] = gnss_meas_buf[i+1];
          gnss_ephem_buf[i] = gnss_ephem_buf[i+1];
          rot_window[i] = rot_window[i+1];
          pos_window[i] = pos_window[i+1];
          vel_window[i] = vel_window[i+1];
          pos_ecef_window[i] = pos_ecef_window[i+1];
          para_rcv_dt[4*i] = para_rcv_dt[4*(i+1)];
        }
        frame_count = wind_size;
        std::vector<ObsPtr> empty_vec_o;
        std::vector<EphemBasePtr> empty_vec_e;
        gnss_meas_buf[frame_count].swap(empty_vec_o);
        gnss_ephem_buf[frame_count].swap(empty_vec_e);
        return false;
    }

    Eigen::Matrix4d sim_trans = icp.getFinalTransformation().cast<double>();
    // Eigen::Vector3d pos_nmea(nmea_meas_[0]->pose.pose.position.x, nmea_meas_[0]->pose.pose.position.y, nmea_meas_[0]->pose.pose.position.z);
    // Eigen::Vector3d pos_nmea(nmea_meas_[0]->pose.pose.position.x, nmea_meas_[0]->pose.pose.position.y, nmea_meas_[0]->pose.pose.position.z);
    anc_ecef = sim_trans.block<3, 1>(0, 3); // icp.getFinalTransformation().template block<3, 1>(0, 3); // pos_trans; // pos_nmea - pos_window[0]; //

    anc_local = Eigen::Vector3d::Zero(); // pos_window[0]; // [WINDOW_SIZE]; // ?
    yaw_enu_local = 0.0; // -2418165.665753, 5385967.410215, 2405315.115443; //
    para_rcv_ddt[0] = 0.0; // 128.0;

    R_ecef_enu = sim_trans.block<3, 3>(0, 0); // ecef2rotation(anc_ecef) * Eigen::AngleAxisd(yaw_enu_local, Eigen::Vector3d::UnitZ()).matrix(); // * yawAngle; // * pitchAngle * rollAngle; //<< 0.772234, 0.501306, -0.390316,
                      // 0.047633, 0.566933, 0.822386,
                      // 0.633550, -0.653666, 0.413926; //
  }
    // Seed graph priors only after a valid anchor/alignment has been found.
    SetInit();
    frame_num = 1; // frame_count;
    last_gnss_time = time2sec(gnss_meas_buf[wind_size][0]->time);
    first_gnss_time = time2sec(gnss_meas_buf[wind_size][0]->time);
    // printf("first gnss time: %f", first_gnss_time);
  // }

  for (uint32_t k_ = 1; k_ < wind_size+1; k_++)
  {
    std::vector<ObsPtr>().swap(gnss_meas_buf[k_]);
    std::vector<EphemBasePtr>().swap(gnss_ephem_buf[k_]);
  }
  runISAM2opt();
  return true;
}

// --------------------------------------------------------------------------
// Satellite-state preparation and per-epoch graph orchestration
// --------------------------------------------------------------------------

void GNSSProcess::updateGNSSStatistics(Eigen::Vector3d &pos) // delete
{
  if (!nolidar)
  {
    Eigen::Vector3d anc_cur;
    Eigen::Matrix3d R_enu_local_;
    // if (frame_num == 1)
    // {
    //   anc_cur = anc_ecef;
    //   R_enu_local_ = R_ecef_enu * Eigen::AngleAxisd(yaw_enu_local, Eigen::Vector3d::UnitZ()) * Rot_gnss_init;
    // }
    // else
    {
      anc_cur = p_assign->isamCurrentEstimate.at<gtsam::Vector3>(E(0));
      R_enu_local_ = p_assign->isamCurrentEstimate.at<gtsam::Rot3>(P(0)).matrix();
    }
    Eigen::Vector3d enu_pos = R_enu_local_ * pos; // - anc_local);
    // R_ecef_enu = ecef2rotation(anc_cur);
    // ecef_pos = anc_cur + R_ecef_enu * enu_pos;
    ecef_pos = anc_cur + enu_pos;
  }
  else
  {
    ecef_pos = pos;
  }
}

void GNSSProcess::GnssPsrDoppMeas(const ObsPtr &obs_, const EphemBasePtr &ephem_)
{
  // Pseudorange/Doppler still follows the project's primary-frequency model.
  // L5 is added separately by addDoubleDifferenceFactors().
  freq = L1_freq(obs_, &freq_idx);
  LOG_IF(FATAL, freq < 0) << "No L1 observation found.";

  uint32_t sys = satsys(obs_->sat, NULL);
  double tof = obs_->psr[freq_idx] / LIGHT_SPEED;
  gtime_t sv_tx = time_add(obs_->time, -tof);

  if (sys == SYS_GLO)
  {
      GloEphemPtr glo_ephem = std::dynamic_pointer_cast<GloEphem>(ephem_);
      svdt = geph2svdt(sv_tx, glo_ephem);
      sv_tx = time_add(sv_tx, -svdt);
      sv_pos = geph2pos(sv_tx, glo_ephem, &svdt);
      sv_vel = geph2vel(sv_tx, glo_ephem, &svddt);
      tgd = 0.0;
      pr_uura = 2.0 * (obs_->psr_std[freq_idx]/0.16);
      dp_uura = 2.0 * (obs_->dopp_std[freq_idx]/0.256);
  }
  else
  {
      EphemPtr eph = std::dynamic_pointer_cast<Ephem>(ephem_);
      svdt = eph2svdt(sv_tx, eph); // used in eva
      sv_tx = time_add(sv_tx, -svdt);
      sv_pos = eph2pos(sv_tx, eph, &svdt); // used in eva
      sv_vel = eph2vel(sv_tx, eph, &svddt); // used in eva
      tgd = eph->tgd[0];
      if (sys == SYS_GAL)
      {
          pr_uura = (eph->ura - 2.0) * (obs_->psr_std[freq_idx]/0.16);
          dp_uura = (eph->ura - 2.0) * (obs_->dopp_std[freq_idx]/0.256);
      }
      else
      {
          pr_uura = (eph->ura - 1.0) * (obs_->psr_std[freq_idx]/0.16);
          dp_uura = (eph->ura - 1.0) * (obs_->dopp_std[freq_idx]/0.256);
      }
  }
  LOG_IF(FATAL, pr_uura <= 0) << "pr_uura is " << pr_uura; // get those parameters mainly, both used in eva
  LOG_IF(FATAL, dp_uura <= 0) << "dp_uura is " << dp_uura;
  // relative_sqrt_info = 10;
}

void GNSSProcess::SvPosCals(const ObsPtr &obs_, const EphemBasePtr &ephem_)
{
  // The primary pseudorange supplies transmit time. The resulting satellite
  // position is shared by all carrier bands at this epoch.
  freq = L1_freq(obs_, &freq_idx);
  LOG_IF(FATAL, freq < 0) << "No L1 observation found.";

  uint32_t sys = satsys(obs_->sat, NULL);
  double tof = obs_->psr[freq_idx] / LIGHT_SPEED;
  gtime_t sv_tx = time_add(obs_->time, -tof);

  if (sys == SYS_GLO)
  {
      GloEphemPtr glo_ephem = std::dynamic_pointer_cast<GloEphem>(ephem_);
      svdt = geph2svdt(sv_tx, glo_ephem);
      sv_tx = time_add(sv_tx, -svdt);
      sv_pos = geph2pos(sv_tx, glo_ephem, &svdt);
      sv_vel = geph2vel(sv_tx, glo_ephem, &svddt);
  }
  else
  {
      EphemPtr eph = std::dynamic_pointer_cast<Ephem>(ephem_);
      svdt = eph2svdt(sv_tx, eph); // used in eva
      sv_tx = time_add(sv_tx, -svdt);
      sv_pos = eph2pos(sv_tx, eph, &svdt); // used in eva
      sv_vel = eph2vel(sv_tx, eph, &svddt); // used in eva
  }
}

bool GNSSProcess::Evaluate(state_output &state)
{
  if (gnss_meas_buf[0].empty()) // ||gnss_meas_buf[0].size() < 4) //
  {
    // cout << "no valid gnss" << endl;
    return false;
  }

  double time_current = time2sec(gnss_meas_buf[0][0]->time);
  double delta_t = time_current - last_gnss_time;

  gtsam::Rot3 rel_rot; // = gtsam::Rot3(pre_integration->delta_q);
  gtsam::Point3 rel_pos, pos, ba, bg, acc, omg;
  gtsam::Vector3 rel_vel, vel;
  Eigen::Matrix3d rot = Eigen::Matrix3d::Identity();
  if (!nolidar) // && !invalid_lidar)
  {
    // LiDAR-enabled mode treats the external local state as this frame's
    // initial estimate and lets GNSS refine the local/global alignment.
    // Eigen::Matrix3d last_rot = p_assign->isamCurrentEstimate.at<gtsam::Rot3>(R(0)).matrix(); // state_const_.rot; //
    // // cout << "check time period" << pre_integration->sum_dt << ";" << time_current - last_gnss_time <<  endl;
    // Eigen::Vector3d last_pos = p_assign->isamCurrentEstimate.at<gtsam::Vector6>(A(0)).segment<3>(0); // state_.pos; //
    // Eigen::Vector3d last_vel = p_assign->isamCurrentEstimate.at<gtsam::Vector6>(A(0)).segment<3>(3); // state_.vel; //
    // Eigen::Vector3d cur_grav = state.rot.transpose() * state.gravity; //
    rot = state.rot; //.normalized().toRotationMatrix(); last_rot.transpose() *
    // rel_rot = gtsam::Rot3(last_rot.transpose() * state.rot.normalized().toRotationMatrix());
    pos = state.pos; // - anc_local; // last_rot.transpose() * (state.pos - last_pos - last_vel * delta_t - 0.5 * state.gravity * delta_t * delta_t);  - last_pos
    // (state.pos - last_pos);
    vel = state.vel; // last_rot.transpose() * (state.vel - last_vel - state.gravity * delta_t); // (state.vel - last_vel);  - last_vel
    ba = state.ba;
    bg = state.bg;
    acc = state.acc;
    omg = state.omg;
  }
  else
  {
    // No-LiDAR mode connects ECEF states using the accumulated IMU delta.
    ba = state.ba;
    bg = state.bg;
    rel_rot = gtsam::Rot3(pre_integration->delta_q);
    rel_pos = pre_integration->delta_p;
    rel_vel = pre_integration->delta_v;
  }

  initializeFrameEstimate(state);
  // rot_pos = state.rot; //.normalized().toRotationMatrix();
  if (AddFactor(rel_rot, rel_pos, rel_vel, state.gravity, delta_t, time_current, ba, bg, pos, vel, acc, omg, rot))
  {
    frame_num ++;
    // First update: obtain the float navigation and ambiguity solution.
    runISAM2opt();
    // Optional second update: validated integer priors immediately refine the
    // state returned from this same epoch.
    if (rtk_dd_added_this_epoch_)
      attemptIntegerAmbiguityResolution(time_current);
    // auto ekfPosNoise = p_assign->isam.marginalCovariance(A(frame_num-1));
    // odo_weight = 60 / (ekfPosNoise(0,0) + ekfPosNoise(1,1) + ekfPosNoise(2,2));
  }
  else
  {
    return false;
  }

  // state.cov.block<3,3>(0, 0) = isam.marginalCovariance(R(frame_num-1));
  // state.cov.block<6,6>(3, 3) = isam.marginalCovariance(F(frame_num-1)).block<6, 6>(0, 0);

  updateStateFromEstimate(state);
  last_gnss_time = time_current;
  pruneCarrierPhaseHistoryByTime(time_current);
  return true;
}

void GNSSProcess::initializeFrameEstimate(const state_output &state)
{
  if (!nolidar)
  {
    // Local mode separates position/velocity (A), motion/bias terms (O),
    // gravity (G), and attitude (R) for the corresponding factor types.
    Eigen::Matrix<double, 6, 1> position_velocity;
    position_velocity << state.pos, state.vel;
    Eigen::Matrix<double, 12, 1> motion_bias;
    motion_bias << state.omg, state.acc, state.bg, state.ba;
    p_assign->initialEstimate.insert(A(frame_num), gtsam::Vector6(position_velocity));
    p_assign->initialEstimate.insert(G(frame_num), gtsam::Vector3(state.gravity));
    p_assign->initialEstimate.insert(O(frame_num), gtsam::Vector12(motion_bias));
  }
  else
  {
    // Direct ECEF mode packs position, velocity, ba, and bg in F(k).
    Eigen::Matrix<double, 12, 1> navigation_state;
    navigation_state << state.pos, state.vel, state.ba, state.bg;
    p_assign->initialEstimate.insert(F(frame_num), gtsam::Vector12(navigation_state));
  }
  p_assign->initialEstimate.insert(R(frame_num), gtsam::Rot3(state.rot));
}

void GNSSProcess::updateStateFromEstimate(state_output &state)
{
  if (nolidar)
  {
    // In direct ECEF mode the optimized graph state is the exported state.
    const auto navigation_state =
        p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num - 1));
    state.rot = p_assign->isamCurrentEstimate.at<gtsam::Rot3>(R(frame_num - 1)).matrix();
    state.pos = navigation_state.segment<3>(0);
    state.vel = navigation_state.segment<3>(3);
    state.ba = navigation_state.segment<3>(6);
    state.bg = navigation_state.segment<3>(9);
    state.gravity = ecef2rotation(state.pos) * gravity_init;
    state.acc = state.rot.transpose() * (-state.gravity);
    return;
  }

  const auto position_velocity =
      p_assign->isamCurrentEstimate.at<gtsam::Vector6>(A(frame_num - 1));
  state_const_.rot = p_assign->isamCurrentEstimate.at<gtsam::Rot3>(R(frame_num - 1)).matrix();
  state_const_.pos = position_velocity.segment<3>(0);
  state_const_.vel = position_velocity.segment<3>(3);
  state.gravity = p_assign->isamCurrentEstimate.at<gtsam::Rot3>(P(0)).matrix().transpose() *
                  ecef2rotation(p_assign->isamCurrentEstimate.at<gtsam::Vector3>(E(0))) *
                  gravity_init;
}

void GNSSProcess::exportOptimizedState(state_output &state)
{
  updateStateFromEstimate(state);
}

void GNSSProcess::normalizeLidarInformation()
{
  if (nolidar)
  {
    return;
  }

  double lidar_weight = 0.0;
  if (p_assign->process_feat_num < 10)
  {
    invalid_lidar = true;
  }
  else
  {
    invalid_lidar = norm_vec_num < 10;
    lidar_weight = 2.0 * static_cast<double>(norm_vec_num) /
                   static_cast<double>(p_assign->process_feat_num);
  }
  norm_vec_num = 0;
  p_assign->process_feat_num = 0;

  const double minimum_information =
      std::min({sqrt_lidar(0, 0), sqrt_lidar(1, 1), sqrt_lidar(2, 2)});
  const double maximum_information =
      std::max({sqrt_lidar(0, 0), sqrt_lidar(1, 1), sqrt_lidar(2, 2)});
  if (minimum_information <= 0.0)
  {
    // Protect normalization from invalid/degenerate information matrices.
    ROS_WARN("Skipping lidar information normalization: non-positive diagonal");
    invalid_lidar = true;
    return;
  }

  const double scale = maximum_information / minimum_information > 3.0
                           ? 0.5
                           : std::max(lidar_weight, 0.5);
  sqrt_lidar *= scale / minimum_information;
}

void GNSSProcess::predictReceiverClock(double delta_t,
                                       double (&rcv_dt)[4],
                                       double &rcv_ddt) const
{
  // Constant-drift prediction: B_k = B_{k-1} + C_{k-1} * dt.
  rcv_ddt = p_assign->isamCurrentEstimate.at<gtsam::Vector1>(C(frame_num - 1))[0];
  const gtsam::Vector4 previous_bias =
      p_assign->isamCurrentEstimate.at<gtsam::Vector4>(B(frame_num - 1));
  for (int system = 0; system < 4; ++system)
  {
    rcv_dt[system] = previous_bias[system] + rcv_ddt * delta_t;
  }
}

void GNSSProcess::updateRoverSatelliteDirections(
    const std::vector<ObsPtr> &raw_observations,
    const std::vector<ObsPtr> &valid_observations,
    const std::vector<EphemBasePtr> &valid_ephemerides,
    const state_output &state)
{
  raw_satellite_directions_.clear();
  valid_satellite_directions_.clear();

  Eigen::Vector3d anchor_ecef = anc_ecef;
  Eigen::Matrix3d ecef_from_local = R_ecef_enu;
  if (p_assign->isamCurrentEstimate.exists(E(0)) &&
      p_assign->isamCurrentEstimate.exists(P(0)))
  {
    anchor_ecef =
        p_assign->isamCurrentEstimate.at<gtsam::Vector3>(E(0));
    ecef_from_local =
        p_assign->isamCurrentEstimate.at<gtsam::Rot3>(P(0)).matrix();
  }

  const Eigen::Vector3d receiver_local =
      state.pos + state.rot * Tex_imu_r;
  const Eigen::Vector3d receiver_ecef =
      anchor_ecef + ecef_from_local * receiver_local;

  const auto append_direction =
      [&](const ObsPtr &observation, const EphemBasePtr &ephemeris,
          std::set<uint32_t> &seen,
          std::vector<RtkSatelliteDirection> &output)
      {
        if (!observation || !seen.insert(observation->sat).second) return;
        const Eigen::Vector3d satellite_ecef =
            visualizationSatelliteEcef(observation, ephemeris);
        if (satellite_ecef.squaredNorm() <= 0.0) return;
        const Eigen::Vector3d local_direction = rtkLocalLineOfSight(
            satellite_ecef, receiver_ecef, ecef_from_local);
        if (local_direction.squaredNorm() <= 0.0) return;
        output.push_back({observation->sat, local_direction});
      };

  std::set<uint32_t> raw_satellites;
  for (const ObsPtr &observation : raw_observations)
  {
    append_direction(observation,
                     closestVisualizationEphemeris(*p_assign, observation),
                     raw_satellites, raw_satellite_directions_);
  }

  std::set<uint32_t> valid_satellites;
  const size_t valid_count =
      std::min(valid_observations.size(), valid_ephemerides.size());
  for (size_t index = 0; index < valid_count; ++index)
  {
    append_direction(valid_observations[index], valid_ephemerides[index],
                     valid_satellites, valid_satellite_directions_);
  }
}

// --------------------------------------------------------------------------
// Synchronized base/rover multi-frequency double-difference carrier processing
// --------------------------------------------------------------------------

void GNSSProcess::addDoubleDifferenceFactors(
    const std::vector<ObsPtr> &observations,
    const std::vector<EphemBasePtr> &ephemerides,
    const Eigen::Matrix3d &rover_rotation,
    double timestamp,
    std::vector<size_t> &factor_ids)
{
  rtk_dd_added_this_epoch_ = false;
  // Missing or unsynchronized base data disable only this optional RTK path;
  // pseudorange, Doppler, motion, and rover-only carrier factors still run.
  if (!use_double_differences || base_station_data_.empty()) return;
  double matched_base_timestamp = 0.0;
  const auto *base_epoch = base_station_data_.epoch(
      timestamp, base_epoch_tolerance, &matched_base_timestamp);
  if (!base_epoch)
  {
    ++rtk_base_epoch_misses;
    if (rtk_debug && isRtkDiagnosticEpoch(rtk_base_epoch_misses,
                                          rtk_debug_epoch_interval))
    {
      const double nearest_delta = matched_base_timestamp > 0.0
          ? matched_base_timestamp - timestamp
          : std::numeric_limits<double>::quiet_NaN();
      const bool outside_coverage =
          timestamp < base_station_data_.firstTimestamp() ||
          timestamp > base_station_data_.lastTimestamp();
      // std::ostringstream diagnostic;
      // diagnostic << "[RTK-TIME] no synchronized base epoch rover_gpst="
      //            << std::fixed << std::setprecision(3) << timestamp
      //            << " nearest_base_gpst=" << matched_base_timestamp
      //            << " delta_base_minus_rover_s=" << nearest_delta
      //            << " tolerance_s=" << base_epoch_tolerance
      //            << " outside_coverage=" << (outside_coverage ? "yes" : "no")
      //            << " cumulative_misses=" << rtk_base_epoch_misses
      //            << "; DD/LAMBDA skipped for this rover epoch";
      // With a 10 Hz rover and 1 Hz RINEX, most rover epochs are expected not
      // to match. Reserve WARN for a gap large enough to break an ambiguity
      // arc or for a true coverage problem.
      // if (!std::isfinite(nearest_delta) ||
      //     std::abs(nearest_delta) > rtk_ambiguity_gap_tolerance)
      //   ROS_WARN_STREAM(diagnostic.str());
      // else
      //   ROS_INFO_STREAM(diagnostic.str());
    }
    return;
  }
  ++rtk_synchronized_epochs;

  struct Candidate
  {
    size_t rover_index;
    int frequency_index;
    const BaseCarrierObservation *base;
    Eigen::Vector3d satellite_ecef;
    uint32_t system;
    RtkSignalBand signal_band;
    double wavelength;
    double single_difference_variance_m2;
    double single_difference_phase_m;
    double frequency_hz;
    bool rover_lock_lost;
    bool base_lock_lost;
  };

  std::vector<Candidate> candidates;
  size_t unsupported_signals = 0;
  size_t invalid_carriers = 0;
  size_t missing_base_carriers = 0;
  size_t half_cycle_rejections = 0;
  size_t rover_primary = 0, rover_secondary = 0, rover_l5 = 0;
  size_t matched_primary = 0, matched_secondary = 0, matched_l5 = 0;
  size_t rover_extra = 0, rover_wide = 0;
  size_t matched_extra = 0, matched_wide = 0;
  std::vector<std::string> unsupported_examples;
  std::vector<std::string> invalid_phase_examples;
  std::vector<std::string> no_base_examples;
  std::vector<std::string> half_cycle_examples;
  for (size_t index = 0; index < observations.size(); ++index)
  {
    const ObsPtr &observation = observations[index];
    SvPosCals(observations[index], ephemerides[index]);
    for (size_t frequency_index = 0;
         frequency_index < observation->freqs.size(); ++frequency_index)
    {
      // Classify every rover carrier rather than calling L1_freq() once. This
      // preserves separate primary, secondary, and L5 observations.
      const double frequency = observation->freqs[frequency_index];
      const double normalized_frequency =
          normalizeRtkFrequency(observation->sat, frequency);
      const RtkSignalBand signal_band =
          classifyRtkSignal(observation, frequency, use_l5, use_secondary);
      if (signal_band == RtkSignalBand::Unsupported || frequency <= 0.0)
      {
        ++unsupported_signals;
        if (unsupported_examples.size() < 8)
        {
          const uint32_t system = satsys(observation->sat, nullptr);
          std::ostringstream sample;
          sample << sat2str(observation->sat)
                 << '/' << rtkSystemName(system)
                 << "/freqMHz=" << std::fixed << std::setprecision(3)
                 << frequency / 1.0e6;
          if (frequency <= 0.0) sample << "/invalid_freq";
          unsupported_examples.push_back(sample.str());
        }
        continue;
      }
      if (signal_band == RtkSignalBand::Primary) ++rover_primary;
      else if (signal_band == RtkSignalBand::Secondary) ++rover_secondary;
      else if (signal_band == RtkSignalBand::L5) ++rover_l5;
      else if (signal_band == RtkSignalBand::Extra) ++rover_extra;
      else if (signal_band == RtkSignalBand::Wide) ++rover_wide;
      if (frequency_index >= observation->cp.size() ||
          frequency_index >= observation->cp_std.size() ||
          observation->cp[frequency_index] == 0.0)
      {
        ++invalid_carriers;
        if (invalid_phase_examples.size() < 8)
        {
          std::ostringstream sample;
          sample << sat2str(observation->sat) << '/' << rtkBandName(signal_band)
                 << "/freqMHz=" << std::fixed << std::setprecision(3)
                 << frequency / 1.0e6;
          if (frequency_index >= observation->cp.size())
            sample << "/missing_cp";
          else if (frequency_index >= observation->cp_std.size())
            sample << "/missing_cp_std";
          else
            sample << "/cp_zero";
          invalid_phase_examples.push_back(sample.str());
        }
        continue;
      }
      const BaseCarrierObservation *base = findMatchingBaseCarrier(
          *base_epoch, observation->sat, normalized_frequency);
      if (!base)
      {
        ++missing_base_carriers;
        if (no_base_examples.size() < 8)
        {
          std::ostringstream sample;
          sample << sat2str(observation->sat) << '/' << rtkBandName(signal_band)
                 << "/freqMHz=" << std::fixed << std::setprecision(3)
                 << frequency / 1.0e6;
          no_base_examples.push_back(sample.str());
        }
        continue;
      }
      const uint8_t rover_lli = frequency_index < observation->LLI.size()
          ? observation->LLI[frequency_index] : 0U;
      if (rtkHalfCycleInvalid(rover_lli) ||
          rtkHalfCycleInvalid(base->loss_of_lock))
      {
        ++half_cycle_rejections;
        if (half_cycle_examples.size() < 8)
        {
          std::ostringstream sample;
          sample << sat2str(observation->sat) << '/' << rtkBandName(signal_band)
                 << "/freqMHz=" << std::fixed << std::setprecision(3)
                 << frequency / 1.0e6
                 << "/rover_lli=" << static_cast<int>(rover_lli)
                 << "/base_lli=" << static_cast<int>(base->loss_of_lock);
          half_cycle_examples.push_back(sample.str());
        }
        continue;
      }
      if (signal_band == RtkSignalBand::Primary) ++matched_primary;
      else if (signal_band == RtkSignalBand::Secondary) ++matched_secondary;
      else if (signal_band == RtkSignalBand::L5) ++matched_l5;
      else if (signal_band == RtkSignalBand::Extra) ++matched_extra;
      else if (signal_band == RtkSignalBand::Wide) ++matched_wide;

      const double wavelength = LIGHT_SPEED / normalized_frequency;
      const double base_wavelength = LIGHT_SPEED / base->frequency_hz;
      const double single_difference_variance_m2 =
          rtkSingleDifferenceVarianceMeters2(
              observation->cp_std[frequency_index], wavelength,
              base->carrier_std_cycles, base_wavelength);
      if (!std::isfinite(single_difference_variance_m2))
      {
        ++invalid_carriers;
        continue;
      }
      candidates.push_back(
          {index, static_cast<int>(frequency_index), base, sv_pos,
           satsys(observation->sat, nullptr), signal_band, wavelength,
           single_difference_variance_m2,
           observation->cp[frequency_index] * wavelength -
               base->carrier_cycles * base_wavelength,
           normalized_frequency, rtkLossOfLock(rover_lli),
           rtkLossOfLock(base->loss_of_lock)});
    }
  }

  // Supplement the raw-frequency factors with primary/secondary wide-lane
  // and narrow-lane combinations.  Both combinations retain one geometry
  // term, so they use the same DD factor with independent ambiguity states.
  const size_t raw_candidate_count = candidates.size();
  size_t lane_candidates = 0;
  for (size_t primary_index = 0; primary_index < raw_candidate_count;
       ++primary_index)
  {
    const Candidate primary = candidates[primary_index];
    if (primary.signal_band != RtkSignalBand::Primary) continue;
    for (size_t secondary_index = 0; secondary_index < raw_candidate_count;
         ++secondary_index)
    {
      const Candidate secondary = candidates[secondary_index];
      if (secondary.rover_index != primary.rover_index ||
          secondary.signal_band != RtkSignalBand::Secondary)
        continue;
      for (const bool wide_lane : {true, false})
      {
        const RtkLaneCombination combination = rtkLaneCombination(
            primary.single_difference_phase_m,
            primary.single_difference_variance_m2, primary.frequency_hz,
            secondary.single_difference_phase_m,
            secondary.single_difference_variance_m2, secondary.frequency_hz,
            wide_lane);
        if (!(combination.wavelength_m > 0.0) ||
            !std::isfinite(combination.phase_m) ||
            !std::isfinite(combination.variance_m2))
          continue;
        candidates.push_back(
            {primary.rover_index, primary.frequency_index, primary.base,
             primary.satellite_ecef, primary.system,
             wide_lane ? RtkSignalBand::WideLane
                       : RtkSignalBand::NarrowLane,
             combination.wavelength_m, combination.variance_m2,
             combination.phase_m, 0.0,
             primary.rover_lock_lost || secondary.rover_lock_lost,
             primary.base_lock_lost || secondary.base_lock_lost});
        ++lane_candidates;
      }
      break;
    }
  }

  Eigen::Vector3d rover_ecef;
  Eigen::Vector3d antenna_offset;
  if (nolidar)
  {
    rover_ecef = p_assign->initialEstimate.at<gtsam::Vector12>(F(frame_num)).head<3>() +
                 rover_rotation * Tex_imu_r;
    antenna_offset = Tex_imu_r;
  }
  else
  {
    const Eigen::Vector3d local =
        p_assign->initialEstimate.at<gtsam::Vector6>(A(frame_num)).head<3>() +
        rover_rotation * Tex_imu_r;
    rover_ecef = p_assign->isamCurrentEstimate.at<gtsam::Vector3>(E(0)) +
                 p_assign->isamCurrentEstimate.at<gtsam::Rot3>(P(0)).matrix() * local;
    antenna_offset = rover_rotation * Tex_imu_r;

    // Publish one direction per satellite, even when several frequencies from
    // that satellite contribute DD factors.  The local direction preserves
    // azimuth/elevation while avoiding an unusable orbital-scale RViz scene.
    rtk_satellite_directions_.clear();
    const Eigen::Matrix3d ecef_from_local =
        p_assign->isamCurrentEstimate.at<gtsam::Rot3>(P(0)).matrix();
    std::set<uint32_t> visualized_satellites;
    for (const Candidate &candidate : candidates)
    {
      const uint32_t satellite = observations[candidate.rover_index]->sat;
      if (!visualized_satellites.insert(satellite).second) continue;
      const Eigen::Vector3d local_direction = rtkLocalLineOfSight(
          candidate.satellite_ecef, rover_ecef, ecef_from_local);
      if (local_direction.squaredNorm() <= 0.0) continue;
      rtk_satellite_directions_.push_back({satellite, local_direction});
    }
  }

  using ReferenceGroup = std::pair<uint32_t, RtkSignalBand>;
  // A reference is chosen independently per constellation and band. Mixing
  // frequencies in one DD would not leave a single integer-wavelength ambiguity.
  std::map<ReferenceGroup, size_t> reference_by_group;
  // Preserve the previous reference while it is still observed. This avoids
  // creating new DD ambiguity definitions because of small quality changes.
  for (size_t index = 0; index < candidates.size(); ++index)
  {
    const ReferenceGroup group(candidates[index].system,
                               candidates[index].signal_band);
    const auto found = reference_by_group.find(group);
    if (found == reference_by_group.end() ||
        candidates[index].single_difference_variance_m2 <
            candidates[found->second].single_difference_variance_m2)
      reference_by_group[group] = index;
  }
  for (size_t index = 0; index < candidates.size(); ++index)
  {
    const ReferenceGroup group(candidates[index].system,
                               candidates[index].signal_band);
    const auto previous_reference = reference_satellites_.find(
        {group.first, static_cast<uint8_t>(group.second)});
    if (previous_reference != reference_satellites_.end() &&
        previous_reference->second ==
            observations[candidates[index].rover_index]->sat)
      reference_by_group[group] = index;
  }
  for (const auto &entry : reference_by_group)
  {
    const Candidate &reference = candidates[entry.second];
    reference_satellites_[{entry.first.first,
                           static_cast<uint8_t>(entry.first.second)}] =
        observations[reference.rover_index]->sat;
  }

  // A code double difference cancels receiver and satellite clock offsets
  // without adding an ambiguity variable. Reuse the primary carrier reference
  // so the code and carrier DD residuals have identical geometry.
  size_t epoch_pseudorange_factors = 0;
  if (use_dd_pseudorange &&
      std::isfinite(double_difference_pseudorange_sigma) &&
      double_difference_pseudorange_sigma > 0.0)
  {
    const auto pseudorange_noise = gtsam::noiseModel::Isotropic::Sigma(
        1, double_difference_pseudorange_sigma);
    for (size_t index = 0; index < candidates.size(); ++index)
    {
      const Candidate &satellite = candidates[index];
      if (satellite.signal_band != RtkSignalBand::Primary) continue;
      const ReferenceGroup group(satellite.system, RtkSignalBand::Primary);
      const auto reference_entry = reference_by_group.find(group);
      if (reference_entry == reference_by_group.end() ||
          reference_entry->second == index)
        continue;
      const Candidate &reference = candidates[reference_entry->second];
      const ObsPtr &satellite_obs = observations[satellite.rover_index];
      const ObsPtr &reference_obs = observations[reference.rover_index];
      if (static_cast<size_t>(satellite.frequency_index) >= satellite_obs->psr.size() ||
          static_cast<size_t>(reference.frequency_index) >= reference_obs->psr.size() ||
          satellite_obs->psr[satellite.frequency_index] <= 0.0 ||
          reference_obs->psr[reference.frequency_index] <= 0.0 ||
          satellite.base->pseudorange_m <= 0.0 ||
          reference.base->pseudorange_m <= 0.0)
        continue;
      const double measured =
          (satellite_obs->psr[satellite.frequency_index] -
           satellite.base->pseudorange_m) -
          (reference_obs->psr[reference.frequency_index] -
           reference.base->pseudorange_m);
      if (nolidar)
        p_assign->gtSAMgraph.add(ligo::DoubleDiffPseudorangeFactorNolidar(
            R(frame_num), F(frame_num), Tex_imu_r,
            base_station_data_.ecefPosition(), satellite.satellite_ecef,
            reference.satellite_ecef, measured, pseudorange_noise));
      else
        p_assign->gtSAMgraph.add(ligo::DoubleDiffPseudorangeFactor(
            E(0), P(0), A(frame_num), antenna_offset,
            base_station_data_.ecefPosition(), satellite.satellite_ecef,
            reference.satellite_ecef, measured, pseudorange_noise));
      factor_ids.push_back(id_accumulate++);
      ++epoch_pseudorange_factors;
      ++rtk_double_difference_pseudorange_factors;
    }
  }

  const size_t factors_before = rtk_double_difference_factors;
  const size_t secondary_before = rtk_secondary_factors;
  const size_t l5_before = rtk_l5_factors;
  size_t new_arcs = 0;
  size_t reset_arcs = 0;
  double minimum_dd_sigma_m = std::numeric_limits<double>::infinity();
  double maximum_dd_sigma_m = 0.0;
  double sum_dd_sigma_m = 0.0;
  std::vector<std::string> dd_sigma_examples;

  const auto ambiguity_prior_noise =
      gtsam::noiseModel::Isotropic::Sigma(1, ambiguity_prior_sigma);
  if (use_dd_carrier)
  for (size_t index = 0; index < candidates.size(); ++index)
  {
    const Candidate &satellite = candidates[index];
    const ReferenceGroup group(satellite.system, satellite.signal_band);
    const Candidate &reference = candidates[reference_by_group.at(group)];
    if (&satellite == &reference) continue;

    const double measurement_sigma = rtkDoubleDifferenceSigmaMeters(
        satellite.single_difference_variance_m2,
        reference.single_difference_variance_m2,
        double_difference_sigma_floor);
    if (!std::isfinite(measurement_sigma) || measurement_sigma <= 0.0)
      continue;
    const auto measurement_noise =
        gtsam::noiseModel::Isotropic::Sigma(1, measurement_sigma);

    const ObsPtr &satellite_obs = observations[satellite.rover_index];
    const ObsPtr &reference_obs = observations[reference.rover_index];
    const double measured = satellite.single_difference_phase_m -
                            reference.single_difference_phase_m;
    // DD geometry cancels receiver/satellite clock terms. The ambiguity state
    // absorbs the remaining integer carrier term (and float atmospheric error).
    const double geometry =
        (satellite.satellite_ecef - rover_ecef).norm() -
        (satellite.satellite_ecef - base_station_data_.ecefPosition()).norm() -
        (reference.satellite_ecef - rover_ecef).norm() +
        (reference.satellite_ecef - base_station_data_.ecefPosition()).norm();

    const AmbiguityId ambiguity_id(
        reference_obs->sat, satellite_obs->sat,
        static_cast<uint8_t>(satellite.signal_band));
    auto ambiguity = ambiguities_.find(ambiguity_id);
    const bool ambiguity_existed = ambiguity != ambiguities_.end();
    const bool satellite_lock_lost = satellite.rover_lock_lost;
    const bool reference_lock_lost = reference.rover_lock_lost;
    const bool satellite_base_lock_lost = satellite.base_lock_lost;
    const bool reference_base_lock_lost = reference.base_lock_lost;
    const double arc_gap = ambiguity_existed
        ? timestamp - ambiguity->second.last_timestamp
        : std::numeric_limits<double>::quiet_NaN();
    const bool discontinuous_gap = ambiguity_existed &&
        !rtkArcIsContinuous(ambiguity->second.last_timestamp, timestamp,
                            rtk_ambiguity_gap_tolerance);
    const bool reset_ambiguity =
        ambiguity == ambiguities_.end() || satellite_lock_lost ||
        reference_lock_lost ||
        satellite_base_lock_lost || reference_base_lock_lost ||
        discontinuous_gap;
    if (reset_ambiguity)
    {
      // A slip, gap, reference change, or band change starts a new ambiguity
      // key so an earlier integer constraint can never contaminate this arc.
      const gtsam::Key key = gtsam::Symbol('n', next_ambiguity_id_++);
      const bool integer_compatible =
          std::abs(satellite.wavelength - reference.wavelength) < 1.0e-8;
      ambiguity = ambiguities_.insert_or_assign(
          ambiguity_id,
          AmbiguityState{key, timestamp,
                         integer_compatible ? satellite.wavelength : 0.0,
                         1, false, 0}).first;
      const double initial_ambiguity = measured - geometry;
      p_assign->initialEstimate.insert(key, gtsam::Vector1(initial_ambiguity));
      p_assign->gtSAMgraph.add(gtsam::PriorFactor<gtsam::Vector1>(
          key, gtsam::Vector1(initial_ambiguity), ambiguity_prior_noise));
      factor_ids.push_back(id_accumulate++);
      ++new_arcs;
      if (ambiguity_existed) ++reset_arcs;

      if (rtk_debug)
      {
        ROS_INFO_STREAM(
            "[RTK-DD] new float ambiguity epoch=" << std::fixed
            << std::setprecision(3) << timestamp
            << " key=" << gtsam::DefaultKeyFormatter(key)
            << " band=" << rtkBandName(satellite.signal_band)
            << " ref=" << sat2str(reference_obs->sat)
            << " sat=" << sat2str(satellite_obs->sat)
            << " measured_m=" << measured
            << " geometry_m=" << geometry
            << " dd_sigma_m=" << measurement_sigma
            << " target_sd_sigma_m="
            << std::sqrt(satellite.single_difference_variance_m2)
            << " ref_sd_sigma_m="
            << std::sqrt(reference.single_difference_variance_m2)
            << " sigma_floor_m=" << double_difference_sigma_floor
            << " float_init_m=" << initial_ambiguity
            << " float_init_cycles="
            << (ambiguity->second.wavelength > 0.0
                    ? initial_ambiguity / ambiguity->second.wavelength
                    : std::numeric_limits<double>::quiet_NaN())
            << " reset_reason="
            << (!ambiguity_existed ? "new" :
                satellite_lock_lost ? "rover_sat_lli" :
                reference_lock_lost ? "rover_ref_lli" :
                satellite_base_lock_lost ? "base_sat_lli" :
                reference_base_lock_lost ? "base_ref_lli" :
                discontinuous_gap ? "epoch_gap" : "unknown")
            << " rover_sat_lli="
            << (satellite.frequency_index < static_cast<int>(satellite_obs->LLI.size())
                    ? static_cast<int>(satellite_obs->LLI[satellite.frequency_index]) : -1)
            << " rover_ref_lli="
            << (reference.frequency_index < static_cast<int>(reference_obs->LLI.size())
                    ? static_cast<int>(reference_obs->LLI[reference.frequency_index]) : -1)
            << " base_sat_lli=" << static_cast<int>(satellite.base->loss_of_lock)
            << " base_ref_lli=" << static_cast<int>(reference.base->loss_of_lock)
            << " epoch_gap_s=" << arc_gap);
      }
    }
    else
    {
      // Consecutive epochs are later used as a lock-duration fixing gate.
      ++ambiguity->second.consecutive_epochs;
    }
    ambiguity->second.last_timestamp = timestamp;
    ++ambiguity->second.dd_observation_count;
    ambiguity->second.latest_measured = measured;
    ambiguity->second.latest_geometry = geometry;

    gtsam::NoiseModelFactor::shared_ptr carrier_factor;
    if (nolidar)
      carrier_factor =
          std::make_shared<ligo::DoubleDiffCarrierFactorNolidar>(
              R(frame_num), F(frame_num), ambiguity->second.key, Tex_imu_r,
              base_station_data_.ecefPosition(), satellite.satellite_ecef,
              reference.satellite_ecef, measured, measurement_noise);
    else
      carrier_factor = std::make_shared<ligo::DoubleDiffCarrierFactor>(
          E(0), P(0), A(frame_num), ambiguity->second.key, antenna_offset,
          base_station_data_.ecefPosition(), satellite.satellite_ecef,
          reference.satellite_ecef, measured, measurement_noise);
    p_assign->gtSAMgraph.add(carrier_factor);
    // Retain the exact current factor so the post-float residual gate uses
    // the optimized navigation state, not the pre-optimization geometry.
    ambiguity->second.latest_sigma_m = measurement_sigma;
    ambiguity->second.latest_factor = carrier_factor;
    factor_ids.push_back(id_accumulate++);
    ++rtk_double_difference_factors;
    minimum_dd_sigma_m = std::min(minimum_dd_sigma_m, measurement_sigma);
    maximum_dd_sigma_m = std::max(maximum_dd_sigma_m, measurement_sigma);
    sum_dd_sigma_m += measurement_sigma;
    if (dd_sigma_examples.size() < 8)
    {
      std::ostringstream sample;
      sample << sat2str(reference_obs->sat) << "->"
             << sat2str(satellite_obs->sat) << '/'
             << rtkBandName(satellite.signal_band)
             << "/sigma_m=" << std::fixed << std::setprecision(6)
             << measurement_sigma
             << "/target_sd_sigma_m="
             << std::sqrt(satellite.single_difference_variance_m2)
             << "/ref_sd_sigma_m="
             << std::sqrt(reference.single_difference_variance_m2)
             << "/floor_m=" << double_difference_sigma_floor;
      dd_sigma_examples.push_back(sample.str());
    }
    if (satellite.signal_band == RtkSignalBand::Secondary)
      ++rtk_secondary_factors;
    if (satellite.signal_band == RtkSignalBand::L5) ++rtk_l5_factors;
  }
  const size_t epoch_factors = rtk_double_difference_factors - factors_before;
  rtk_dd_added_this_epoch_ = epoch_factors > 0;
  if (epoch_factors == 0) ++rtk_zero_factor_epochs;
  if (rtk_debug &&
      (epoch_factors == 0 || isRtkDiagnosticEpoch(
                                 rtk_synchronized_epochs,
                                 rtk_debug_epoch_interval)))
  {
    std::ostringstream references;
    bool first = true;
    for (const auto &entry : reference_by_group)
    {
      if (!first) references << ',';
      first = false;
      references << sat2str(candidates[entry.second].base->satellite)
                 << '/' << rtkBandName(entry.first.second);
    }
    std::ostringstream unsupported_detail;
    for (size_t index = 0; index < unsupported_examples.size(); ++index)
    {
      if (index > 0) unsupported_detail << ',';
      unsupported_detail << unsupported_examples[index];
    }
    std::ostringstream invalid_phase_detail;
    for (size_t index = 0; index < invalid_phase_examples.size(); ++index)
    {
      if (index > 0) invalid_phase_detail << ',';
      invalid_phase_detail << invalid_phase_examples[index];
    }
    std::ostringstream no_base_detail;
    for (size_t index = 0; index < no_base_examples.size(); ++index)
    {
      if (index > 0) no_base_detail << ',';
      no_base_detail << no_base_examples[index];
    }
    std::ostringstream half_cycle_detail;
    for (size_t index = 0; index < half_cycle_examples.size(); ++index)
    {
      if (index > 0) half_cycle_detail << ',';
      half_cycle_detail << half_cycle_examples[index];
    }
    std::ostringstream dd_sigma_detail;
    for (size_t index = 0; index < dd_sigma_examples.size(); ++index)
    {
      if (index > 0) dd_sigma_detail << ',';
      dd_sigma_detail << dd_sigma_examples[index];
    }
    ROS_INFO_STREAM(
        "[RTK-DD] epoch=" << std::fixed << std::setprecision(3) << timestamp
        << " base_epoch=" << matched_base_timestamp
        << " base_minus_rover_ms="
        << (matched_base_timestamp - timestamp) * 1000.0
        << " synchronized_epoch=" << rtk_synchronized_epochs
        << " candidates=" << candidates.size()
        << " rover_signals=[primary:" << rover_primary
        << ",secondary:" << rover_secondary << ",l5:" << rover_l5
        << ",extra:" << rover_extra << ",wide:" << rover_wide << ']'
        << " matched_signals=[primary:" << matched_primary
        << ",secondary:" << matched_secondary << ",l5:" << matched_l5
        << ",extra:" << matched_extra << ",wide:" << matched_wide << ']'
        << " rejected=[unsupported:" << unsupported_signals
        << ",invalid_phase:" << invalid_carriers
        << ",no_base_signal:" << missing_base_carriers
        << ",half_cycle:" << half_cycle_rejections << ']'
        << " unsupported_examples=[" << unsupported_detail.str() << ']'
        << " invalid_phase_examples=[" << invalid_phase_detail.str() << ']'
        << " no_base_examples=[" << no_base_detail.str() << ']'
        << " half_cycle_examples=[" << half_cycle_detail.str() << ']'
        << " groups=" << reference_by_group.size()
        << " lane_candidates=" << lane_candidates
        << " references=[" << references.str() << ']'
        << " factors_added=" << epoch_factors
        << " pseudorange_factors_added=" << epoch_pseudorange_factors
        << " dd_sigma_m=[min:"
        << (epoch_factors > 0 ? minimum_dd_sigma_m : 0.0)
        << ",mean:"
        << (epoch_factors > 0 ? sum_dd_sigma_m / epoch_factors : 0.0)
        << ",max:"
        << (epoch_factors > 0 ? maximum_dd_sigma_m : 0.0)
        << ",floor:" << double_difference_sigma_floor << ']'
        << " dd_sigma_examples=[" << dd_sigma_detail.str() << ']'
        << " secondary_added="
        << (rtk_secondary_factors - secondary_before)
        << " l5_added=" << (rtk_l5_factors - l5_before)
        << " new_arcs=" << new_arcs
        << " reset_arcs=" << reset_arcs
        << " active_arcs=" << ambiguities_.size()
        << " zero_factor_epochs=" << rtk_zero_factor_epochs
        << " total_DD=" << rtk_double_difference_factors);
  }
}

// --------------------------------------------------------------------------
// Float-to-fixed ambiguity resolution
// --------------------------------------------------------------------------

bool GNSSProcess::attemptIntegerAmbiguityResolution(double timestamp)
{
  if (!enable_integer_fixing || !use_double_differences || ambiguities_.empty())
  {
    integer_solution_available = false;
    rtk_float_ambiguity_count = 0;
    rtk_fixed_ambiguity_count = 0;
    rtk_fix_status = !enable_integer_fixing ? "DISABLED" :
                     !use_double_differences ? "DD_DISABLED" : "NO_AMBIGUITIES";
    return false;
  }

  integer_solution_available = false;
  struct EligibleAmbiguity
  {
    AmbiguityId id;
    AmbiguityState *state;
  };
  // Evaluate the exact DD factors at the freshly optimized float solution.
  // Rejection is satellite-wide: a bad raw or combined-frequency factor
  // removes every ambiguity belonging to that DD target from this LAMBDA run.
  std::set<uint32_t> residual_rejected_satellites;
  size_t residual_factors_tested = 0;
  double maximum_absolute_normalized_residual = 0.0;
  std::vector<std::string> residual_rejection_examples;
  for (const auto &entry : ambiguities_)
  {
    const AmbiguityState &ambiguity = entry.second;
    if (std::abs(ambiguity.last_timestamp - timestamp) > 1.0e-6 ||
        !ambiguity.latest_factor || ambiguity.latest_sigma_m <= 0.0)
      continue;
    try
    {
      const gtsam::Vector error = ambiguity.latest_factor->unwhitenedError(
          p_assign->isamCurrentEstimate);
      if (error.size() != 1) continue;
      const double normalized = rtkStandardizedResidual(
          error[0], ambiguity.latest_sigma_m);
      ++residual_factors_tested;
      maximum_absolute_normalized_residual = std::max(
          maximum_absolute_normalized_residual, std::abs(normalized));
      if (std::abs(normalized) > lambda_max_normalized_residual)
      {
        residual_rejected_satellites.insert(std::get<1>(entry.first));
        if (residual_rejection_examples.size() < 8)
        {
          std::ostringstream sample;
          sample << sat2str(std::get<0>(entry.first)) << "->"
                 << sat2str(std::get<1>(entry.first)) << '/'
                 << rtkBandName(static_cast<RtkSignalBand>(
                        std::get<2>(entry.first)))
                 << "/residual_m=" << std::fixed << std::setprecision(4)
                 << error[0]
                 << "/sigma_m=" << ambiguity.latest_sigma_m
                 << "/normalized=" << normalized;
          residual_rejection_examples.push_back(sample.str());
        }
      }
    }
    catch (const std::exception &exception)
    {
      // If a current DD factor cannot be evaluated, its target is not safe to
      // include in integer validation. The float graph solution is retained.
      residual_rejected_satellites.insert(std::get<1>(entry.first));
      if (residual_rejection_examples.size() < 8)
      {
        std::ostringstream sample;
        sample << sat2str(std::get<0>(entry.first)) << "->"
               << sat2str(std::get<1>(entry.first))
               << "/factor_unavailable=" << exception.what();
        residual_rejection_examples.push_back(sample.str());
      }
    }
  }
  if (rtk_debug &&
      (!residual_rejected_satellites.empty() ||
       isRtkDiagnosticEpoch(rtk_synchronized_epochs,
                            rtk_debug_epoch_interval)))
  {
    std::ostringstream rejected_satellites, examples;
    for (const uint32_t satellite : residual_rejected_satellites)
    {
      if (rejected_satellites.tellp() > 0) rejected_satellites << ',';
      rejected_satellites << sat2str(satellite);
    }
    for (size_t index = 0; index < residual_rejection_examples.size(); ++index)
    {
      if (index > 0) examples << ',';
      examples << residual_rejection_examples[index];
    }
    ROS_INFO_STREAM("[RTK-RESIDUAL] tested=" << residual_factors_tested
                    << " threshold=" << lambda_max_normalized_residual
                    << " max_abs_normalized="
                    << maximum_absolute_normalized_residual
                    << " rejected_satellites=["
                    << rejected_satellites.str() << ']'
                    << " examples=[" << examples.str() << ']');
  }
  std::vector<EligibleAmbiguity> eligible;
  size_t active_float = 0;
  size_t active_fixed = 0;
  size_t waiting_for_lock = 0;
  for (auto &entry : ambiguities_)
  {
    AmbiguityState &ambiguity = entry.second;
    const bool active = rtkArcIsContinuous(
        ambiguity.last_timestamp, timestamp, rtk_ambiguity_gap_tolerance);
    if (!active) continue;
    if (ambiguity.fixed)
    {
      ++active_fixed;
      integer_solution_available = true;
      continue;
    }
    ++active_float;
    if (ambiguity.consecutive_epochs < lambda_min_lock_epochs)
      ++waiting_for_lock;
    if (ambiguity.wavelength > 0.0 &&
        ambiguity.consecutive_epochs >= lambda_min_lock_epochs &&
        std::abs(ambiguity.last_timestamp - timestamp) <= 1.0e-6 &&
        residual_rejected_satellites.count(std::get<1>(entry.first)) == 0 &&
        p_assign->isamCurrentEstimate.exists(ambiguity.key))
      eligible.push_back({entry.first, &ambiguity});
  }
  rtk_float_ambiguity_count = active_float;
  rtk_fixed_ambiguity_count = active_fixed;

  const size_t minimum = static_cast<size_t>(std::max(2, lambda_min_ambiguities));
  if (eligible.size() < minimum)
  {
    rtk_fix_status = waiting_for_lock > 0 ? "FLOAT_WAITING_FOR_LOCK" :
                     integer_solution_available ? "FIXED_HOLD" :
                     "FLOAT_TOO_FEW_ELIGIBLE";
    if (rtk_debug && isRtkDiagnosticEpoch(rtk_synchronized_epochs,
                                          rtk_debug_epoch_interval))
      ROS_INFO_STREAM(
          "[RTK-FLOAT] status=" << rtk_fix_status
          << " active_float=" << active_float
          << " active_fixed=" << active_fixed
          << " waiting_for_lock=" << waiting_for_lock
          << " lambda_eligible=" << eligible.size()
          << " residual_rejected_satellites="
          << residual_rejected_satellites.size()
          << " required=" << minimum);
    return false;
  }

  ++rtk_fix_attempts;
  const bool verbose_attempt =
      rtk_debug && isRtkDiagnosticEpoch(rtk_fix_attempts,
                                        rtk_debug_epoch_interval);
  rtk_fix_status = "FLOAT_LAMBDA_TESTING";
  if (verbose_attempt)
    ROS_INFO_STREAM("[RTK-LAMBDA] attempt=" << rtk_fix_attempts
                    << " eligible=" << eligible.size()
                    << " minimum=" << minimum
                    << " ratio_threshold=" << lambda_ratio_threshold
                    << " max_std_cycles=" << lambda_max_std_cycles);

  // Partial ambiguity resolution: remove the least precise ambiguity and retry
  // rather than rejecting an otherwise strong subset.
  while (eligible.size() >= minimum)
  {
    try
    {
      gtsam::KeyVector keys;
      keys.reserve(eligible.size());
      for (const EligibleAmbiguity &ambiguity : eligible)
        keys.push_back(ambiguity.state->key);
      const gtsam::JointMarginal joint =
          p_assign->isam.jointMarginalCovariance(keys);

      // Ambiguities are stored in metres in the graph. LAMBDA requires cycles,
      // so scale both the mean and every joint-covariance block by wavelength.
      Eigen::VectorXd floating(eligible.size());
      Eigen::MatrixXd covariance(eligible.size(), eligible.size());
      size_t worst_index = 0;
      double worst_std = -1.0;
      for (size_t row = 0; row < eligible.size(); ++row)
      {
        floating[row] =
            p_assign->isamCurrentEstimate.at<gtsam::Vector1>(eligible[row].state->key)[0] /
            eligible[row].state->wavelength;
        for (size_t column = 0; column < eligible.size(); ++column)
          covariance(row, column) =
              joint.at(eligible[row].state->key,
                       eligible[column].state->key)(0, 0) / (eligible[row].state->wavelength * eligible[column].state->wavelength);
        const double standard_deviation =
          covariance(row, row) > 0.0 ? std::sqrt(covariance(row, row)) : std::numeric_limits<double>::infinity();
        if (standard_deviation > worst_std)
        {
          worst_std = standard_deviation;
          worst_index = row;
        }
        if (verbose_attempt)
          ROS_INFO_STREAM(
              "[RTK-FLOAT] key="
              << gtsam::DefaultKeyFormatter(eligible[row].state->key)
              << " ref=" << sat2str(std::get<0>(eligible[row].id))
              << " sat=" << sat2str(std::get<1>(eligible[row].id))
              << " band=" << rtkBandName(
                     static_cast<RtkSignalBand>(std::get<2>(eligible[row].id)))
              << " lock_epochs=" << eligible[row].state->consecutive_epochs
              << " dd_observations="
              << eligible[row].state->dd_observation_count
              << " float_cycles=" << std::fixed << std::setprecision(4)
              << floating[row]
              << " nearest_integer_distance_cycles="
              << std::abs(floating[row] - std::round(floating[row]))
              << " std_cycles=" << standard_deviation
              << " latest_prefit_residual_m="
              << (eligible[row].state->latest_geometry +
                  floating[row] * eligible[row].state->wavelength -
                  eligible[row].state->latest_measured));
      }

      double maximum_absolute_correlation = 0.0;
      for (size_t row = 0; row < eligible.size(); ++row)
        for (size_t column = row + 1; column < eligible.size(); ++column)
        {
          const double denominator =
              std::sqrt(std::max(covariance(row, row), 0.0) *
                        std::max(covariance(column, column), 0.0));
          if (denominator > 0.0)
            maximum_absolute_correlation = std::max(
                maximum_absolute_correlation,
                std::abs(covariance(row, column) / denominator));
        }
      if (verbose_attempt)
        ROS_INFO_STREAM("[RTK-FLOAT] subset=" << eligible.size()
                        << " max_abs_correlation="
                        << maximum_absolute_correlation
                        << " worst_std_cycles=" << worst_std);

      if (!covariance.allFinite() || worst_std > lambda_max_std_cycles)
      {
        // Partial ambiguity resolution: retry without the weakest float state.
        if (verbose_attempt)
          ROS_WARN_STREAM(
              "[RTK-LAMBDA] precision gate rejected key="
              << gtsam::DefaultKeyFormatter(eligible[worst_index].state->key)
              << " std_cycles=" << worst_std
              << " limit=" << lambda_max_std_cycles
              << " subset_before=" << eligible.size());
        rtk_fix_status = "FLOAT_REJECTED_PRECISION";
        eligible.erase(eligible.begin() + worst_index);
        continue;
      }

      const LambdaResult result =
          LambdaAmbiguityResolver::solve(floating, covariance, 2);
      if (!result.valid || result.squared_norms.size() < 2)
      {
        if (verbose_attempt)
          ROS_WARN_STREAM("[RTK-LAMBDA] solver rejected subset="
                          << eligible.size() << "; removing weakest key="
                          << gtsam::DefaultKeyFormatter(
                                 eligible[worst_index].state->key));
        rtk_fix_status = "FLOAT_REJECTED_SOLVER";
        eligible.erase(eligible.begin() + worst_index);
        continue;
      }

      last_lambda_ratio = result.squared_norms[1] /
                          std::max(result.squared_norms[0], 1.0e-12);
      if (!std::isfinite(last_lambda_ratio) ||
          last_lambda_ratio < lambda_ratio_threshold)
      {
        if (verbose_attempt)
          ROS_WARN_STREAM(
              "[RTK-LAMBDA] ratio gate rejected subset=" << eligible.size()
              << " best_norm=" << result.squared_norms[0]
              << " second_norm=" << result.squared_norms[1]
              << " ratio=" << last_lambda_ratio
              << " required=" << lambda_ratio_threshold
              << "; removing weakest key="
              << gtsam::DefaultKeyFormatter(
                     eligible[worst_index].state->key));
        rtk_fix_status = "FLOAT_REJECTED_RATIO";
        eligible.erase(eligible.begin() + worst_index);
        continue;
      }

      for (size_t index = 0; index < eligible.size(); ++index)
      {
        AmbiguityState &ambiguity = *eligible[index].state;
        ambiguity.fixed_integer =
            static_cast<long long>(std::llround(result.candidates(index, 0)));
        ambiguity.fixed = true;
        const double fixed_value =
            static_cast<double>(ambiguity.fixed_integer) * ambiguity.wavelength;
        const double sigma = std::max(
            fixed_ambiguity_sigma_cycles * ambiguity.wavelength, 1.0e-6);
        p_assign->gtSAMgraph.add(gtsam::PriorFactor<gtsam::Vector1>(
            ambiguity.key, gtsam::Vector1(fixed_value),
            gtsam::noiseModel::Isotropic::Sigma(1, sigma)));
        ROS_INFO_STREAM(
            "[RTK-LAMBDA] FIXED key=" << gtsam::DefaultKeyFormatter(ambiguity.key)
            << " ref=" << sat2str(std::get<0>(eligible[index].id))
            << " sat=" << sat2str(std::get<1>(eligible[index].id))
            << " band=" << rtkBandName(
                   static_cast<RtkSignalBand>(std::get<2>(eligible[index].id)))
            << " float_cycles=" << std::fixed << std::setprecision(4)
            << floating[index]
            << " integer_cycles=" << ambiguity.fixed_integer
            << " correction_cycles="
            << (static_cast<double>(ambiguity.fixed_integer) - floating[index])
            << " fixed_m=" << fixed_value);
      }
      // Fix-and-hold: commit integer priors and immediately recalculate the
      // navigation estimate before Evaluate() exports it.
      p_assign->isam.update(p_assign->gtSAMgraph, gtsam::Values());
      p_assign->gtSAMgraph.resize(0);
      p_assign->isam.update();
      p_assign->isamCurrentEstimate = p_assign->isam.calculateEstimate();
      integer_solution_available = true;
      ++rtk_fix_successes;
      rtk_fixed_ambiguity_count += eligible.size();
      rtk_float_ambiguity_count -= std::min(rtk_float_ambiguity_count,
                                            eligible.size());
      rtk_fix_status = "FIXED";
      ROS_INFO_STREAM("[RTK-LAMBDA] FIX SUCCESS attempt=" << rtk_fix_attempts
                      << " fixed=" << eligible.size()
                      << " ratio=" << last_lambda_ratio
                      << " best_norm=" << result.squared_norms[0]
                      << " second_norm=" << result.squared_norms[1]
                      << " successful_fixes=" << rtk_fix_successes);
      return true;
    }
    catch (const std::exception &exception)
    {
      ROS_WARN_STREAM_THROTTLE(
          5.0, "[RTK-LAMBDA] covariance unavailable; FLOAT retained: "
                   << exception.what());
      rtk_fix_status = "FLOAT_COVARIANCE_UNAVAILABLE";
      return false;
    }
  }
  if (verbose_attempt)
    ROS_WARN_STREAM("[RTK-LAMBDA] NO FIX: status=" << rtk_fix_status
                    << " attempt=" << rtk_fix_attempts
                    << " last_ratio=" << last_lambda_ratio
                    << " active_float=" << rtk_float_ambiguity_count
                    << " active_fixed=" << rtk_fixed_ambiguity_count);
  return false;
}

bool GNSSProcess::AddFactor(gtsam::Rot3 rel_rot, gtsam::Point3 rel_pos, gtsam::Vector3 rel_vel, Eigen::Vector3d state_gravity, double delta_t, double time_current,
                Eigen::Vector3d ba, Eigen::Vector3d bg, Eigen::Vector3d pos, Eigen::Vector3d vel, Eigen::Vector3d acc, Eigen::Vector3d omg, Eigen::Matrix3d rot)
{
  // This method only assembles the current graph update. Evaluate() owns the
  // subsequent iSAM2 float update and optional integer-fix update.
  double rcv_dt[4];
  bool rcv_sys[4];
  rcv_sys[0] = false; rcv_sys[1] = false; rcv_sys[2] = false; rcv_sys[3] = false;
  double rcv_ddt;
  invalid_lidar = false;
  normalizeLidarInformation();
  // sqrt_lidar.block<6, 6>(3, 3) *= 2;
  if (nolidar_cur && !nolidar) nolidar_cur = false;
  predictReceiverClock(delta_t, rcv_dt, rcv_ddt);

  const std::vector<ObsPtr> &curr_obs = gnss_meas_buf[0];
  const std::vector<EphemBasePtr> &curr_ephem = gnss_ephem_buf[0];

  // Stage 1: find an earlier continuous carrier observation for each current
  // satellite. These pairs feed the legacy time-differenced carrier factors.
  std::map<sat_first, std::map<uint32_t, double[6]> >::reverse_iterator it;

  std::deque<uint32_t> pair_sat_copy;
  // std::deque<uint32_t>().swap(pair_sat_copy);

  std::deque<double> meas_sats, meas_sats_final;
  // std::deque<double> meas_cov_sats; //, meas_cov_sats_final;
  // std::deque<double> meas_time_sats, meas_time_sats_final;
  std::deque<int> meas_index_sats, meas_index_sats_final;
  std::deque<Eigen::Vector3d> meas_RTex_sats, meas_RTex_sats_final;
  std::deque<Eigen::Vector3d> meas_svpos_sats, meas_svpos_sats_final;

  for (uint32_t j = 0; j < curr_obs.size(); j++) //   && j < 10
  {
    std::map<uint32_t, double[6]>::iterator it_old; // t_old_best,
    double meas;
    // double meas_cov;
    double meas_time;
    int meas_index;
    Eigen::Vector3d meas_svpos;
    Eigen::Vector3d RTex_sats;
    // Eigen::Vector3d best_svpos;
    // if (pair_sat.size() > 0)
    // if (curr_obs[j]->cp_std[freq_idx] < p_assign->gnss_cp_std_threshold)
    {
      // if (j == pair_sat.front())
      // if (curr_obs[j]->cp[freq_idx] > 10)
      {
        bool cp_found = false;
        for (it = sat2cp.rbegin(); it != sat2cp.rend(); it++)
        {
          it_old = it->second.find(curr_obs[j]->sat); // the same satellite
          // it_old_best = it->second.find(curr_obs[best_sat]->sat);
          if (it_old != it->second.end()) // && it_old_best != it->second.end())
          {
            if (it->first.timecur >= p_assign->sat_track_time[curr_obs[j]->sat])
            {
            // if ((time_current - it->first.timecur) / gnss_sample_period <= p_assign->sat_track_status[curr_obs[best_sat]->sat] - p_assign->gnss_track_num_threshold &&
            // if ((time_current - it->first.timecur) / gnss_sample_period <= p_assign->sat_track_status[curr_obs[j]->sat]) //- p_assign->gnss_track_num_threshold)
            if (time_current > p_assign->sat_track_time[curr_obs[j]->sat] && p_assign->sat_track_status[curr_obs[j]->sat] > 0) //- p_assign->gnss_track_num_threshold)
            {
              cp_found = true;
              meas = it_old->second[0]; // - it_old_best->second[1] + it_old->second[1]); it_old_best->second[0] -
              // meas_cov = (it_old->second[2] * it_old->second[2]); // it_old_best->second[2] * it_old_best->second[2] +
              meas_time = it->first.timecur;
              meas_index = it->first.frame_num;
              RTex_sats << it->first.RTex[0], it->first.RTex[1], it->first.RTex[2];
              meas_svpos << it_old->second[3], it_old->second[4], it_old->second[5];
              // best_svpos << it_old_best->second[3], it_old_best->second[4], it_old_best->second[5];
              break;
            }
            }
          }
        }

        if (cp_found)
        {
          meas_sats.push_back(meas);
          // meas_cov_sats.push_back(meas_cov);
          // meas_time_sats.push_back(meas_time);
          meas_index_sats.push_back(meas_index);
          meas_svpos_sats.push_back(meas_svpos);
          meas_RTex_sats.push_back(RTex_sats);
          // meas_svpos_best.push_back(best_svpos);
          pair_sat_copy.push_back(j);
        }
        // pair_sat.pop_front();
      }
    }
  }

  std::map<uint32_t, double[6]> curr_cp_map;
  std::vector<double> meas_cp;
  // std::vector<double> cov_cp;
  std::vector<Eigen::Vector3d> sv_pos_pair, sat_svpos;
  double cov_cp_best, meas_cp_best; //, esti_cp_best,
  // Eigen::Vector3d sv_pos_best;
  std::vector<size_t> factor_id_cur, sys_idx_cp;
  M3D omg_skew;
  omg_skew << SKEW_SYM_MATRX(omg);
  Eigen::Vector3d hat_omg_T = omg_skew * Tex_imu_r;
  for (uint32_t j = 0; j < curr_obs.size(); j++) //   && j < 10
  {
    // Stage 2: calculate primary-frequency satellite corrections, cache the
    // current carrier sample, and add one pseudorange/Doppler factor.
    bool balance = false;
    if (j > curr_obs.size() / 2)
    {
      balance = true;
    }
    const uint32_t sys = satsys(curr_obs[j]->sat, NULL);
    const uint32_t sys_idx = gnss_comm::sys2idx.at(sys);
    GnssPsrDoppMeas(curr_obs[j], curr_ephem[j]); //, latest_gnss_iono_params);
    freq = L1_freq(curr_obs[j], &freq_idx); // save
    const double wavelength = LIGHT_SPEED / freq; // save
    if (curr_obs[j]->cp_std[freq_idx] < p_assign->gnss_cp_std_threshold)
    {
      if (curr_obs[j]->cp[freq_idx] * wavelength > 100)
      {
        curr_cp_map[curr_obs[j]->sat][0] = curr_obs[j]->cp[freq_idx] * wavelength + svdt * LIGHT_SPEED - tgd * LIGHT_SPEED;
        // curr_cp_map[curr_obs[j]->sat][2] = curr_obs[j]->cp_std[freq_idx] * 0.004;
        curr_cp_map[curr_obs[j]->sat][3] = sv_pos[0];
        curr_cp_map[curr_obs[j]->sat][4] = sv_pos[1];
        curr_cp_map[curr_obs[j]->sat][5] = sv_pos[2];

        if (pair_sat_copy.size() > 0)
        {
          for (size_t k = 0; k < pair_sat_copy.size(); k++)
          {
            if (j == pair_sat_copy[k])
            {
              meas_cp.push_back(curr_obs[j]->cp[freq_idx] * wavelength + svdt * LIGHT_SPEED - tgd * LIGHT_SPEED);
              sys_idx_cp.push_back(sys_idx);
              meas_sats_final.push_back(meas_sats[k]);
              // cov_cp.push_back(curr_obs[j]->cp_std[freq_idx] * curr_obs[j]->cp_std[freq_idx] * 0.004 * 0.004);
              // meas_cov_sats_final.push_back(meas_cov_sats[k]);
              sv_pos_pair.push_back(sv_pos);
              meas_svpos_sats_final.push_back(meas_svpos_sats[k]);
              // meas_time_sats_final.push_back(meas_time_sats[k]);
              meas_index_sats_final.push_back(meas_index_sats[k]);
              meas_RTex_sats_final.push_back(meas_RTex_sats[k]);
              break;
              // pair_sat_copy.pop_front();
            }
          }
        }
      }
    }
    /////////////////////////////////
    double values[27];
    values[0] = Tex_imu_r[0]; values[1] = Tex_imu_r[1]; values[2] = Tex_imu_r[2]; //values[3] = anc_local[0]; values[4] = anc_local[1]; values[5] = anc_local[2];
    values[3] = sv_pos[0]; values[4] = sv_pos[1]; values[5] = sv_pos[2]; values[6] = sv_vel[0]; values[7] = sv_vel[1]; values[8] = sv_vel[2];
    values[9] = svdt; values[10] = tgd; values[11] = svddt; values[12] = pr_uura; values[13] = dp_uura; values[14] = relative_sqrt_info; // psr_weight_adjust;
    values[15] = p_assign->latest_gnss_iono_params[0]; values[16] = p_assign->latest_gnss_iono_params[1]; values[17] = p_assign->latest_gnss_iono_params[2]; values[18] = p_assign->latest_gnss_iono_params[3];
    values[19] = p_assign->latest_gnss_iono_params[4]; values[20] = p_assign->latest_gnss_iono_params[5]; values[21] = p_assign->latest_gnss_iono_params[6]; values[22] = p_assign->latest_gnss_iono_params[7];
    values[23] = time_current; values[24] = freq; values[25] = psr_meas_hatch_filter[j]; values[26] = curr_obs[j]->dopp[freq_idx]; //curr_obs[j]->psr[freq_idx];
    rcv_sys[sys_idx] = true;
    if (!nolidar)
    {
      // Local navigation variables plus the persistent ECEF alignment.
      Eigen::Vector3d RTex = rot * Tex_imu_r;
      values[0] = RTex[0]; values[1] = RTex[1]; values[2] = RTex[2];
      if (frame_num < delete_thred)
      {
        p_assign->gtSAMgraph.add(ligo::GnssPsrDoppFactorNoR(A(frame_num), B(frame_num), C(frame_num), E(0), P(0), balance, values, sys_idx, rot * hat_omg_T, p_assign->robustpsrdoppNoise_init));
      }
      else
      {
        p_assign->gtSAMgraph.add(ligo::GnssPsrDoppFactorNoR(A(frame_num), B(frame_num), C(frame_num), E(0), P(0), balance, values, sys_idx, rot * hat_omg_T, p_assign->robustpsrdoppNoise));
      }
      // p_assign->gtSAMgraph.add(ligo::GnssPsrDoppFactorNoR(A(frame_num), B(frame_num), C(frame_num), E(0), P(0), invalid_lidar, values, sys_idx, rot * hat_omg_T, p_assign->robustpsrdoppNoise));
      // p_assign->gtSAMgraph.add(ligo::GnssPsrDoppFactor(R(frame_num), A(frame_num), B(frame_num), C(frame_num), E(0), P(0), invalid_lidar, values, sys_idx, hat_omg_T, p_assign->robustpsrdoppNoise));
      // p_assign->gtSAMgraph.add(ligo::GnssPsrDoppFactorPos(A(frame_num), B(frame_num), C(frame_num), E(0), P(0), invalid_lidar, values, sys_idx, rot_pos, hat_omg_T, p_assign->robustpsrdoppNoise));
    }
    else
    {
      // Direct ECEF navigation variables when LiDAR is unavailable.
      if (frame_num < delete_thred)
      {
        p_assign->gtSAMgraph.add(ligo::GnssPsrDoppFactorNolidar(R(frame_num), F(frame_num), B(frame_num), C(frame_num), values, sys_idx, hat_omg_T, p_assign->robustpsrdoppNoise_init)); // not work
      }
      else
      {
        p_assign->gtSAMgraph.add(ligo::GnssPsrDoppFactorNolidar(R(frame_num), F(frame_num), B(frame_num), C(frame_num), values, sys_idx, hat_omg_T, p_assign->robustpsrdoppNoise)); // not work
      }
    }
    factor_id_cur.push_back(id_accumulate);
    id_accumulate += 1;
  }
  sat_first cur_key;
  cur_key.RTex = rot * Tex_imu_r;
  cur_key.timecur = time2sec(curr_obs[0]->time);
  cur_key.frame_num = frame_num;
  sat2cp[cur_key] = curr_cp_map;

  // Stage 3: add synchronized base/rover primary and L5 DD carrier factors.
  addDoubleDifferenceFactors(curr_obs, curr_ephem, rot, time_current, factor_id_cur);
  // if (frame_num < delete_thred)
  // {
  //   p_assign->gtSAMgraph.add(ligo::DdtSmoothFactor(C(frame_num-1), C(frame_num), p_assign->ddtNoise_init));
  //   // p_assign->gtSAMgraph.add(gtsam::PriorFactor<gtsam::Vector1>(C(frame_num), gtsam::Vector1(rcv_ddt), p_assign->ddtNoise));
  //   p_assign->gtSAMgraph.add(ligo::DtDdtFactor(B(frame_num-1), B(frame_num), C(frame_num-1), C(frame_num), rcv_sys, delta_t, p_assign->dtNoise_init)); // not work
  // }
  // else
  // {
    // Stage 4: connect constellation clock biases and common clock drift.
    p_assign->gtSAMgraph.add(ligo::DdtSmoothFactor(C(frame_num-1), C(frame_num), p_assign->ddtNoise));
    // p_assign->gtSAMgraph.add(gtsam::PriorFactor<gtsam::Vector1>(C(frame_num), gtsam::Vector1(rcv_ddt), p_assign->ddtNoise));

  p_assign->gtSAMgraph.add(ligo::DtDdtFactor(B(frame_num-1), B(frame_num), C(frame_num-1), C(frame_num), rcv_sys, delta_t, p_assign->dtNoise)); // not work
  // }
  {
  // if (frame_num > 1)
  {
    p_assign->factor_id_frame[frame_num-1-frame_delete].push_back(id_accumulate+1);
    p_assign->factor_id_frame[frame_num-1-frame_delete].push_back(id_accumulate);
  }
  // else
  // {
    // p_assign->factor_id_frame.push_back(std::vector<size_t>(id_accumulate, id_accumulate+1));
  // }
  }
  // else
  // {
    // p_assign->factor_id_frame[frame_num-1-frame_delete].push_back(id_accumulate+1);
    // p_assign->factor_id_frame[frame_num-1-frame_delete].push_back(id_accumulate);
  // }
  id_accumulate += 2;
  if (!nolidar)
  {
    // Stage 5: constrain the graph with the current LiDAR/IMU local state.
    bool no_weight = false;
    // if (frame_num < delete_thred)
    // {
      // p_assign->gtSAMgraph.add(ligo::GnssLioFactor(P(0), E(0), R(0), A(0), R(frame_num), A(frame_num), gravity_init, state_gravity, ba, bg, rot, sqrt_lidar, p_assign->odomaNoise)); //LioNoise)); // odomNoiseIMU));
    // }
    // else
    {
      // p_assign->gtSAMgraph.add(ligo::GnssLioHardFactorNoR(A(frame_num), ba, bg, sqrt_lidar, no_weight, p_assign->odomNoise)); //LioNoise)); // odomNoiseIMU));
      p_assign->gtSAMgraph.add(ligo::GnssLioFactor(P(0), E(0), R(frame_num), A(frame_num), O(frame_num), G(frame_num), gravity_init, state_gravity, pos, vel, rot, ba, bg, acc, omg, sqrt_lidar, p_assign->odomNoise)); //LioNoise)); // odomNoiseIMU));
    }
      // p_assign->gtSAMgraph.add(ligo::GnssLioHardFactor(R(frame_num), A(frame_num), ba, bg, rot, sqrt_lidar, no_weight, p_assign->odomNoise)); //LioNoise)); // odomNoiseIMU));
    factor_id_cur.push_back(id_accumulate);
    id_accumulate += 1;
  }
  else
  {
    // In no-LiDAR mode use IMU preintegration between consecutive ECEF states.
    p_assign->gtSAMgraph.add(ligo::GnssLioFactorNolidar(R(frame_num-1), F(frame_num-1), R(frame_num), F(frame_num), rel_rot, rel_pos, rel_vel,
                  state_gravity, delta_t, ba, bg, pre_integration, p_assign->odomNoiseIMU));
    p_assign->factor_id_frame[frame_num-1-frame_delete].push_back(id_accumulate);
    id_accumulate += 1;
  }
  p_assign->initialEstimate.insert(C(frame_num), gtsam::Vector1(rcv_ddt));
  p_assign->initialEstimate.insert(B(frame_num), gtsam::Vector4(rcv_dt[0], rcv_dt[1], rcv_dt[2], rcv_dt[3]));

  // Stage 6: add time-differenced rover carrier factors. A continuous carrier
  // ambiguity cancels between the stored and current observations.
  for (uint32_t j = 0; j < meas_index_sats_final.size(); j++)
  {
    double values[11];
    values[0] = Tex_imu_r[0]; values[1] = Tex_imu_r[1]; values[2] = Tex_imu_r[2]; // values[3] = anc_local[0]; values[4] = anc_local[1]; values[5] = anc_local[2];
    values[3] = meas_svpos_sats_final[j][0]; values[4] = meas_svpos_sats_final[j][1]; values[5] = meas_svpos_sats_final[j][2];
    // values[9] = meas_svpos_sats_final[j][0]; values[10] = meas_svpos_sats_final[j][1]; values[11] = meas_svpos_sats_final[j][2];
    // values[12] = sv_pos_best[0]; values[13] = sv_pos_best[1]; values[14] = sv_pos_best[2];
    values[6] = sv_pos_pair[j][0]; values[7] = sv_pos_pair[j][1]; values[8] = sv_pos_pair[j][2];
    values[9] = meas_cp[j] - meas_sats_final[j]; values[10] = cp_weight; //_adjust;
    // values[14] = rcv_dt[0] - p_assign->isamCurrentEstimate.at<gtsam::Vector4>(B(meas_index_sats_final[j]))[0] + dt_com;
    if (!nolidar)
    {
      Eigen::Vector3d RTex1 = rot * Tex_imu_r;
      values[0] = RTex1[0]; values[1] = RTex1[1]; values[2] = RTex1[2];
      if (frame_num < delete_thred)
      {
        p_assign->gtSAMgraph.add(ligo::GnssCpFactorNoR(E(0), P(0), A(meas_index_sats_final[j]), A(frame_num), B(meas_index_sats_final[j]), B(frame_num), sys_idx_cp[j], invalid_lidar, values, meas_RTex_sats_final[j], p_assign->robustcpNoise_init));
      }
      else
      {
        p_assign->gtSAMgraph.add(ligo::GnssCpFactorNoR(E(0), P(0), A(meas_index_sats_final[j]), A(frame_num), B(meas_index_sats_final[j]), B(frame_num), sys_idx_cp[j], invalid_lidar, values, meas_RTex_sats_final[j], p_assign->robustcpNoise));
      }
      // p_assign->gtSAMgraph.add(ligo::GnssCpFactorNoR(E(0), P(0), A(meas_index_sats_final[j]), A(frame_num), B(meas_index_sats_final[j]), B(frame_num), sys_idx_cp[j], invalid_lidar, values, meas_RTex_sats_final[j], p_assign->robustcpNoise));
    }
    else
    {
      if (frame_num < delete_thred)
      {
        p_assign->gtSAMgraph.add(ligo::GnssCpFactorNolidar(R(meas_index_sats_final[j]), F(meas_index_sats_final[j]), R(frame_num), F(frame_num), B(meas_index_sats_final[j]), B(frame_num), sys_idx_cp[j], values, p_assign->robustcpNoise_init)); // not work
      }
      else
      {// p_assign->gtSAMgraph.add(ligo::GnssCpFactorNolidar(R(meas_index_sats_final[j]), F(meas_index_sats_final[j]), R(frame_num), F(frame_num), sys_idx_cp[j], values, p_assign->robustcpNoise)); // not work
        p_assign->gtSAMgraph.add(ligo::GnssCpFactorNolidar(R(meas_index_sats_final[j]), F(meas_index_sats_final[j]), R(frame_num), F(frame_num), B(meas_index_sats_final[j]), B(frame_num), sys_idx_cp[j], values, p_assign->robustcpNoise)); // not work
      }// Eigen::Matrix3d rot_before = p_assign->isamCurrentEstimate.at<gtsam::Rot3>(R(meas_index_sats_final[j])).matrix();
      // p_assign->gtSAMgraph.add(ligo::GnssCpFactorNolidarPos(F(meas_index_sats_final[j]), F(frame_num), values, rot_before, rot_pos, p_assign->robustcpNoise)); // not work
    }
    // factor_id_cur.push_back(id_accumulate);
    p_assign->factor_id_frame[meas_index_sats_final[j]-frame_delete].push_back(id_accumulate);
    id_accumulate += 1;
  }

  {
    // Ownership lists allow runISAM2opt() to remove factors with old frames.
    p_assign->factor_id_frame.push_back(factor_id_cur);
    std::vector<size_t>().swap(factor_id_cur);
  }
  // if (meas_index_sats_final.size() < 4)
  // {
  //   runISAM2opt();
  //   frame_num ++;
  //   return false;
  // }
  return true;
}

// --------------------------------------------------------------------------
// Initial graph priors
// --------------------------------------------------------------------------

void GNSSProcess::SetInit()
{
  if (!nolidar)
  {
    // Local mode anchors the global transform and initializes local navigation,
    // gravity, IMU terms, and receiver-clock states.
    // Eigen::Matrix3d R_enu_local_;
    // R_enu_local_ = R_ecef_enu; // * Rot_gnss_init; // * Eigen::AngleAxisd(yaw_enu_local, Eigen::Vector3d::UnitZ())
    // prior factor
    Eigen::Matrix<double, 6, 1> init_vel_bias_vector;
    Eigen::Matrix<double, 12, 1> init_others_vector;
    // init_vel_bias_vector.block<3,1>(0,0) = Rot_gnss_init.transpose() * pos_window[wind_size];
    init_vel_bias_vector.block<3,1>(0,0) = Eigen::Vector3d::Zero();
    init_vel_bias_vector.block<3,1>(3,0) = Eigen::Vector3d::Zero(); // vel_window[wind_size];
    init_others_vector.block<3,1>(0,0) = Eigen::Vector3d::Zero(); // vel_window[wind_size];
    init_others_vector.block<3,1>(3,0) = Eigen::Vector3d::Zero(); // vel_window[wind_size];
    init_others_vector.block<3,1>(6,0) = Eigen::Vector3d::Zero(); // vel_window[wind_size];
    init_others_vector.block<3,1>(9,0) = Eigen::Vector3d::Zero(); // vel_window[wind_size];
    // dt[0] = para_rcv_dt[wind_size*4]; dt[1] = para_rcv_dt[wind_size*4+1], dt[2] = para_rcv_dt[wind_size*4+2], dt[3] = para_rcv_dt[wind_size*4+3];
    // ddt = para_rcv_ddt[wind_size];
    p_assign->initialEstimate.insert(R(0), gtsam::Rot3(Rot_gnss_init)); //.transpose() * rot_window[wind_size]));
    p_assign->initialEstimate.insert(G(0), gtsam::Vector3(gravity_init)); //.transpose() * rot_window[wind_size]));
    // p_assign->initialEstimate.insert(F(0), gtsam::Vector12(init_vel_bias_vector));
    p_assign->initialEstimate.insert(A(0), gtsam::Vector6(init_vel_bias_vector));
    p_assign->initialEstimate.insert(O(0), gtsam::Vector12(init_others_vector));
    // p_assign->initialEstimate.insert(B(0), gtsam::Vector4(para_rcv_dt[wind_size*4], para_rcv_dt[wind_size*4+1], para_rcv_dt[wind_size*4+2], para_rcv_dt[wind_size*4+3]));
    p_assign->initialEstimate.insert(B(0), gtsam::Vector4(para_rcv_dt[4*wind_size], para_rcv_dt[4*wind_size], para_rcv_dt[4*wind_size], para_rcv_dt[4*wind_size])); //(1429495.922912-134967.935, 1429510.167255-134967.935, 1429516.520987-134967.935, 1429082.399893-134967.935)); //
    // p_assign->initialEstimate.insert(C(0), gtsam::Vector1(para_rcv_ddt[wind_size]));
    p_assign->initialEstimate.insert(C(0), gtsam::Vector1(para_rcv_ddt[0])); //(163.119147)); //(161.874045)
    // p_assign->initialEstimate.insert(Y(0), gtsam::Vector1(yaw_enu_local));
    p_assign->initialEstimate.insert(E(0), gtsam::Vector3(anc_ecef[0], anc_ecef[1], anc_ecef[2]));
    // cout << anc_ecef.transpose() << endl;
    p_assign->initialEstimate.insert(P(0), gtsam::Rot3(R_ecef_enu));

    gtsam::PriorFactor<gtsam::Rot3> init_rot_ext(P(0), gtsam::Rot3(gtsam::Rot3(R_ecef_enu)), p_assign->priorextrotNoise);
    gtsam::PriorFactor<gtsam::Vector3> init_pos_ext(E(0), gtsam::Vector3(anc_ecef[0], anc_ecef[1], anc_ecef[2]), p_assign->priorextposNoise);
    // gtsam::PriorFactor<gtsam::Vector4> init_dt(B(0), gtsam::Vector4(para_rcv_dt[wind_size*4], para_rcv_dt[wind_size*4+1], para_rcv_dt[wind_size*4+2], para_rcv_dt[wind_size*4+3]), p_assign->priordtNoise);
    gtsam::PriorFactor<gtsam::Vector4> init_dt(B(0), gtsam::Vector4(para_rcv_dt[4*wind_size], para_rcv_dt[4*wind_size], para_rcv_dt[4*wind_size], para_rcv_dt[4*wind_size]), p_assign->priordtNoise);
    // gtsam::PriorFactor<gtsam::Vector1> init_ddt(C(0), gtsam::Vector1(para_rcv_ddt[wind_size]), p_assign->priorddtNoise);
    gtsam::PriorFactor<gtsam::Vector1> init_ddt(C(0), gtsam::Vector1(para_rcv_ddt[0]), p_assign->priorddtNoise); // (161.874045) 163.119147
    gtsam::PriorFactor<gtsam::Rot3> init_rot_(R(0), gtsam::Rot3(Rot_gnss_init), p_assign->priorrotNoise);
    gtsam::PriorFactor<gtsam::Vector6> init_vel_(A(0), gtsam::Vector6(init_vel_bias_vector), p_assign->priorNoise); // priorposNoise);
    gtsam::PriorFactor<gtsam::Vector12> init_bias_(O(0), gtsam::Vector12(init_others_vector), p_assign->priorBiasNoise); // priorposNoise);
    gtsam::PriorFactor<gtsam::Vector3> init_grav_(G(0), gtsam::Vector3(gravity_init), p_assign->priorGravNoise);
    p_assign->gtSAMgraph.add(init_rot_ext);
    p_assign->gtSAMgraph.add(init_pos_ext);
    p_assign->gtSAMgraph.add(init_dt);
    p_assign->gtSAMgraph.add(init_ddt);
    p_assign->gtSAMgraph.add(init_rot_);
    p_assign->gtSAMgraph.add(init_vel_);
    p_assign->gtSAMgraph.add(init_bias_);
    p_assign->gtSAMgraph.add(init_grav_);
    p_assign->factor_id_frame.push_back(std::vector<size_t>{0, 1, 2, 3, 4, 5, 6, 7});
    id_accumulate += 8;
  }
  else
  {
    // No-LiDAR mode initializes a direct ECEF navigation trajectory.
  //   Eigen::Matrix3d R_enu_local_;
  //   R_enu_local_ = Eigen::AngleAxisd(yaw_enu_local, Eigen::Vector3d::UnitZ());
    // dt[0] = para_rcv_dt[wind_size*4], dt[1] = para_rcv_dt[wind_size*4+1], dt[2] = para_rcv_dt[wind_size*4+2], dt[3] = para_rcv_dt[wind_size*4+3];
    // ddt = para_rcv_ddt[wind_size];
    gtsam::PriorFactor<gtsam::Rot3> init_rot(R(0), gtsam::Rot3(R_ecef_enu * rot_window[wind_size]), p_assign->priorrotNoise); //  * R_enu_local_
    Eigen::Matrix<double, 12, 1> init_vel_bias_vector;
    init_vel_bias_vector.block<3,1>(0,0) = anc_ecef + R_ecef_enu * (pos_window[wind_size] - rot_window[wind_size] * Tex_imu_r); //  * R_enu_local_- pos_window[0]
    init_vel_bias_vector.block<3,1>(3,0) = R_ecef_enu * vel_window[wind_size]; // R_enu_local_ *
    init_vel_bias_vector.block<6,1>(6,0) = Eigen::Matrix<double, 6, 1>::Zero();
    gtsam::PriorFactor<gtsam::Vector12> init_vel_bias(F(0), gtsam::Vector12(init_vel_bias_vector), p_assign->priorposNoise);
    // gtsam::PriorFactor<gtsam::Vector4> init_dt(B(0), gtsam::Vector4(para_rcv_dt[wind_size*4], para_rcv_dt[wind_size*4+1], para_rcv_dt[wind_size*4+2], para_rcv_dt[wind_size*4+3]), p_assign->priordtNoise);
    gtsam::PriorFactor<gtsam::Vector4> init_dt(B(0), gtsam::Vector4(para_rcv_dt[0], para_rcv_dt[0], para_rcv_dt[0], para_rcv_dt[0]), p_assign->priordtNoise);
    gtsam::PriorFactor<gtsam::Vector1> init_ddt(C(0), gtsam::Vector1(para_rcv_ddt[0]), p_assign->priorddtNoise); // para_rcv_ddt[wind_size]
    p_assign->gtSAMgraph.add(init_rot);
    p_assign->gtSAMgraph.add(init_vel_bias);
    p_assign->gtSAMgraph.add(init_dt);
    p_assign->gtSAMgraph.add(init_ddt);
    p_assign->factor_id_frame.push_back(std::vector<size_t>{0, 1, 2, 3}); //{i * 4, i * 4 + 1, i * 4  + 2, i * 4 + 3});
    p_assign->initialEstimate.insert(R(0), gtsam::Rot3(R_ecef_enu * rot_window[wind_size])); // R_enu_local_ *
    p_assign->initialEstimate.insert(F(0), gtsam::Vector12(init_vel_bias_vector));
    // p_assign->initialEstimate.insert(B(0), gtsam::Vector4(para_rcv_dt[wind_size*4], para_rcv_dt[wind_size*4+1], para_rcv_dt[wind_size*4+2], para_rcv_dt[wind_size*4+3]));
    p_assign->initialEstimate.insert(B(0), gtsam::Vector4(para_rcv_dt[0], para_rcv_dt[0], para_rcv_dt[0], para_rcv_dt[0]));
    p_assign->initialEstimate.insert(C(0), gtsam::Vector1(para_rcv_ddt[0])); // para_rcv_ddt[wind_size]
    id_accumulate += 4;
  }
}
