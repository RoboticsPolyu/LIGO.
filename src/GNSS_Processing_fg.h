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

#pragma once
#include "GNSS_Initialization.h"
#include "GNSS_Assignment.h"
#include "BaseStationData.h"

#include <string>
#include <tuple>

#include <pcl/registration/icp.h>
using namespace gnss_comm;

#define WINDOW_SIZE (10) // should be 0

/**
 * GNSS/LiDAR/IMU factor-graph front end.
 *
 * Lifecycle:
 *   1. processGNSS() filters rover observations and fills the initialization
 *      window until GNSSLIAlign() establishes the local-to-ECEF transform.
 *   2. Evaluate() creates one graph state, calls AddFactor(), runs iSAM2, and
 *      writes the optimized result back to state_output.
 *   3. If synchronized base data are enabled, multi-frequency double-difference
 *      ambiguities are estimated as floats and conditionally fixed by LAMBDA.
 *
 * Graph-key convention used by this class:
 *   A(k): local position/velocity       R(k): attitude
 *   O(k): motion and IMU biases         G(k): gravity
 *   F(k): no-LiDAR ECEF navigation      B(k): constellation clock biases
 *   C(k): receiver clock drift          E(0): ECEF anchor
 *   P(0): local-to-ECEF rotation        n(j): DD carrier ambiguity
 */
class GNSSProcess
{
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  struct RtkSatelliteDirection
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    uint32_t satellite = 0;
    Eigen::Vector3d local_unit_direction = Eigen::Vector3d::Zero();
  };

  GNSSProcess();
  ~GNSSProcess();
  
  /// Clear all buffers, graph state, ambiguity arcs, and initialization state.
  void Reset();
  /// Supply the eight broadcast ionosphere coefficients used by GNSS factors.
  void inputIonoParams(double ts, const std::vector<double> &iono_params);
  void inputGNSSTimeDiff(const double t_diff);
  /// Load surveyed base coordinates and frequency-tagged carrier observations.
  bool loadBaseStationFile(const std::string &path);
  /// Latest raw rover satellite directions that have usable broadcast ephemerides.
  const std::vector<RtkSatelliteDirection> &rawSatelliteDirections() const
  {
    return raw_satellite_directions_;
  }
  /// Latest front-end-valid rover satellite directions in the local map frame.
  const std::vector<RtkSatelliteDirection> &validSatelliteDirections() const
  {
    return valid_satellite_directions_;
  }
  /// Latest synchronized RTK-candidate directions in the local map frame.
  const std::vector<RtkSatelliteDirection> &rtkSatelliteDirections() const
  {
    return rtk_satellite_directions_;
  }
  /// Filter an incoming rover epoch and either initialize or buffer it for Evaluate().
  void processGNSS(const std::vector<ObsPtr> &gnss_meas, state_output &state);
  /// Estimate the local-to-ECEF alignment from the initialization window.
  bool GNSSLIAlign();
  void updateGNSSStatistics(Eigen::Vector3d &pos);
  void inputpvt(double ts, double lat, double lon, double alt, int float_sol, int diff_sol);
  void inputgt(double ts, double lat, double lon, double alt);
  void inputlla(double ts, double lat, double lon, double alt);
  void processIMUOutput(double dt, const Vector3d &linear_acceleration, const Vector3d &angular_velocity);
  void processIMU(double dt, const Vector3d &linear_acceleration, const Vector3d &angular_velocity);
  Eigen::Vector3d local2enu(Eigen::Matrix3d enu_rot, Eigen::Vector3d anc, Eigen::Vector3d &pos);
  /// Insert priors for the first graph frame after local/ECEF alignment.
  void SetInit();
  /// Assemble all measurement, clock, motion, carrier, and RTK factors for one epoch.
  bool AddFactor(gtsam::Rot3 rel_rot_, gtsam::Point3 rel_pos_, gtsam::Vector3 rel_v_, Eigen::Vector3d state_gravity, double delta_t, double time_current,
                Eigen::Vector3d ba, Eigen::Vector3d bg,  Eigen::Vector3d pos, Eigen::Vector3d vel, Eigen::Vector3d acc, Eigen::Vector3d omg, Eigen::Matrix3d rot);
  std::map<sat_first, std::map<uint32_t, double[6]>> sat2cp; // 
  std::vector<ObsPtr> gnss_meas_buf[WINDOW_SIZE+1]; //
  std::vector<double> psr_meas_hatch_filter;
  std::vector<EphemBasePtr> gnss_ephem_buf[WINDOW_SIZE+1]; // 
  Eigen::Matrix3d rot_window[WINDOW_SIZE+1]; //
  Eigen::Vector3d pos_window[WINDOW_SIZE+1]; //
  Eigen::Vector3d pos_ecef_window[WINDOW_SIZE+1]; //
  Eigen::Vector3d vel_window[WINDOW_SIZE+1]; //
  Eigen::Vector3d Tex_imu_r;
  std::vector<double> pvt_time;
  std::vector<double> lla_time;
  std::vector<Eigen::Vector3d> pvt_holder;
  std::vector<int> diff_holder;
  std::vector<int> float_holder;
  std::vector<Eigen::Vector3d> lla_holder;
  std::queue<std::vector<ObsPtr>> gnss_msg;

  bool invalid_lidar = false;
  // double dt[4];
  // double ddt; 
  size_t id_accumulate = 0; // 
  size_t frame_delete = 0; // 

  int frame_num = 0; // 
  double last_gnss_time = 0.0; //
  double first_gnss_time = 0.0; //
  double gnss_sample_period = 0.1;

  double diff_t_gnss_local = 0.0;
  // double gnss_cp_std_threshold = 30;
  double gnss_cp_time_threshold = 30;
  Eigen::Vector3d ecef_pos, first_xyz_ecef_pvt, first_xyz_ecef_lla, first_lla_pvt, first_lla_lla;
  Eigen::Matrix3d Rot_gnss_init = Eigen::Matrix3d::Identity();
  bool gnss_ready = false;
  int frame_count = 0; //
  int delete_thred = 0;
  int wind_size = WINDOW_SIZE;
  int norm_vec_num = 0;
  bool nolidar = false;
  bool nolidar_cur = false;
  std::vector<Eigen::Vector3d> norm_vec_holder;
  // double para_yaw_enu_local[1];
  double para_rcv_dt[(WINDOW_SIZE+1)*4] = {0}; //
  double para_rcv_ddt[WINDOW_SIZE+1] = {0}; //
  Eigen::Vector3d anc_ecef, gravity_init;
  Eigen::Vector3d anc_local = Eigen::Vector3d::Zero();
  Eigen::Matrix3d R_ecef_enu;
  double yaw_enu_local = 0.0;
  
  /// Commit pending factors, enforce the fixed-lag window, and calculate an estimate.
  void runISAM2opt(void);
  /// Compute primary-frequency pseudorange/Doppler satellite terms and weights.
  void GnssPsrDoppMeas(const ObsPtr &obs_, const EphemBasePtr &ephem_);
  /// Compute satellite position/velocity for carrier-factor geometry.
  void SvPosCals(const ObsPtr &obs_, const EphemBasePtr &ephem_);
  /// Process one initialized GNSS epoch and export the optimized state.
  bool Evaluate(state_output &state);
  state_output state_const_;
  state_output state_const_last;
  double relative_sqrt_info = 10;
  double cp_weight = 1.0;
  // Base/rover RTK controls. The legacy pseudorange/Doppler path remains on
  // the primary signal; use_secondary/use_l5 apply to the DD carrier path.
  bool use_double_differences = false;
  bool use_secondary = true;
  bool use_l5 = true;
  double base_epoch_tolerance = 0.05;
  double rtk_ambiguity_gap_tolerance = 1.5;
  double base_carrier_std_cycles = 0.01;
  // Residual modelling floor only. Each DD factor also propagates the four
  // rover/base carrier-phase standard deviations from its observations.
  double double_difference_sigma_floor = 0.003;
  double double_difference_pseudorange_sigma = 2.0;
  double ambiguity_prior_sigma = 100.0;
  // Conservative LAMBDA validation and fix-and-hold thresholds.
  bool enable_integer_fixing = true;
  int lambda_min_ambiguities = 4;
  int lambda_min_lock_epochs = 10;
  double lambda_ratio_threshold = 3.0;
  double lambda_max_std_cycles = 0.25;
  double fixed_ambiguity_sigma_cycles = 0.001;
  bool integer_solution_available = false;
  double last_lambda_ratio = 0.0;
  size_t rtk_synchronized_epochs = 0;
  size_t rtk_base_epoch_misses = 0;
  size_t rtk_zero_factor_epochs = 0;
  size_t rtk_double_difference_factors = 0;
  size_t rtk_double_difference_pseudorange_factors = 0;
  size_t rtk_secondary_factors = 0;
  size_t rtk_l5_factors = 0;
  // Runtime observability for distinguishing synchronized, float, and fixed RTK.
  bool rtk_debug = true;
  int rtk_debug_epoch_interval = 10;
  bool rtk_satellite_visualization = true;
  double rtk_satellite_display_radius = 50.0;
  size_t rtk_fix_attempts = 0;
  size_t rtk_fix_successes = 0;
  size_t rtk_float_ambiguity_count = 0;
  size_t rtk_fixed_ambiguity_count = 0;
  std::string rtk_fix_status = "NOT_ATTEMPTED";
  Eigen::Matrix<double, 24, 24> sqrt_lidar;
  double odo_weight1 = 1.0;
  double odo_weight2 = 1.0;
  double odo_weight3 = 1.0;
  // Eigen::Matrix<double, 3, 3> rot_weight = Eigen::Matrix3d::Identity();
  double odo_weight4 = 2.0;
  double odo_weight5 = 2.0;
  double odo_weight6 = 2.0;
  IntegrationBase* pre_integration = new IntegrationBase{Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero()};
  GNSSAssignment* p_assign = new GNSSAssignment();
  private:
    // Small lifecycle helpers extracted from the original monolithic loop.
    void clearGnssBuffer(size_t index);
    void shiftInitializationWindow(size_t first_index);
    int firstValidClockBias(const Eigen::Matrix<double, 7, 1> &rough_xyzt) const;
    void initializeFrameEstimate(const state_output &state);
    void updateStateFromEstimate(state_output &state);
    void pruneCarrierPhaseHistoryByFrame();
    void pruneCarrierPhaseHistoryByTime(double time_current);
    void normalizeLidarInformation();
    void predictReceiverClock(double delta_t, double (&rcv_dt)[4], double &rcv_ddt) const;
    /// Refresh raw and front-end-valid sky-view layers for the current rover epoch.
    void updateRoverSatelliteDirections(
        const std::vector<ObsPtr> &raw_observations,
        const std::vector<ObsPtr> &valid_observations,
        const std::vector<EphemBasePtr> &valid_ephemerides,
        const state_output &state);
    /// Add raw-frequency and supplementary wide-/narrow-lane carrier DD factors.
    void addDoubleDifferenceFactors(const std::vector<ObsPtr> &observations,
                                    const std::vector<EphemBasePtr> &ephemerides,
                                    const Eigen::Matrix3d &rover_rotation,
                                    double timestamp,
                                    std::vector<size_t> &factor_ids);
    /// Resolve eligible float ambiguities and re-optimize after an accepted fix.
    bool attemptIntegerAmbiguityResolution(double timestamp);

    /// State belonging to one continuous reference/subject/frequency ambiguity arc.
    struct AmbiguityState
    {
      gtsam::Key key = 0;
      double last_timestamp = 0.0;
      double wavelength = 0.0;
      int consecutive_epochs = 1;
      bool fixed = false;
      long long fixed_integer = 0;
      size_t dd_observation_count = 0;
      double latest_measured = 0.0;
      double latest_geometry = 0.0;
    };
    BaseStationData base_station_data_;
    // (reference satellite, subject satellite, RtkSignalBand value).
    using AmbiguityId = std::tuple<uint32_t, uint32_t, uint8_t>;
    std::map<AmbiguityId, AmbiguityState> ambiguities_;
    // Keep one reference satellite per constellation/band until it disappears.
    std::map<std::pair<uint32_t, uint8_t>, uint32_t> reference_satellites_;
    size_t next_ambiguity_id_ = 0;
    bool rtk_dd_added_this_epoch_ = false;
    std::vector<RtkSatelliteDirection> raw_satellite_directions_;
    std::vector<RtkSatelliteDirection> valid_satellite_directions_;
    std::vector<RtkSatelliteDirection> rtk_satellite_directions_;

    const ObsPtr obs;
    const EphemBasePtr ephem;
    const std::vector<double> iono_paras;
    int freq_idx = 0;
    double freq = 0.0;
    Eigen::Vector3d sv_pos;
    Eigen::Vector3d sv_vel;
    double svdt, svddt, tgd;
    double pr_uura, dp_uura;
    Eigen::Matrix3d rot_pos;
};

// # endif
