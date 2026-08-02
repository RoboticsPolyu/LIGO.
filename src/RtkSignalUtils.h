#pragma once

#include "BaseStationData.h"

#include <gnss_comm/gnss_constant.hpp>
#include <gnss_comm/gnss_utility.hpp>

#include <Eigen/Core>

#include <cmath>
#include <cstdint>
#include <vector>

enum class RtkSignalBand : uint8_t
{
  Unsupported = 0,
  Primary = 1,
  Secondary = 2,
  L5 = 5,
  Extra = 6,
  Wide = 8
};

constexpr double kRtkFrequencyToleranceHz = 1.0e5;

inline bool rtkFrequencyNear(double frequency_hz, double nominal_hz)
{
  return std::abs(frequency_hz - nominal_hz) <= kRtkFrequencyToleranceHz;
}

inline bool rtkGlonassFrequencyNear(double frequency_hz, double base_hz,
                                    double step_hz)
{
  if (frequency_hz <= 0.0) return false;
  const double channel = std::round((frequency_hz - base_hz) / step_hz);
  if (channel < -7.0 || channel > 6.0) return false;
  return std::abs(frequency_hz - (base_hz + channel * step_hz)) <=
         kRtkFrequencyToleranceHz;
}

// Some historical ublox_driver bags tag Galileo E5b observations with the
// GPS-L2 centre frequency (1227.60 MHz).  The carrier values in those records
// are Galileo E5b measurements, so normalize the metadata to the standard
// E5b frequency before base matching or wavelength conversion.  Keeping this
// compatibility rule constellation-specific avoids reinterpreting real GPS L2.
inline double normalizeRtkFrequency(uint32_t satellite, double frequency_hz)
{
  const uint32_t system = gnss_comm::satsys(satellite, nullptr);
  if (system == SYS_GAL && rtkFrequencyNear(frequency_hz, FREQ2))
    return FREQ7;
  return frequency_hz;
}

inline RtkSignalBand classifyRtkSignal(const gnss_comm::ObsPtr &observation,
                                       double frequency_hz, bool enable_l5,
                                       bool enable_secondary)
{
  const double normalized_frequency =
      normalizeRtkFrequency(observation->sat, frequency_hz);
  const double primary_frequency = gnss_comm::L1_freq(observation, nullptr);
  if (primary_frequency > 0.0 &&
      rtkFrequencyNear(normalized_frequency, primary_frequency))
    return RtkSignalBand::Primary;

  const uint32_t system = gnss_comm::satsys(observation->sat, nullptr);
  const bool is_secondary =
      (system == SYS_GPS && rtkFrequencyNear(normalized_frequency, FREQ2)) ||
      (system == SYS_GLO &&
       rtkGlonassFrequencyNear(normalized_frequency, FREQ2_GLO, DFRQ2_GLO)) ||
      (system == SYS_GAL && rtkFrequencyNear(normalized_frequency, FREQ7)) ||
      (system == SYS_BDS &&
       rtkFrequencyNear(normalized_frequency, FREQ2_BDS));
  if (enable_secondary && is_secondary) return RtkSignalBand::Secondary;

  const bool system_has_l5 =
      system == SYS_GPS || system == SYS_GAL || system == SYS_BDS;
  if (enable_l5 && system_has_l5 &&
      rtkFrequencyNear(normalized_frequency, FREQ5))
    return RtkSignalBand::L5;

  const bool is_extra =
      (system == SYS_GAL && rtkFrequencyNear(normalized_frequency, FREQ6)) ||
      (system == SYS_BDS &&
       rtkFrequencyNear(normalized_frequency, FREQ3_BDS)) ||
      (system == SYS_GLO &&
       rtkFrequencyNear(normalized_frequency, FREQ3_GLO));
  if (enable_secondary && is_extra) return RtkSignalBand::Extra;

  const bool is_wide =
      ((system == SYS_GAL || system == SYS_BDS) &&
       rtkFrequencyNear(normalized_frequency, FREQ8));
  if (enable_secondary && is_wide) return RtkSignalBand::Wide;

  return RtkSignalBand::Unsupported;
}

// RINEX 3.02 LLI is a bit mask. Bits 0 and 1 invalidate an integer arc;
// Galileo's bit 2 only identifies BOC tracking and must not cause a reset.
inline bool rtkLossOfLock(uint8_t lli)
{
  return (lli & 0x03U) != 0U;
}

inline bool rtkCycleSlip(uint8_t lli)
{
  return (lli & 0x01U) != 0U;
}

// Bit 1 means that the receiver cannot guarantee a full-cycle carrier phase.
// Such a sample must not define an integer ambiguity or be selected as a DD
// reference. In contrast, a bit-0 cycle-slip sample can start a new valid arc.
inline bool rtkHalfCycleInvalid(uint8_t lli)
{
  return (lli & 0x02U) != 0U;
}

inline bool rtkArcIsContinuous(double previous_timestamp,
                               double current_timestamp,
                               double gap_tolerance_seconds)
{
  return previous_timestamp > 0.0 && current_timestamp >= previous_timestamp &&
         current_timestamp - previous_timestamp <= gap_tolerance_seconds;
}

// Convert the receiver and base carrier-phase uncertainties to the variance
// of one rover-minus-base single difference. Inputs are standard deviations
// in cycles and wavelengths in metres/cycle.
inline double rtkSingleDifferenceVarianceMeters2(
    double rover_std_cycles, double rover_wavelength_m,
    double base_std_cycles, double base_wavelength_m)
{
  if (!std::isfinite(rover_std_cycles) || rover_std_cycles < 0.0 ||
      !std::isfinite(rover_wavelength_m) || rover_wavelength_m <= 0.0 ||
      !std::isfinite(base_std_cycles) || base_std_cycles < 0.0 ||
      !std::isfinite(base_wavelength_m) || base_wavelength_m <= 0.0)
    return std::numeric_limits<double>::quiet_NaN();
  const double rover_sigma_m = rover_std_cycles * rover_wavelength_m;
  const double base_sigma_m = base_std_cycles * base_wavelength_m;
  return rover_sigma_m * rover_sigma_m + base_sigma_m * base_sigma_m;
}

// A carrier double difference contains four measurements: rover/base for the
// target satellite and rover/base for the reference. Their independent
// variances add. The floor represents residual atmosphere, multipath and
// ephemeris modelling uncertainty; it is not the complete DD noise anymore.
inline double rtkDoubleDifferenceSigmaMeters(
    double satellite_single_difference_variance_m2,
    double reference_single_difference_variance_m2,
    double modelling_sigma_floor_m)
{
  if (!std::isfinite(satellite_single_difference_variance_m2) ||
      satellite_single_difference_variance_m2 < 0.0 ||
      !std::isfinite(reference_single_difference_variance_m2) ||
      reference_single_difference_variance_m2 < 0.0 ||
      !std::isfinite(modelling_sigma_floor_m) ||
      modelling_sigma_floor_m < 0.0)
    return std::numeric_limits<double>::quiet_NaN();
  return std::sqrt(satellite_single_difference_variance_m2 +
                   reference_single_difference_variance_m2 +
                   modelling_sigma_floor_m * modelling_sigma_floor_m);
}

// Convert an ECEF satellite position into a unit line-of-sight vector in the
// local mapping frame.  RViz cannot use the roughly 20,000 km ECEF position
// directly alongside a local LiDAR map, so only the physically meaningful
// direction is retained for visualization.
inline Eigen::Vector3d rtkLocalLineOfSight(
    const Eigen::Vector3d &satellite_ecef,
    const Eigen::Vector3d &receiver_ecef,
    const Eigen::Matrix3d &ecef_from_local)
{
  const Eigen::Vector3d ecef_line_of_sight = satellite_ecef - receiver_ecef;
  const double norm = ecef_line_of_sight.norm();
  if (!ecef_line_of_sight.allFinite() || !ecef_from_local.allFinite() ||
      !std::isfinite(norm) || norm <= 0.0)
    return Eigen::Vector3d::Zero();
  return ecef_from_local.transpose() * (ecef_line_of_sight / norm);
}

// Project a satellite direction onto a small sphere centred on the rover.
// The resulting point is an equivalent display position, not satellite ECEF.
inline Eigen::Vector3d rtkEquivalentSatellitePoint(
    const Eigen::Vector3d &receiver_local,
    const Eigen::Vector3d &local_line_of_sight,
    double display_radius)
{
  const double norm = local_line_of_sight.norm();
  if (!receiver_local.allFinite() || !local_line_of_sight.allFinite() ||
      !std::isfinite(norm) || !std::isfinite(display_radius) || norm <= 0.0 ||
      display_radius <= 0.0)
    return receiver_local;
  return receiver_local + display_radius * local_line_of_sight / norm;
}

inline const BaseCarrierObservation *findMatchingBaseCarrier(
    const std::vector<BaseCarrierObservation> &base_epoch, uint32_t satellite,
    double frequency_hz)
{
  const double normalized_frequency =
      normalizeRtkFrequency(satellite, frequency_hz);
  const BaseCarrierObservation *best = nullptr;
  double best_frequency_error = kRtkFrequencyToleranceHz;
  for (const BaseCarrierObservation &observation : base_epoch)
  {
    if (observation.satellite != satellite) continue;
    const double frequency_error =
        std::abs(observation.frequency_hz - normalized_frequency);
    if (frequency_error <= best_frequency_error)
    {
      best = &observation;
      best_frequency_error = frequency_error;
    }
  }
  return best;
}
