#pragma once

#include <Eigen/Core>

#include <cstdint>
#include <map>
#include <string>
#include <vector>

struct BaseCarrierObservation
{
  uint32_t satellite = 0;
  double frequency_hz = 0.0;
  double carrier_cycles = 0.0;
  double carrier_std_cycles = 0.0;
  uint8_t loss_of_lock = 0;
};

class BaseStationData
{
 public:
  /** Load a RINEX 3.x observation file, a ZIP containing observation files,
   * or a directory of either. RINEX 3.02 is the primary supported format.
   */
  bool load(const std::string &path, std::string &error,
            double default_carrier_std_cycles = 0.01);
  /** Return the nearest epoch when it is within tolerance.  matched_timestamp
   * is filled with the nearest available base timestamp even when the match
   * fails, which lets runtime diagnostics distinguish a normal 1 Hz sampling
   * miss from an out-of-coverage or time-system error.
   */
  const std::vector<BaseCarrierObservation> *epoch(double timestamp,
                                                    double tolerance,
                                                    double *matched_timestamp = nullptr) const;
  bool empty() const { return epochs_.empty(); }
  const Eigen::Vector3d &ecefPosition() const { return ecef_position_; }
  size_t epochCount() const { return epochs_.size(); }
  const std::string &timeSystem() const { return time_system_; }
  double firstTimestamp() const;
  double lastTimestamp() const;

 private:
  Eigen::Vector3d ecef_position_ = Eigen::Vector3d::Zero();
  std::map<double, std::vector<BaseCarrierObservation>> epochs_;
  std::string time_system_ = "UNKNOWN";
};
