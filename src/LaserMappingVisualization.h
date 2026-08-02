#pragma once

#include "GNSS_Processing_fg.h"

#include <Eigen/Core>
#include <ros/publisher.h>

#include <cstdint>
#include <map>
#include <vector>

enum class SatelliteMarkerStage
{
  Raw,
  Valid,
  RtkCandidate
};

struct StagedSatellite
{
  GNSSProcess::RtkSatelliteDirection direction;
  SatelliteMarkerStage stage = SatelliteMarkerStage::Raw;
};

// Merge the three processing layers. RTK status upgrades only a satellite
// that is valid in the current rover epoch, preventing stale RTK markers.
std::map<uint32_t, StagedSatellite> stageRtkSatellites(
    const std::vector<GNSSProcess::RtkSatelliteDirection> &raw,
    const std::vector<GNSSProcess::RtkSatelliteDirection> &valid,
    const std::vector<GNSSProcess::RtkSatelliteDirection> &rtk);

class RtkSatelliteVisualizer
{
 public:
  void publish(const ros::Publisher &publisher, const GNSSProcess &gnss,
               const Eigen::Vector3d &receiver_local, double timestamp);

 private:
  bool previously_published_ = false;
};
