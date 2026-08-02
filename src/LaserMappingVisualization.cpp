#include "LaserMappingVisualization.h"

#include "RtkSignalUtils.h"

#include <geometry_msgs/Point.h>
#include <std_msgs/ColorRGBA.h>
#include <visualization_msgs/Marker.h>

#include <string>

namespace
{
std_msgs::ColorRGBA stageColor(SatelliteMarkerStage stage, float alpha = 1.0f)
{
  std_msgs::ColorRGBA color;
  color.a = alpha;
  switch (stage)
  {
    case SatelliteMarkerStage::Raw:
      color.r = color.g = color.b = 0.55f;
      break;
    case SatelliteMarkerStage::Valid:
      color.r = 1.0f; color.g = 0.72f; color.b = 0.05f;
      break;
    case SatelliteMarkerStage::RtkCandidate:
      color.r = 0.1f; color.g = 0.95f; color.b = 0.25f;
      break;
  }
  return color;
}

const char *stageName(SatelliteMarkerStage stage)
{
  switch (stage)
  {
    case SatelliteMarkerStage::Raw: return "RAW";
    case SatelliteMarkerStage::Valid: return "VALID";
    case SatelliteMarkerStage::RtkCandidate: return "RTK";
  }
  return "UNKNOWN";
}

geometry_msgs::Point rosPoint(const Eigen::Vector3d &point)
{
  geometry_msgs::Point output;
  output.x = point.x(); output.y = point.y(); output.z = point.z();
  return output;
}
}  // namespace

std::map<uint32_t, StagedSatellite> stageRtkSatellites(
    const std::vector<GNSSProcess::RtkSatelliteDirection> &raw,
    const std::vector<GNSSProcess::RtkSatelliteDirection> &valid,
    const std::vector<GNSSProcess::RtkSatelliteDirection> &rtk)
{
  std::map<uint32_t, StagedSatellite> satellites;
  for (const auto &direction : raw)
    satellites[direction.satellite] = {direction, SatelliteMarkerStage::Raw};
  for (const auto &direction : valid)
    satellites[direction.satellite] = {direction, SatelliteMarkerStage::Valid};
  for (const auto &direction : rtk)
  {
    const auto found = satellites.find(direction.satellite);
    if (found == satellites.end() ||
        found->second.stage != SatelliteMarkerStage::Valid)
      continue;
    found->second = {direction, SatelliteMarkerStage::RtkCandidate};
  }
  return satellites;
}

void RtkSatelliteVisualizer::publish(
    const ros::Publisher &publisher, const GNSSProcess &gnss,
    const Eigen::Vector3d &receiver_local, double timestamp)
{
  const auto satellites = stageRtkSatellites(
      gnss.rawSatelliteDirections(), gnss.validSatelliteDirections(),
      gnss.rtkSatelliteDirections());
  if (satellites.empty() && !previously_published_) return;

  visualization_msgs::Marker clear;
  clear.header.frame_id = "camera_init";
  clear.header.stamp = ros::Time().fromSec(timestamp);
  clear.action = visualization_msgs::Marker::DELETEALL;
  publisher.publish(clear);
  previously_published_ = !satellites.empty();
  if (satellites.empty()) return;

  const auto make_points = [&](const std::string &name,
                               SatelliteMarkerStage stage) {
    visualization_msgs::Marker marker;
    marker.header = clear.header; marker.ns = name; marker.id = 0;
    marker.type = visualization_msgs::Marker::SPHERE_LIST;
    marker.action = visualization_msgs::Marker::ADD;
    marker.pose.orientation.w = 1.0;
    marker.scale.x = marker.scale.y = marker.scale.z = 1.8;
    marker.color = stageColor(stage);
    return marker;
  };
  const auto make_rays = [&](const std::string &name,
                             SatelliteMarkerStage stage) {
    visualization_msgs::Marker marker;
    marker.header = clear.header; marker.ns = name; marker.id = 0;
    marker.type = visualization_msgs::Marker::LINE_LIST;
    marker.action = visualization_msgs::Marker::ADD;
    marker.pose.orientation.w = 1.0; marker.scale.x = 0.14;
    marker.color = stageColor(stage, 0.45f);
    return marker;
  };

  visualization_msgs::Marker raw_points =
      make_points("satellites_raw_only", SatelliteMarkerStage::Raw);
  visualization_msgs::Marker valid_points =
      make_points("satellites_valid", SatelliteMarkerStage::Valid);
  visualization_msgs::Marker rtk_points =
      make_points("satellites_rtk_candidates", SatelliteMarkerStage::RtkCandidate);
  visualization_msgs::Marker raw_rays =
      make_rays("satellite_rays_raw_only", SatelliteMarkerStage::Raw);
  visualization_msgs::Marker valid_rays =
      make_rays("satellite_rays_valid", SatelliteMarkerStage::Valid);
  visualization_msgs::Marker rtk_rays =
      make_rays("satellite_rays_rtk_candidates", SatelliteMarkerStage::RtkCandidate);

  const geometry_msgs::Point receiver_point = rosPoint(receiver_local);
  size_t label_id = 100, raw_count = 0, valid_count = 0, rtk_count = 0;
  for (const auto &entry : satellites)
  {
    const auto &direction = entry.second.direction;
    const SatelliteMarkerStage stage = entry.second.stage;
    const geometry_msgs::Point equivalent_point = rosPoint(
        rtkEquivalentSatellitePoint(receiver_local,
                                    direction.local_unit_direction,
                                    gnss.rtk_satellite_display_radius));
    visualization_msgs::Marker *points = nullptr, *rays = nullptr;
    if (stage == SatelliteMarkerStage::Raw)
    { points = &raw_points; rays = &raw_rays; ++raw_count; }
    else if (stage == SatelliteMarkerStage::Valid)
    { points = &valid_points; rays = &valid_rays; ++valid_count; }
    else
    { points = &rtk_points; rays = &rtk_rays; ++rtk_count; }
    points->points.push_back(equivalent_point);
    rays->points.push_back(receiver_point);
    rays->points.push_back(equivalent_point);

    visualization_msgs::Marker label;
    label.header = clear.header; label.ns = "satellite_stage_labels";
    label.id = static_cast<int>(label_id++);
    label.type = visualization_msgs::Marker::TEXT_VIEW_FACING;
    label.action = visualization_msgs::Marker::ADD;
    label.pose.position = equivalent_point; label.pose.position.z += 1.8;
    label.pose.orientation.w = 1.0; label.scale.z = 1.8;
    label.color = stageColor(stage);
    label.text = gnss_comm::sat2str(direction.satellite) +
                 std::string(" [") + stageName(stage) + "]";
    publisher.publish(label);
  }

  if (!raw_points.points.empty())
  { publisher.publish(raw_points); publisher.publish(raw_rays); }
  if (!valid_points.points.empty())
  { publisher.publish(valid_points); publisher.publish(valid_rays); }
  if (!rtk_points.points.empty())
  { publisher.publish(rtk_points); publisher.publish(rtk_rays); }

  visualization_msgs::Marker legend;
  legend.header = clear.header; legend.ns = "satellite_stage_legend";
  legend.id = 0; legend.type = visualization_msgs::Marker::TEXT_VIEW_FACING;
  legend.action = visualization_msgs::Marker::ADD;
  legend.pose.position = receiver_point; legend.pose.position.z += 4.0;
  legend.pose.orientation.w = 1.0; legend.scale.z = 1.5;
  legend.color.r = legend.color.g = legend.color.b = 1.0f;
  legend.color.a = 0.95f;
  legend.text = "RAW-only(gray): " + std::to_string(raw_count) +
                "  VALID(amber): " + std::to_string(valid_count) +
                "  RTK(green): " + std::to_string(rtk_count);
  publisher.publish(legend);
}
