#pragma once

#include <ros/ros.h>

class LaserMappingRosInterface
{
 public:
  explicit LaserMappingRosInterface(ros::NodeHandle &node);

  ros::Subscriber lidar;
  ros::Subscriber imu;
  ros::Subscriber ephemeris;
  ros::Subscriber glonass_ephemeris;
  ros::Subscriber gnss_measurement;
  ros::Subscriber ionosphere;
  ros::Subscriber gnss_time_pulse;
  ros::Subscriber local_trigger;
  ros::Subscriber rtk_pvt;
  ros::Subscriber rtk_lla;

  ros::Publisher registered_cloud;
  ros::Publisher registered_body_cloud;
  ros::Publisher effective_cloud;
  ros::Publisher laser_map;
  ros::Publisher mapped_odometry;
  ros::Publisher path;
  ros::Publisher plane_marker;
  ros::Publisher rtk_satellites;

 private:
  void configureSubscribers(ros::NodeHandle &node);
  void configurePublishers(ros::NodeHandle &node);
  void loadGroundTruth();
};
