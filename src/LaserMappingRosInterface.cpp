#include "LaserMappingRosInterface.h"

#include "li_initialization.h"
#include "parameters.h"
#include "Urbannav_process/handler.h"

#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <sensor_msgs/PointCloud2.h>
#include <visualization_msgs/Marker.h>

#include <iostream>
#include <vector>

LaserMappingRosInterface::LaserMappingRosInterface(ros::NodeHandle &node)
{
  configureSubscribers(node);
  configurePublishers(node);
}

void LaserMappingRosInterface::loadGroundTruth()
{
  std::vector<Eigen::Vector4d> ground_truth;
  if (gt_file_type == LIVOX)
    GtfromTXT_LIVOX(LOCAL_FILE_DIR(gt_fname), ground_truth);
  else if (gt_file_type == URBAN)
    GtfromTXT_URBAN(LOCAL_FILE_DIR(gt_fname), ground_truth);
  else if (gt_file_type == M2DGR)
    GtfromTXT_M2DGR(LOCAL_FILE_DIR(gt_fname), ground_truth);

  std::cout << "ground-truth samples: " << ground_truth.size() << std::endl;
  for (const Eigen::Vector4d &sample : ground_truth)
  {
    if (gt_file_type == M2DGR)
      inputpvt_ecef(sample[0], sample[1], sample[2], sample[3],
                    p_gnss->first_lla_pvt, p_gnss->first_xyz_ecef_pvt,
                    p_gnss->pvt_time, p_gnss->pvt_holder,
                    p_gnss->diff_holder, p_gnss->float_holder);
    else
      inputpvt_lla(sample[0], sample[1], sample[2], sample[3],
                   p_gnss->first_lla_pvt, p_gnss->first_xyz_ecef_pvt,
                   p_gnss->pvt_time, p_gnss->pvt_holder,
                   p_gnss->diff_holder, p_gnss->float_holder);
  }
}

void LaserMappingRosInterface::configureSubscribers(ros::NodeHandle &node)
{
  lidar = p_pre->lidar_type == AVIA
      ? node.subscribe(lid_topic, 200000, livox_pcl_cbk)
      : node.subscribe(lid_topic, 200000, standard_pcl_cbk);
  imu = node.subscribe(imu_topic, 200000, imu_cbk);

  ephemeris = node.subscribe(gnss_ephem_topic, 10000, gnss_ephem_callback);
    glonass_ephemeris = node.subscribe(
        gnss_glo_ephem_topic, 10000, gnss_glo_ephem_callback);
    if (p_gnss->p_assign->obs_from_rinex)
    {
      gnss_measurement = node.subscribe(
          "/gnss_preprocessor_node/GNSSPsrCarRov1", 200,
          gnss_meas_callback_urbannav);
      rtk_pvt = node.subscribe(
          "/gnss_preprocessor_node/ECEFSolutionRTK", 500,
          rtklibOdomHandler);
    }
    else
    {
      gnss_measurement = node.subscribe(
          gnss_meas_topic, 10000, gnss_meas_callback);
      rtk_lla = node.subscribe(rtk_lla_topic, 1000, rtk_lla_callback);
    }

    if (p_gnss->p_assign->pvt_is_gt)
      rtk_pvt = node.subscribe(rtk_pvt_topic, 1000, rtk_pvt_callback);
    else
      loadGroundTruth();

    ionosphere = node.subscribe(
        gnss_iono_params_topic, 10000, gnss_iono_params_callback);
    if (gnss_local_online_sync)
    {
      gnss_time_pulse = node.subscribe(
          gnss_tp_info_topic, 100, gnss_tp_info_callback);
      local_trigger = node.subscribe(
          local_trigger_info_topic, 100, local_trigger_info_callback);
    }
    else
    {
      time_diff_gnss_local = gnss_local_time_diff;
      p_gnss->inputGNSSTimeDiff(time_diff_gnss_local);
      time_diff_valid = true;
  }
}

void LaserMappingRosInterface::configurePublishers(ros::NodeHandle &node)
{
  registered_cloud = node.advertise<sensor_msgs::PointCloud2>(
      "/cloud_registered", 1000);
  registered_body_cloud = node.advertise<sensor_msgs::PointCloud2>(
      "/cloud_registered_body", 1000);
  effective_cloud = node.advertise<sensor_msgs::PointCloud2>(
      "/cloud_effected", 1000);
  laser_map = node.advertise<sensor_msgs::PointCloud2>("/Laser_map", 1000);
  mapped_odometry = node.advertise<nav_msgs::Odometry>(
      "/aft_mapped_to_init", 1000);
  path = node.advertise<nav_msgs::Path>("/path", 1000);
  plane_marker = node.advertise<visualization_msgs::Marker>(
      "/planner_normal", 1000);
  rtk_satellites = node.advertise<visualization_msgs::Marker>(
      "/rtk_satellites", 100);
}
