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

// #include <so3_math.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <visualization_msgs/Marker.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/io/pcd_io.h>
#include <tf/transform_datatypes.h>
#include <tf/transform_broadcaster.h>
#include "li_initialization.h"
#include <malloc.h>
// #include <cv_bridge/cv_bridge.h>
#include <opencv2/opencv.hpp>
#include "chi-square.h"
#include "LaserMappingVisualization.h"
#include "LaserMappingRosInterface.h"
// #include <ros/console.h>

// This translation unit is the top-level orchestration layer for LIGO mapping.
// Sensor callbacks and ROS endpoint ownership live in LaserMappingRosInterface;
// this file consumes synchronized LiDAR, IMU, and GNSS data, advances the
// estimator in timestamp order, maintains the incremental voxel map, and
// publishes the resulting trajectory and point clouds.
//
// Processing lifecycle:
//   1. Configure ROS, sensor models, estimator covariance, and map storage.
//   2. Wait for a synchronized measurement package.
//   3. Initialize IMU attitude and build the first voxel map.
//   4. Interleave IMU/GNSS propagation with point-wise LiDAR corrections.
//   5. Insert accepted scan points and publish/log the current solution.

#define PUBFRAME_PERIOD     (20)

const float MOV_THRESHOLD = 1.5f;

string root_dir = ROOT_DIR;

int time_log_counter = 0;

// Lifecycle flags and the currently selected navigation observations.
bool init_map = false, flg_first_scan = true;
std::vector<ObsPtr> gnss_cur;
nav_msgs::OdometryPtr nmea_cur;
Eigen::Vector3d first_pvt_anc, first_lla_anc;
Eigen::Vector3d first_pvt_used, first_lla_used;

// These flags are also set by callbacks or signal handling. The mapping loop
// performs the actual reset/exit at a safe boundary.
bool  flg_reset = false, flg_exit = false;

// Scan buffers are expressed in the frame indicated by their suffix. The
// undistorted and body-frame clouds belong to the current scan; the world-frame
// initialization buffer persists until the initial map is large enough.
PointCloudXYZI::Ptr feats_undistort(new PointCloudXYZI());
PointCloudXYZI::Ptr feats_down_body_space(new PointCloudXYZI());
PointCloudXYZI::Ptr init_feats_world(new PointCloudXYZI());
std::deque<PointCloudXYZI::Ptr> depth_feats_world;
pcl::VoxelGrid<PointType> downSizeFilterSurf;

V3D euler_cur;

nav_msgs::Path path;
nav_msgs::Odometry odomAftMapped;
geometry_msgs::PoseStamped msg_body_pose;
RtkSatelliteVisualizer rtk_satellite_visualizer;

void SigHandle(int sig)
{
    flg_exit = true;
    ROS_WARN("catch sig %d", sig);
    sig_buffer.notify_all();
}

// Apply the calibrated LiDAR-to-IMU rigid transform to one point. Intensity is
// copied unchanged because it is not part of the geometric transformation.
void pointBodyLidarToIMU(PointType const * const pi, PointType * const po)
{
    V3D p_body_lidar(pi->x, pi->y, pi->z);
    V3D p_body_imu;
    {
        p_body_imu = Lidar_R_wrt_IMU * p_body_lidar + Lidar_T_wrt_IMU;
    }
    po->x = p_body_imu(0);
    po->y = p_body_imu(1);
    po->z = p_body_imu(2);
    po->intensity = pi->intensity;
}

// Insert only points that add new spatial support to the incremental voxel
// map. A point is skipped when a neighbor already occupies the same map voxel.
void MapIncremental() {
    PointVector points_to_add;
    int cur_pts = feats_down_world->size();
    points_to_add.reserve(cur_pts);

    for (size_t i = 0; i < cur_pts; ++i) {
        /* decide if need add to map */
        PointType &point_world = feats_down_world->points[i];
        if (!Nearest_Points[i].empty()) {
            const PointVector &points_near = Nearest_Points[i];

            Eigen::Vector3f center =
                ((point_world.getVector3fMap() / filter_size_map_min).array().floor() + 0.5) * filter_size_map_min;
            bool need_add = true;
            for (int readd_i = 0; readd_i < points_near.size(); readd_i++) {
                Eigen::Vector3f dis_2_center = points_near[readd_i].getVector3fMap() - center;
                if (fabs(dis_2_center.x()) < 0.5 * filter_size_map_min &&
                    fabs(dis_2_center.y()) < 0.5 * filter_size_map_min &&
                    fabs(dis_2_center.z()) < 0.5 * filter_size_map_min) {
                    need_add = false;
                    break;
                }
            }
            if (need_add) {
                points_to_add.emplace_back(point_world);
            }
        } else {
            points_to_add.emplace_back(point_world);
        }
    }
    ivox_->AddPoints(points_to_add);
}

// Publish the accumulated bootstrap cloud before it is transferred to IVox.
void publish_init_map(const ros::Publisher & pubLaserCloudFullRes)
{
    int size_init_map = init_feats_world->size();

    sensor_msgs::PointCloud2 laserCloudmsg;

    pcl::toROSMsg(*init_feats_world, laserCloudmsg);

    laserCloudmsg.header.stamp = ros::Time().fromSec(lidar_end_time);
    laserCloudmsg.header.frame_id = "camera_init";
    pubLaserCloudFullRes.publish(laserCloudmsg);
}

PointCloudXYZI::Ptr pcl_wait_pub(new PointCloudXYZI(500000, 1));
PointCloudXYZI::Ptr pcl_wait_save(new PointCloudXYZI());

// Publish the registered scan in the fixed map frame and, when requested,
// accumulate scans for periodic PCD output.
void publish_frame_world(const ros::Publisher & pubLaserCloudFullRes)
{
    if (scan_pub_en)
    {
        PointCloudXYZI::Ptr laserCloudFullRes(feats_down_body); // (points_num); //
        int size = laserCloudFullRes->points.size();

        PointCloudXYZI::Ptr   laserCloudWorld(new PointCloudXYZI(size, 1));

        for (int i = 0; i < size; i++)
        {
            // if (i % 3 == 0)
            {
            laserCloudWorld->points[i].x = feats_down_world->points[i].x; // updatedmap[i / 3](0); //
            laserCloudWorld->points[i].y = feats_down_world->points[i].y; // updatedmap[i / 3](1); //
            laserCloudWorld->points[i].z = feats_down_world->points[i].z; // updatedmap[i / 3](2); //
            laserCloudWorld->points[i].intensity = feats_down_world->points[i].intensity; // feats_down_world->points[i].y; // updatedmap[i / 3](2); //feats_down_world->points[i].z; //
            }
        }
        sensor_msgs::PointCloud2 laserCloudmsg;
        pcl::toROSMsg(*laserCloudWorld, laserCloudmsg);

        laserCloudmsg.header.stamp = ros::Time().fromSec(lidar_end_time); // (map_time); //
        laserCloudmsg.header.frame_id = "camera_init";
        pubLaserCloudFullRes.publish(laserCloudmsg);
        // publish_count -= PUBFRAME_PERIOD;
    }

    /**************** save map ****************/
    /* 1. make sure you have enough memories
    /* 2. noted that pcd save will influence the real-time performences **/
    if (pcd_save_en)
    {
        int size = points_num; // feats_down_world->points.size();
        PointCloudXYZI::Ptr   laserCloudWorld(new PointCloudXYZI(size, 1));

        for (int i = 0; i < size; i++)
        {
            laserCloudWorld->points[i].x = feats_down_world->points[i].x; // updatedmap[i](0); //
            laserCloudWorld->points[i].y = feats_down_world->points[i].y; // updatedmap[i](1); //
            laserCloudWorld->points[i].z = feats_down_world->points[i].z; // updatedmap[i](2); //
            laserCloudWorld->points[i].intensity = feats_down_world->points[i].intensity; // updatedmap[i](2); //
        }

        *pcl_wait_save += *laserCloudWorld;

        static int scan_wait_num = 0;
        scan_wait_num ++;
        if (pcl_wait_save->size() > 0 && pcd_save_interval > 0  && scan_wait_num >= pcd_save_interval)
        {
            pcd_index ++;
            string all_points_dir(string(string(ROOT_DIR) + "PCD/scans_") + to_string(pcd_index) + string(".pcd"));
            pcl::PCDWriter pcd_writer;
            cout << "current scan saved to /PCD/" << all_points_dir << endl;
            pcd_writer.writeBinary(all_points_dir, *pcl_wait_save);
            pcl_wait_save->clear();
            scan_wait_num = 0;
        }
    }
}

// Publish the undistorted scan in the IMU body frame for diagnostics.
void publish_frame_body(const ros::Publisher & pubLaserCloudFull_body)
{
    int size = feats_undistort->points.size();
    PointCloudXYZI::Ptr laserCloudIMUBody(new PointCloudXYZI(size, 1));

    for (int i = 0; i < size; i++)
    {
        pointBodyLidarToIMU(&feats_undistort->points[i], \
                            &laserCloudIMUBody->points[i]);
    }

    sensor_msgs::PointCloud2 laserCloudmsg;
    pcl::toROSMsg(*laserCloudIMUBody, laserCloudmsg);
    laserCloudmsg.header.stamp = ros::Time().fromSec(lidar_end_time);
    laserCloudmsg.header.frame_id = "body";
    pubLaserCloudFull_body.publish(laserCloudmsg);
}

template<typename T>
// Copy the current filter pose into any ROS pose-like message type.
void set_posestamp(T & out)
{
    {
        out.position.x = kf_output.x_.pos(0);
        out.position.y = kf_output.x_.pos(1);
        out.position.z = kf_output.x_.pos(2);
        Eigen::Quaterniond q(kf_output.x_.rot);
        out.orientation.x = q.coeffs()[0];
        out.orientation.y = q.coeffs()[1];
        out.orientation.z = q.coeffs()[2];
        out.orientation.w = q.coeffs()[3];
    }
}

// Publish both ROS odometry and the matching map-to-body TF transform. The
// timestamp policy follows publish_odometry_without_downsample.
void publish_odometry(const ros::Publisher & pubOdomAftMapped)
{
    odomAftMapped.header.frame_id = "camera_init";
    odomAftMapped.child_frame_id = "aft_mapped";
    if (publish_odometry_without_downsample)
    {
        odomAftMapped.header.stamp = ros::Time().fromSec(time_current);
    }
    else
    {
        odomAftMapped.header.stamp = ros::Time().fromSec(lidar_end_time);
    }
    set_posestamp(odomAftMapped.pose.pose);

    pubOdomAftMapped.publish(odomAftMapped);

    static tf::TransformBroadcaster br;
    tf::Transform                   transform;
    tf::Quaternion                  q;
    transform.setOrigin(tf::Vector3(odomAftMapped.pose.pose.position.x, \
                                    odomAftMapped.pose.pose.position.y, \
                                    odomAftMapped.pose.pose.position.z));
    q.setW(odomAftMapped.pose.pose.orientation.w);
    q.setX(odomAftMapped.pose.pose.orientation.x);
    q.setY(odomAftMapped.pose.pose.orientation.y);
    q.setZ(odomAftMapped.pose.pose.orientation.z);
    transform.setRotation( q );
    br.sendTransform( tf::StampedTransform( transform, odomAftMapped.header.stamp, "camera_init", "aft_mapped" ) );
}

// Append the latest filter pose to the complete path and publish it.
void publish_path(const ros::Publisher pubPath)
{
    set_posestamp(msg_body_pose.pose);
    // msg_body_pose.header.stamp = ros::Time::now();
    msg_body_pose.header.stamp = ros::Time().fromSec(lidar_end_time);
    msg_body_pose.header.frame_id = "camera_init";
    static int jjj = 0;
    jjj++;
    // if (jjj % 2 == 0) // if path is too large, the rvis will crash
    {
        path.poses.emplace_back(msg_body_pose);
        pubPath.publish(path);
    }
}

// Restore all stateful mapping components after a bag replay/time reset.
void resetMappingRuntime(Eigen::Matrix<double, 24, 24> &initial_covariance)
{
    ROS_WARN("reset when rosbag play back");
    p_imu->Reset();
    feats_undistort.reset(new PointCloudXYZI());
    state_out = state_output();
    kf_output.change_P(initial_covariance);
    is_first_gnss = true;
    flg_first_scan = true;
    is_first_frame = true;
    flg_reset = false;
    init_map = false;
    ivox_.reset(new IVoxType(ivox_options_));
    ivox_last_.reset(new IVoxType(ivox_options_));
    traj_manager.reset(new curvefitter::TrajectoryManager<4>());
}

// Establish the common time origin and discard IMU samples that precede the
// first LiDAR scan. When IMU processing is disabled, use configured gravity.
void initializeFirstScan()
{
    first_lidar_time = Measures.lidar_beg_time;
    flg_first_scan = false;
    if (first_imu_time < 1.0)
        first_imu_time = imu_next.header.stamp.toSec();
    time_current = 0.0;
    kf_output.x_.gravity << VEC_FROM_ARRAY(gravity);

    if (imu_en)
    {
        if (!nolidar && !imu_deque.empty())
        {
            while (Measures.lidar_beg_time > imu_next.header.stamp.toSec())
            {
                imu_deque.pop_front();
                if (imu_deque.empty()) break;
                imu_last = imu_next;
                imu_next = *imu_deque.front();
            }
        }
    }
    else
    {
        kf_output.x_.acc << VEC_FROM_ARRAY(gravity);
        kf_output.x_.acc *= -1;
        p_imu->imu_need_init_ = false;
    }
    G_m_s2 = std::sqrt(gravity[0] * gravity[0] +
                       gravity[1] * gravity[1] +
                       gravity[2] * gravity[2]);
}

// Undistort and optionally voxel-filter the current scan, then group its points
// by relative timestamp for chronological point-wise filter updates.
void preprocessCurrentScan()
{
    p_imu->Process(Measures, feats_undistort);
    if (space_down_sample)
    {
        downSizeFilterSurf.setInputCloud(feats_undistort);
        downSizeFilterSurf.filter(*feats_down_body);
    }
    else
    {
        feats_down_body = Measures.lidar;
    }
    std::sort(feats_down_body->points.begin(), feats_down_body->points.end(),
              time_list);
    if (!nolidar)
    {
        time_seq = time_compressing<int>(feats_down_body);
        feats_down_size = feats_down_body->points.size();
    }
    else
    {
        time_seq.clear();
    }
}

// Return false while the IMU initializer still needs observations. Once ready,
// initialize attitude from measured gravity (or the configured fallback).
bool initializeImuStateIfReady()
{
    if (p_imu->after_imu_init_) return true;
    if (p_imu->imu_need_init_) return false;

    V3D initial_gravity;
    if (init_with_imu && imu_en)
        initial_gravity = -p_imu->mean_acc / p_imu->mean_acc.norm() * G_m_s2;
    else
    {
        initial_gravity << VEC_FROM_ARRAY(gravity_init);
        p_imu->after_imu_init_ = true;
    }
    M3D initial_rotation;
    p_imu->Set_init(initial_gravity, initial_rotation);
    kf_output.x_.rot = initial_rotation;
    kf_output.x_.acc =
        -initial_rotation.transpose() * kf_output.x_.gravity;
    return true;
}

// Returns true when the current scan was consumed by map initialization and
// the caller should wait for the next synchronized package.
bool updateInitialMap(const ros::Publisher &map_publisher)
{
    if (init_map || nolidar || lose_lid) return false;

    feats_down_world->resize(feats_undistort->size());
    for (size_t index = 0; index < feats_undistort->size(); ++index)
        pointBodyToWorld(&feats_undistort->points[index],
                         &feats_down_world->points[index]);
    init_feats_world->points.insert(init_feats_world->points.end(),
                                    feats_down_world->points.begin(),
                                    feats_down_world->points.end());
    if (init_feats_world->size() >= init_map_size)
    {
        ivox_->AddPoints(init_feats_world->points);
        publish_init_map(map_publisher);
        init_feats_world.reset(new PointCloudXYZI());
        init_map = true;
        traj_manager->ResetTrajectory(pose_graph_key_pose, pose_time_vector,
                                      LiDAR_points, points_num);
    }
    return true;
}

// Application runner. The process entry point below intentionally remains a
// small composition root; this function owns the mapping lifecycle and keeps
// algorithm execution separate from process startup.
int runLaserMappingApplication(int argc, char** argv)
{
    // ---- Process and algorithm configuration ----
    ros::init(argc, argv, "laserMapping");
    ros::NodeHandle nh("~");
    ros::AsyncSpinner spinner(0);
    spinner.start();
    readParameters(nh);
    cout<<"lidar_type: "<<lidar_type<<endl;
    ivox_ = std::make_shared<IVoxType>(ivox_options_);
    ivox_last_ = std::make_shared<IVoxType>(ivox_options_); //(*ivox_);

    path.header.stamp    = ros::Time().fromSec(lidar_end_time);
    path.header.frame_id ="camera_init";

    // Running timing statistics are maintained only for diagnostic logging.
    int frame_num = 0;
    double aver_time_consu = 0, aver_time_icp = 0, aver_time_match = 0, aver_time_incre = 0, aver_time_solve = 0, aver_time_propag = 0;

    memset(point_selected_surf, true, sizeof(point_selected_surf));
    downSizeFilterSurf.setLeafSize(filter_size_surf_min, filter_size_surf_min, filter_size_surf_min);
    // Convert configured extrinsics once; all scan processing below uses these
    // matrix/vector forms.
    {
        Lidar_T_wrt_IMU<<VEC_FROM_ARRAY(extrinT);
        Lidar_R_wrt_IMU<<MAT_FROM_ARRAY(extrinR);
    }

    // ---- Sensor-model initialization ----
    p_imu->lidar_type = p_pre->lidar_type = lidar_type;
    p_imu->imu_en = imu_en;
    // Register process and measurement models with the shared error-state
    // Kalman filter, then install its initial covariance and process noise.
    {
        std::copy(default_gnss_iono_params.begin(), default_gnss_iono_params.end(),
            std::back_inserter(p_gnss->p_assign->latest_gnss_iono_params));
        p_gnss->Tex_imu_r << VEC_FROM_ARRAY(extrinT_gnss);
        p_gnss->gnss_ready = false; // gnss_quick_init; // edit
        p_gnss->nolidar = nolidar; // edit
        p_gnss->pre_integration->setnoise();

        if (p_gnss->p_assign->ephem_from_rinex)
        {
            p_gnss->p_assign->Ephemfromrinex(LOCAL_FILE_DIR(ephem_fname));
        }
    }

    {
        kf_output.init_dyn_share_modified_3h(get_f_output, df_dx_output, h_model_output, h_model_IMU_output, h_model_GNSS_output);
    }
    Eigen::Matrix<double, 24, 24> P_init_output; // = MD(24, 24)::Identity() * 0.01;
    reset_cov_output(P_init_output);
    kf_output.change_P(P_init_output);
    Eigen::Matrix<double, 24, 24> Q_output = process_noise_cov_output();
    open_file();

    // Own every subscriber and publisher in one lifecycle object. Keeping the
    // object in main preserves all callback connections for the process life.
    LaserMappingRosInterface ros_interface(nh);
    #ifdef process_ppp

    PPPfromTXT(LOCAL_FILE_DIR(ppp_fname), ppp_sol, ppp_ecef);

    #endif

    // Short aliases keep the numerical loop readable while ownership remains
    // explicit in ros_interface.
    const ros::Publisher &pubLaserCloudFullRes = ros_interface.registered_cloud;
    const ros::Publisher &pubLaserCloudFullRes_body =
        ros_interface.registered_body_cloud;
    const ros::Publisher &pubLaserCloudEffect = ros_interface.effective_cloud;
    const ros::Publisher &pubLaserCloudMap = ros_interface.laser_map;
    const ros::Publisher &pubOdomAftMapped = ros_interface.mapped_odometry;
    const ros::Publisher &pubPath = ros_interface.path;
    const ros::Publisher &plane_pub = ros_interface.plane_marker;
    const ros::Publisher &rtk_satellite_pub = ros_interface.rtk_satellites;

    // ---- Synchronized mapping loop ----
    signal(SIGINT, SigHandle);
    ros::Rate loop_rate(500);
    bool status = ros::ok();
    while (status)
    {
        if (flg_exit) break;
        ros::spinOnce();
        // sync_packages returns only when a processable LiDAR/IMU time window
        // (plus any GNSS/NMEA observations in that window) is available.
        if(sync_packages(Measures, p_gnss->gnss_msg, p_nmea->nmea_msg))
        {
            if (flg_reset)
            {
                resetMappingRuntime(P_init_output);
            }

            if (flg_first_scan)
            {
                initializeFirstScan();
            }

            double t0, t5;
            t0 = omp_get_wtime();

            // Initialization stages consume the current package and defer ICP
            // until gravity and the bootstrap map are both available.
            preprocessCurrentScan();
            if (!initializeImuStateIfReady())
            {
                continue;
            }
            if (updateInitialMap(pubLaserCloudMap))
            {
                continue;
            }

            // ---- Prepare point-wise LiDAR corrections ----
            // Allocate output/neighbor buffers and cache each point in the IMU
            // frame together with its skew matrix for measurement Jacobians.
            normvec->resize(feats_down_size);
            feats_down_world->resize(feats_down_size);

            Nearest_Points.resize(feats_down_size);
            // t2 = omp_get_wtime();

            crossmat_list.reserve(feats_down_size);
            pbody_list.reserve(feats_down_size);
            pimu_list.reserve(feats_down_size);
            // pbody_ext_list.reserve(feats_down_size);

            for (size_t i = 0; i < feats_down_body->size(); i++)
            {
                V3D point_this(feats_down_body->points[i].x,
                            feats_down_body->points[i].y,
                            feats_down_body->points[i].z);
                pbody_list[i]=point_this;
                {
                    point_this = Lidar_R_wrt_IMU * point_this + Lidar_T_wrt_IMU;
                    pimu_list[i] = point_this;
                }
                M3D point_crossmat;
                point_crossmat << SKEW_SYM_MATRX(point_this);
                crossmat_list[i]=point_crossmat;
            }

            effct_feat_num = 0;
            // A non-empty time sequence selects LiDAR-aided processing. Each
            // group contains points with the same relative scan timestamp.
            if (time_seq.size() > 0)
            {
                p_gnss->p_assign->process_feat_num += time_seq.size();
                p_gnss->nolidar_cur = false;

                double pcl_beg_time = Measures.lidar_beg_time;
                idx = -1;
                for (k = 0; k < time_seq.size(); k++)
                {
                    PointType &point_body  = feats_down_body->points[idx+time_seq[k]];

                    time_current = point_body.curvature / 1000.0 + pcl_beg_time;
                    if (time_current < time_predict_last_const)
                    {
                        continue;
                    }

                    // Seed propagation timestamps and inertial inputs before
                    // processing the first point group of the first scan.
                    if (is_first_frame)
                    {
                        if(imu_en && !imu_deque.empty())
                        {
                            while (time_current > imu_next.header.stamp.toSec())
                            {
                                imu_deque.pop_front();
                                if (imu_deque.empty()) break;
                                imu_last = imu_next;
                                imu_next = *(imu_deque.front());
                            }
                            angvel_avr<<imu_last.angular_velocity.x, imu_last.angular_velocity.y, imu_last.angular_velocity.z;
                            acc_avr   <<imu_last.linear_acceleration.x, imu_last.linear_acceleration.y, imu_last.linear_acceleration.z;
                            if (imu_deque.empty()) break;
                        }

                        // std::vector<Eigen::Vector3d>().swap(p_gnss->norm_vec_holder);
                        p_gnss->p_assign->process_feat_num = 0;
                        p_gnss->norm_vec_num = 0;
                        // acc_avr_norm = acc_avr * G_m_s2 / acc_norm;
                        // p_gnss->pre_integration->repropagate(kf_output.x_.ba, kf_output.x_.bg);
                        // p_gnss->pre_integration->setacc0gyr0(acc_avr_norm, angvel_avr);

                        is_first_frame = false;
                        time_update_last = time_current;
                        time_predict_last_const = time_current;
                    }
                    if(imu_en && !imu_deque.empty())
                    {
                        // Advance all IMU samples up to this LiDAR point time.
                        // GNSS epochs between two IMU samples are handled first,
                        // ensuring that measurements are never applied out of order.
                        bool last_imu = imu_next.header.stamp.toSec() == imu_deque.front()->header.stamp.toSec();
                        while (imu_next.header.stamp.toSec() < time_predict_last_const && !imu_deque.empty())
                        {
                            if (!last_imu)
                            {
                                imu_last = imu_next;
                                imu_next = *(imu_deque.front());
                                break;
                            }
                            else
                            {
                                imu_deque.pop_front();
                                if (imu_deque.empty()) break;
                                imu_last = imu_next;
                                imu_next = *(imu_deque.front());
                            }
                            if (imu_deque.empty()) break;
                        }
                        bool imu_comes = time_current >= imu_next.header.stamp.toSec();
                        while (imu_comes)
                        {
                            if (!p_gnss->gnss_msg.empty())
                            {
                                gnss_cur = p_gnss->gnss_msg.front();
                                // printf("%f, %f, %f\n", time2sec(gnss_cur[0]->time), time_diff_gnss_local, time_predict_last_const);
                                while (time2sec(gnss_cur[0]->time) - time_diff_gnss_local < time_predict_last_const)
                                {
                                    p_gnss->gnss_msg.pop();
                                    if(!p_gnss->gnss_msg.empty())
                                    {
                                        gnss_cur = p_gnss->gnss_msg.front();
                                    }
                                    else
                                    {
                                        break;
                                    }
                                }
                                if (p_gnss->gnss_msg.empty()) break;
                                while ((imu_next.header.stamp.toSec() >= time2sec(gnss_cur[0]->time) - time_diff_gnss_local) && (time2sec(gnss_cur[0]->time) - time_diff_gnss_local >= time_predict_last_const))
                                {
                                    double dt = time2sec(gnss_cur[0]->time) - time_diff_gnss_local - time_predict_last_const;
                                    double dt_cov = time2sec(gnss_cur[0]->time) - time_diff_gnss_local - time_update_last;

                                    // Before initialization, processGNSS builds the
                                    // GNSS anchor/state. Once ready, Evaluate decides
                                    // whether this epoch contributes a correction.
                                    if (p_gnss->gnss_ready)
                                    {
                                        if (dt_cov > 0.0)
                                        {
                                            kf_output.predict(dt_cov, Q_output, input_in, false, true);
                                        }
                                        kf_output.predict(dt, Q_output, input_in, true, false);
                                        // p_gnss->pre_integration->push_back(dt, kf_output.x_.acc + kf_output.x_.ba, kf_output.x_.omg + kf_output.x_.bg); // acc_avr, angvel_avr);
                                        // p_gnss->processIMUOutput(dt, kf_output.x_.acc, kf_output.x_.omg);
                                        time_predict_last_const = time2sec(gnss_cur[0]->time) - time_diff_gnss_local;
                                        time_update_last = time_predict_last_const;
                                        p_gnss->processGNSS(gnss_cur, kf_output.x_);
                                        p_gnss->sqrt_lidar = Eigen::LLT<Eigen::Matrix<double, 24, 24>>(kf_output.P_.inverse()).matrixL().transpose();
                                        // p_gnss->sqrt_lidar *= 0.002;
                                        update_gnss = p_gnss->Evaluate(kf_output.x_);
                                        if (!p_gnss->gnss_ready)
                                        {
                                            flg_reset = true;
                                            p_gnss->gnss_msg.pop();
                                            if(!p_gnss->gnss_msg.empty())
                                            {
                                                gnss_cur = p_gnss->gnss_msg.front();
                                            }
                                            break; // ?
                                        }

                                        if (update_gnss)
                                        {
                                            state_output out_state = kf_output.x_;
                                            kf_output.update_iterated_dyn_share_GNSS();
                                            Eigen::Vector3d pos_enu;
                                            if (!runtime_pos_log) cout_state_to_file(pos_enu);
                                            // sensor_msgs::NavSatFix gnss_lla_msg;
                                            // gnss_lla_msg.header.stamp = ros::Time().fromSec(time_current);
                                            // gnss_lla_msg.header.frame_id = "camera_init";
                                            // gnss_lla_msg.latitude = pos_enu(0);
                                            // gnss_lla_msg.longitude = pos_enu(1);
                                            // gnss_lla_msg.altitude = pos_enu(2);
                                            // pub_gnss_lla.publish(gnss_lla_msg);
                                            // A material GNSS correction invalidates
                                            // recent registered points. Refit their
                                            // trajectory and rebuild the affected map.
                                            if ((out_state.pos - kf_output.x_.pos).norm() > 0.1 && pose_graph_key_pose.size() > 4)
                                            {
                                                curvefitter::PoseData pose_data;
                                                pose_data.timestamp = time2sec(gnss_cur[0]->time) - time_diff_gnss_local;
                                                map_time = pose_data.timestamp;
                                                pose_data.orientation = Sophus::SO3d(Eigen::Quaterniond(kf_output.x_.rot).normalized().toRotationMatrix());
                                                pose_data.position = kf_output.x_.pos;
                                                if (map_time > pose_graph_key_pose.back().timestamp) // + 1e-9)
                                                {
                                                    pose_time_vector.push_back(pose_data.timestamp);
                                                    pose_graph_key_pose.emplace_back(pose_data);
                                                }
                                                else
                                                // else if (map_time == pose_time_vector.back())
                                                {
                                                    pose_data.timestamp = pose_graph_key_pose.back().timestamp;
                                                    pose_graph_key_pose.back() = pose_data;
                                                }
                                                // curvefitter::Trajectory<4> traj(0.1);
                                                // std::shared_ptr<curvefitter::Trajectory<4> > Traj_ptr = std::make_shared<curvefitter::Trajectory<4> >(traj);
                                                traj_manager->SetTrajectory(std::make_shared<curvefitter::Trajectory<4> >(0.025));
                                                traj_manager->FitCurve(pose_graph_key_pose[0].orientation.unit_quaternion(), pose_graph_key_pose[0].position, pose_time_vector[0], pose_time_vector.back(), pose_graph_key_pose);
                                                updatedmap.resize(points_num);
                                                updatedmap = traj_manager->GetUpdatedMapPoints(pose_time_vector, LiDAR_points);
                                                ivox_last_->AddPoints(updatedmap);
                                                ivox_->grids_map_ = ivox_last_->grids_map_;
                                                // for (auto &t : ivox_last_->grids_map_)
                                                // {
                                                    // ivox_->grids_map_[t.first] = (t.second);
                                                // }
                                                // ivox_ = std::make_shared<IVoxType>(*ivox_last_);
                                            }
                                            else
                                            {
                                                ivox_last_->grids_map_ = ivox_->grids_map_;
                                            }
                                            // reset_cov_output(kf_output.P_);
                                            traj_manager->ResetTrajectory(pose_graph_key_pose, pose_time_vector, LiDAR_points, points_num);
                                        }
                                    }
                                    else
                                    {
                                        if (dt_cov > 0.0)
                                        {
                                            kf_output.predict(dt_cov, Q_output, input_in, false, true);
                                        }

                                        kf_output.predict(dt, Q_output, input_in, true, false);

                                        time_predict_last_const = time2sec(gnss_cur[0]->time) - time_diff_gnss_local;
                                        time_update_last = time_predict_last_const;
                                        state_out = kf_output.x_;
                                        // state_out.rot = state_out.rot; //.normalized().toRotationMatrix();
                                        // state_out.rot.normalize();
                                        // state_out.pos = state_out.pos;
                                        // state_out.vel = state_out.vel;
                                        p_gnss->processGNSS(gnss_cur, state_out);
                                        if (p_gnss->gnss_ready)
                                        {
                                            // printf("time gnss ready: %f \n", time_predict_last_const);
                                            Eigen::Vector3d pos_enu;
                                            if (!runtime_pos_log) cout_state_to_file(pos_enu);
                                            // sensor_msgs::NavSatFix gnss_lla_msg;
                                            // gnss_lla_msg.header.stamp = ros::Time().fromSec(time_current);
                                            // gnss_lla_msg.header.frame_id = "camera_init";
                                            // gnss_lla_msg.latitude = pos_enu(0);
                                            // gnss_lla_msg.longitude = pos_enu(1);
                                            // gnss_lla_msg.altitude = pos_enu(2);
                                            // pub_gnss_lla.publish(gnss_lla_msg);
                                        }
                                    }
                                    p_gnss->gnss_msg.pop();
                                    if(!p_gnss->gnss_msg.empty())
                                    {
                                        gnss_cur = p_gnss->gnss_msg.front();
                                    }
                                    else
                                    {
                                        break;
                                    }
                                }
                            }

                            if (flg_reset)
                            {
                                break;
                            }
                            angvel_avr<<imu_next.angular_velocity.x, imu_next.angular_velocity.y, imu_next.angular_velocity.z;
                            acc_avr   <<imu_next.linear_acceleration.x, imu_next.linear_acceleration.y, imu_next.linear_acceleration.z;

                            // Propagate both covariance and nominal state to this
                            // IMU timestamp before applying its measurement update.
                            double dt = imu_next.header.stamp.toSec() - time_predict_last_const;
                            time_predict_last_const = imu_next.header.stamp.toSec();
                            double dt_cov = imu_next.header.stamp.toSec() - time_update_last;

                            if (dt_cov > 0.0)
                            {
                                time_update_last = imu_next.header.stamp.toSec();

                                kf_output.predict(dt_cov, Q_output, input_in, false, true);
                            }
                            kf_output.predict(dt, Q_output, input_in, true, false);
                            kf_output.update_iterated_dyn_share_IMU();
                            imu_deque.pop_front();
                            if (imu_deque.empty()) break;
                            imu_last = imu_next;
                            imu_next = *(imu_deque.front());
                            imu_comes = time_current >= imu_next.header.stamp.toSec();
                        }
                    }
                    if (flg_reset)
                    {
                        break;
                    }

                    if (!p_gnss->gnss_msg.empty())
                    {
                        gnss_cur = p_gnss->gnss_msg.front();
                        // printf("%f, %f, %f\n", time2sec(gnss_cur[0]->time), time_diff_gnss_local, time_predict_last_const);
                        while ( time2sec(gnss_cur[0]->time) - time_diff_gnss_local < time_predict_last_const)
                        {
                            p_gnss->gnss_msg.pop();
                            if(!p_gnss->gnss_msg.empty())
                            {
                                gnss_cur = p_gnss->gnss_msg.front();
                            }
                            else
                            {
                                break;
                            }
                        }
                        if (p_gnss->gnss_msg.empty()) break;
                        while (time_current >= time2sec(gnss_cur[0]->time) - time_diff_gnss_local && time2sec(gnss_cur[0]->time) - time_diff_gnss_local >= time_predict_last_const)
                        {
                            double dt = time2sec(gnss_cur[0]->time) - time_diff_gnss_local - time_predict_last_const;
                            double dt_cov = time2sec(gnss_cur[0]->time) - time_diff_gnss_local - time_update_last;
                            // cout << "check gnss ready:" << p_gnss->gnss_ready << endl;
                            if (p_gnss->gnss_ready)
                            {
                                if (dt_cov > 0.0)
                                {
                                    kf_output.predict(dt_cov, Q_output, input_in, false, true);
                                }
                                kf_output.predict(dt, Q_output, input_in, true, false);

                                // p_gnss->pre_integration->push_back(dt, kf_output.x_.acc + kf_output.x_.ba, kf_output.x_.omg + kf_output.x_.bg); // acc_avr, angvel_avr);
                                // p_gnss->processIMUOutput(dt, kf_output.x_.acc, kf_output.x_.omg);

                                time_predict_last_const = time2sec(gnss_cur[0]->time) - time_diff_gnss_local;
                                time_update_last = time_predict_last_const;
                                p_gnss->processGNSS(gnss_cur, kf_output.x_);
                                p_gnss->sqrt_lidar = Eigen::LLT<Eigen::Matrix<double, 24, 24>>(kf_output.P_.inverse()).matrixL().transpose();
                                // p_gnss->sqrt_lidar *= 0.002;
                                update_gnss = p_gnss->Evaluate(kf_output.x_);
                                if (!p_gnss->gnss_ready)
                                {
                                    flg_reset = true;
                                    p_gnss->gnss_msg.pop();
                                    if(!p_gnss->gnss_msg.empty())
                                    {
                                        gnss_cur = p_gnss->gnss_msg.front();
                                    }
                                    break; // ?
                                }

                                if (update_gnss)
                                {
                                    state_output out_state = kf_output.x_;
                                    kf_output.update_iterated_dyn_share_GNSS();
                                    // reset_cov_output(kf_output.P_);
                                    Eigen::Vector3d pos_enu;
                                    if (!runtime_pos_log) cout_state_to_file(pos_enu);
                                    // sensor_msgs::NavSatFix gnss_lla_msg;
                                    // gnss_lla_msg.header.stamp = ros::Time().fromSec(time_current);
                                    // gnss_lla_msg.header.frame_id = "camera_init";
                                    // gnss_lla_msg.latitude = pos_enu(0);
                                    // gnss_lla_msg.longitude = pos_enu(1);
                                    // gnss_lla_msg.altitude = pos_enu(2);
                                    // pub_gnss_lla.publish(gnss_lla_msg);
                                    if ((out_state.pos - kf_output.x_.pos).norm() > 0.1 && pose_graph_key_pose.size() > 4)
                                    {
                                        curvefitter::PoseData pose_data;
                                        pose_data.timestamp = time2sec(gnss_cur[0]->time) - time_diff_gnss_local;
                                        map_time = pose_data.timestamp;
                                        // pose_time_vector.push_back(pose_data.timestamp);
                                        pose_data.orientation = Sophus::SO3d(Eigen::Quaterniond(kf_output.x_.rot).normalized().toRotationMatrix());
                                        pose_data.position = kf_output.x_.pos;
                                        if (map_time > pose_graph_key_pose.back().timestamp) // + 1e-9)
                                        {
                                            pose_time_vector.push_back(pose_data.timestamp);
                                            pose_graph_key_pose.emplace_back(pose_data);
                                        }
                                        else
                                        // else if (map_time == pose_time_vector.back())
                                        {
                                            pose_data.timestamp = pose_graph_key_pose.back().timestamp;
                                            pose_graph_key_pose.back() = pose_data;
                                        }
                                        // pose_graph_key_pose.emplace_back(pose_data);
                                        traj_manager->SetTrajectory(std::make_shared<curvefitter::Trajectory<4> >(0.025));
                                        traj_manager->FitCurve(pose_graph_key_pose[0].orientation.unit_quaternion(), pose_graph_key_pose[0].position, pose_time_vector[0], pose_time_vector.back(), pose_graph_key_pose);
                                        updatedmap.resize(points_num);
                                        updatedmap = traj_manager->GetUpdatedMapPoints(pose_time_vector, LiDAR_points);
                                        ivox_last_->AddPoints(updatedmap);
                                        ivox_->grids_map_ = ivox_last_->grids_map_;
                                    }
                                    else
                                    {
                                        ivox_last_->grids_map_ = ivox_->grids_map_;
                                    }
                                    traj_manager->ResetTrajectory(pose_graph_key_pose, pose_time_vector, LiDAR_points, points_num);
                                }
                            }
                            else
                            {
                                if (dt_cov > 0.0)
                                {
                                    kf_output.predict(dt_cov, Q_output, input_in, false, true);
                                }
                                kf_output.predict(dt, Q_output, input_in, true, false);
                                time_predict_last_const = time2sec(gnss_cur[0]->time) - time_diff_gnss_local;
                                time_update_last = time_predict_last_const;
                                state_out = kf_output.x_;
                                // state_out.rot = state_out.rot; //.normalized().toRotationMatrix();
                                // state_out.rot.normalize();
                                // state_out.pos = state_out.pos;
                                // state_out.vel = state_out.vel;
                                p_gnss->processGNSS(gnss_cur, state_out);
                                if (p_gnss->gnss_ready)
                                {
                                    // printf("time gnss ready: %f \n", time_predict_last_const);
                                    Eigen::Vector3d pos_enu;
                                    if (!runtime_pos_log) cout_state_to_file(pos_enu);
                                    // sensor_msgs::NavSatFix gnss_lla_msg;
                                    // gnss_lla_msg.header.stamp = ros::Time().fromSec(time_current);
                                    // gnss_lla_msg.header.frame_id = "camera_init";
                                    // gnss_lla_msg.latitude = pos_enu(0);
                                    // gnss_lla_msg.longitude = pos_enu(1);
                                    // gnss_lla_msg.altitude = pos_enu(2);
                                    // pub_gnss_lla.publish(gnss_lla_msg);
                                }
                            }
                            p_gnss->gnss_msg.pop();
                            if(!p_gnss->gnss_msg.empty())
                            {
                                gnss_cur = p_gnss->gnss_msg.front();
                            }
                            else
                            {
                                break;
                            }
                        }
                    }

                    if (flg_reset)
                    {
                        break;
                    }
                    // Bring the state exactly to the LiDAR point-group timestamp
                    // before evaluating the point-to-map measurement model.
                    double dt = time_current - time_predict_last_const;
                    // double propag_state_start = omp_get_wtime();
                    if(!prop_at_freq_of_imu)
                    {
                        double dt_cov = time_current - time_update_last;
                        if (dt_cov > 0.0)
                        {
                            kf_output.predict(dt_cov, Q_output, input_in, false, true);
                            time_update_last = time_current;
                        }
                    }
                    // if (dt > 0.0)

                    kf_output.predict(dt, Q_output, input_in, true, false);
                    time_predict_last_const = time_current;
                    if (feats_down_size < 1)
                    {
                        ROS_WARN("No point, skip this scan!\n");
                        idx += time_seq[k];
                        continue;
                    }
                    if (!kf_output.update_iterated_dyn_share_modified())
                    {
                        idx = idx+time_seq[k];
                        continue;
                    }

                    // else
                    // {
                    //     idx = idx+time_seq[k];
                    //     continue;
                    // }

                    // solve_start = omp_get_wtime();

                    if (publish_odometry_without_downsample)
                    {
                        /******* Publish odometry *******/

                        publish_odometry(pubOdomAftMapped);
                        if (runtime_pos_log)
                        {
                            euler_cur = SO3ToEuler(kf_output.x_.rot);
                            fout_out << setw(20) << Measures.lidar_beg_time - first_lidar_time << " " << kf_output.x_.pos.transpose() << " " << euler_cur.transpose() << " " << kf_output.x_.vel.transpose() \
                            <<" "<<kf_output.x_.omg.transpose()<<" "<<kf_output.x_.acc.transpose()<<" "<<kf_output.x_.gravity.transpose()<<" "<<kf_output.x_.bg.transpose()<<" "<<kf_output.x_.ba.transpose() << endl;
                        }
                    }
                    // Retain corrected poses and points so a later GNSS update can
                    // refit the local trajectory and rebuild the corresponding map.
                    std::vector<Eigen::Vector3d> lidarpoints;
                    for (int j = 0; j < time_seq[k]; j++)
                    {
                        PointType &point_body_j  = feats_down_body->points[idx+j+1];
                        PointType &point_world_j = feats_down_world->points[idx+j+1];
                        pointBodyToWorld(&point_body_j, &point_world_j);
                        lidarpoints.push_back(pimu_list[idx+j+1]); // (Eigen::Vector3d(point_body_j.x, point_body_j.y, point_body_j.z));
                    }

                    if (pose_graph_key_pose.empty()){
                        traj_manager->AddGraphPose(Eigen::Quaterniond(kf_output.x_.rot).normalized(), kf_output.x_.pos, lidarpoints, time_current, pose_graph_key_pose, pose_time_vector, LiDAR_points, points_num);
                    }
                    else
                    {
                        if (time_current > pose_graph_key_pose.back().timestamp && lidarpoints.size() > 0)
                            traj_manager->AddGraphPose(Eigen::Quaterniond(kf_output.x_.rot).normalized(), kf_output.x_.pos, lidarpoints, time_current, pose_graph_key_pose, pose_time_vector, LiDAR_points, points_num);
                    }

                    idx += time_seq[k];
                }
            }
            else
            {
                // Without LiDAR point groups, continue chronological IMU/GNSS
                // processing for GNSS-only operation or temporary LiDAR loss.
                p_gnss->nolidar_cur = true;

                if (!imu_deque.empty())
                {
                    imu_last = imu_next;
                    imu_next = *(imu_deque.front());

                while (imu_next.header.stamp.toSec() > time_current && ((imu_next.header.stamp.toSec() < imu_first_time + lidar_time_inte && nolidar) || (imu_next.header.stamp.toSec() < Measures.lidar_beg_time + lidar_time_inte && !nolidar)))
                { // >= ?
                    if (is_first_frame)
                    {
                        if (!nolidar) //std::vector<Eigen::Vector3d>().swap(p_gnss->norm_vec_holder);
                        {p_gnss->p_assign->process_feat_num = 0;
                        p_gnss->norm_vec_num = 0;}

                        if (!p_gnss->gnss_msg.empty())
                        {
                            gnss_cur = p_gnss->gnss_msg.front();
                            double front_gnss_ts = time2sec(gnss_cur[0]->time); // take time
                            time_current = front_gnss_ts - time_diff_gnss_local;
                            while (imu_next.header.stamp.toSec() < time_current) // 0.05
                            {
                                ROS_WARN("throw IMU, only should happen at the beginning 2510");
                                imu_deque.pop_front();
                                if (imu_deque.empty()) break;
                                imu_last = imu_next;
                                imu_next = *(imu_deque.front()); // could be used to initialize
                            }
                            if (imu_deque.empty()) break;
                        }
                        else
                        {
                            if (nolidar)
                            {
                                while (imu_next.header.stamp.toSec() < imu_first_time + lidar_time_inte)
                                {
                                    // meas.imu.emplace_back(imu_deque.front()); should add to initialization
                                    imu_deque.pop_front();
                                    if(imu_deque.empty()) break;
                                    imu_last = imu_next;
                                    imu_next = *(imu_deque.front()); // could be used to initialize
                                }
                                // if (imu_deque.empty()) break;
                            }
                            else
                            {
                                while (imu_next.header.stamp.toSec() < Measures.lidar_beg_time + lidar_time_inte)
                                {
                                    // meas.imu.emplace_back(imu_deque.front()); should add to initialization
                                    imu_deque.pop_front();
                                    if(imu_deque.empty()) break;
                                    imu_last = imu_next;
                                    imu_next = *(imu_deque.front());
                                }
                            }
                            break;
                        }
                        angvel_avr<<imu_last.angular_velocity.x, imu_last.angular_velocity.y, imu_last.angular_velocity.z;
                        if (nolidar) kf_output.x_.omg = angvel_avr;

                        acc_avr   <<imu_last.linear_acceleration.x, imu_last.linear_acceleration.y, imu_last.linear_acceleration.z;
                        time_current = imu_next.header.stamp.toSec();

                        time_update_last = time_current;
                        time_predict_last_const = time_current;
                        acc_avr_norm = acc_avr * G_m_s2 / acc_norm;
                        {
                        p_gnss->pre_integration->repropagate(kf_output.x_.ba, kf_output.x_.bg);
                        p_gnss->pre_integration->setacc0gyr0(acc_avr_norm, angvel_avr);
                        }

                        {
                            is_first_frame = false;
                        }
                    }
                    time_current = imu_next.header.stamp.toSec();

                    if (!is_first_frame)
                    {
                    if (!p_gnss->gnss_msg.empty())
                    {
                        gnss_cur = p_gnss->gnss_msg.front();
                        while (time2sec(gnss_cur[0]->time) - time_diff_gnss_local <= time_predict_last_const)
                        {
                            p_gnss->gnss_msg.pop();
                            if(!p_gnss->gnss_msg.empty())
                            {
                                gnss_cur = p_gnss->gnss_msg.front();
                            }
                            else
                            {
                                break;
                            }
                        }
                        if (p_gnss->gnss_msg.empty()) break;
                    while ((time_current > time2sec(gnss_cur[0]->time) - time_diff_gnss_local) && (time2sec(gnss_cur[0]->time) - time_diff_gnss_local > time_predict_last_const))
                    {
                        double dt = time2sec(gnss_cur[0]->time) - time_diff_gnss_local - time_predict_last_const;
                        double dt_cov = time2sec(gnss_cur[0]->time) - time_diff_gnss_local - time_update_last;

                        if (p_gnss->gnss_ready)
                        {
                            if (dt_cov > 0.0)
                            {
                                // kf_output.predict(dt_cov, Q_output, input_in, false, true);
                                time_update_last = time2sec(gnss_cur[0]->time) - time_diff_gnss_local;
                            }
                            // kf_output.predict(dt, Q_output, input_in, true, false);
                            p_gnss->pre_integration->push_back(dt, acc_avr_norm, angvel_avr); //acc_avr_norm, angvel_avr);
                            // change to state_const.omg and state_const.acc?
                            time_predict_last_const = time2sec(gnss_cur[0]->time) - time_diff_gnss_local;
                            p_gnss->processGNSS(gnss_cur, kf_output.x_);
                            if (!nolidar)
                            {
                                p_gnss->sqrt_lidar = Eigen::LLT<Eigen::Matrix<double, 24, 24>>(kf_output.P_.inverse()).matrixL().transpose();
                            }
                            update_gnss = p_gnss->Evaluate(kf_output.x_);
                            if (!p_gnss->gnss_ready)
                            {
                                flg_reset = true;
                                p_gnss->gnss_msg.pop();
                                if(!p_gnss->gnss_msg.empty())
                                {
                                    gnss_cur = p_gnss->gnss_msg.front();
                                }
                                break; // ?
                            }
                            if (update_gnss)
                            {
                                if (!nolidar)
                                {
                                    state_output out_state = kf_output.x_;
                                    kf_output.update_iterated_dyn_share_GNSS();
                                    // reset_cov_output(kf_output.P_);
                                    if ((out_state.pos - kf_output.x_.pos).norm() > 0.1 && pose_graph_key_pose.size() > 4)
                                    {
                                        curvefitter::PoseData pose_data;
                                        pose_data.timestamp = time2sec(gnss_cur[0]->time) - time_diff_gnss_local;
                                        map_time = pose_data.timestamp;
                                        // pose_time_vector.push_back(pose_data.timestamp);
                                        pose_data.orientation = Sophus::SO3d(Eigen::Quaterniond(kf_output.x_.rot).normalized().toRotationMatrix());
                                        pose_data.position = kf_output.x_.pos;
                                        if (map_time > pose_graph_key_pose.back().timestamp) // + 1e-9)
                                        {
                                            pose_time_vector.push_back(pose_data.timestamp);
                                            pose_graph_key_pose.emplace_back(pose_data);
                                        }
                                        // else if (map_time == pose_time_vector.back())
                                        else
                                        {
                                            pose_data.timestamp = pose_graph_key_pose.back().timestamp;
                                            pose_graph_key_pose.back() = pose_data;
                                        }
                                        // pose_graph_key_pose.emplace_back(pose_data);
                                        // curvefitter::Trajectory<4> traj(0.1);
                                        // std::shared_ptr<curvefitter::Trajectory<4> > Traj_ptr = std::make_shared<curvefitter::Trajectory<4> >(traj);
                                        traj_manager->SetTrajectory(std::make_shared<curvefitter::Trajectory<4> >(0.025));
                                        traj_manager->FitCurve(pose_graph_key_pose[0].orientation.unit_quaternion(), pose_graph_key_pose[0].position, pose_time_vector[0], pose_time_vector.back(), pose_graph_key_pose);
                                        updatedmap.resize(points_num);
                                        updatedmap = traj_manager->GetUpdatedMapPoints(pose_time_vector, LiDAR_points);
                                        ivox_last_->AddPoints(updatedmap);
                                        ivox_->grids_map_ = ivox_last_->grids_map_;
                                        // for (auto &t : ivox_last_->grids_map_)
                                        // {
                                            // (ivox_->grids_map_[t.first]) = (t.second);
                                        // }
                                        // ivox_ = std::make_shared<IVoxType>(*ivox_last_);
                                    }
                                    else
                                    {
                                        ivox_last_->grids_map_ = ivox_->grids_map_;
                                    }
                                    traj_manager->ResetTrajectory(pose_graph_key_pose, pose_time_vector, LiDAR_points, points_num);
                                }
                                Eigen::Vector3d pos_enu;
                                if (!runtime_pos_log) cout_state_to_file(pos_enu);
                                // sensor_msgs::NavSatFix gnss_lla_msg;
                                // gnss_lla_msg.header.stamp = ros::Time().fromSec(time_current);
                                // gnss_lla_msg.header.frame_id = "camera_init";
                                // gnss_lla_msg.latitude = pos_enu(0);
                                // gnss_lla_msg.longitude = pos_enu(1);
                                // gnss_lla_msg.altitude = pos_enu(2);
                                // pub_gnss_lla.publish(gnss_lla_msg);
                            }
                        }
                        else
                        {
                            if (dt_cov > 0.0)
                            {
                                // kf_output.predict(dt_cov, Q_output, input_in, false, true);
                                time_update_last = time2sec(gnss_cur[0]->time) - time_diff_gnss_local;
                            }
                            // kf_output.predict(dt, Q_output, input_in, true, false);
                            time_predict_last_const = time2sec(gnss_cur[0]->time) - time_diff_gnss_local;
                            p_gnss->processGNSS(gnss_cur, kf_output.x_);
                            if (p_gnss->gnss_ready)
                            {
                                Eigen::Vector3d pos_enu;
                                if (!runtime_pos_log) cout_state_to_file(pos_enu);
                                // printf("time gnss ready: %f \n", time_predict_last_const);
                                // sensor_msgs::NavSatFix gnss_lla_msg;
                                // gnss_lla_msg.header.stamp = ros::Time().fromSec(time_current);
                                // gnss_lla_msg.header.frame_id = "camera_init";
                                // gnss_lla_msg.latitude = pos_enu(0);
                                // gnss_lla_msg.longitude = pos_enu(1);
                                // gnss_lla_msg.altitude = pos_enu(2);
                                // pub_gnss_lla.publish(gnss_lla_msg);
                                if (nolidar)
                                {
                                    // Eigen::Matrix3d R_enu_local_;
                                    // R_enu_local_ = Eigen::AngleAxisd(p_gnss->yaw_enu_local, Eigen::Vector3d::UnitZ());
                                    kf_output.x_.pos = p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(p_gnss->frame_num-1)).segment<3>(0); // p_gnss->anc_ecef - p_gnss->R_ecef_enu * R_enu_local_ * state_const.rot_end * p_gnss->Tex_imu_r;
                                    kf_output.x_.rot = p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Rot3>(R(p_gnss->frame_num-1)).matrix(); // p_gnss->R_ecef_enu * R_enu_local_ * state_const.rot_end;
                                    // kf_output.x_.rot.normalize();
                                    kf_output.x_.vel = p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(p_gnss->frame_num-1)).segment<3>(3); // p_gnss->R_ecef_enu * R_enu_local_ * state_const.vel_end; // Eigen::Vector3d::Zero(); // R_ecef_enu * state.vel_end;
                                    kf_output.x_.ba = Eigen::Vector3d::Zero(); // R_ecef_enu * state.vel_end;
                                    kf_output.x_.bg = Eigen::Vector3d::Zero(); // R_ecef_enu * state.vel_end;
                                    kf_output.x_.omg = Eigen::Vector3d::Zero(); // R_ecef_enu * state.vel_end;
                                    kf_output.x_.gravity = p_gnss->R_ecef_enu * kf_output.x_.gravity; // * R_enu_local_
                                    kf_output.x_.acc = kf_output.x_.rot.transpose() * (-kf_output.x_.gravity); // R_ecef_enu * state.vel_end;.conjugate().normalized()

                                    kf_output.P_ = MD(24,24)::Identity() * INIT_COV;
                                }
                            }
                        }
                        p_gnss->gnss_msg.pop();
                        if(!p_gnss->gnss_msg.empty())
                        {
                            gnss_cur = p_gnss->gnss_msg.front();
                        }
                        else
                        {
                            break;
                        }
                    }
                    }

                    if (flg_reset)
                    {
                        break;
                    }
                    double dt = time_current - time_predict_last_const;
                    {
                        double dt_cov = time_current - time_update_last;
                        if (dt_cov > 0.0)
                        {
                            // kf_output.predict(dt_cov, Q_output, input_in, false, true);
                            time_update_last = time_current;
                        }
                        // kf_output.predict(dt, Q_output, input_in, true, false);
                        p_gnss->pre_integration->push_back(dt, acc_avr_norm, angvel_avr); // acc_avr_norm, angvel_avr); //
                    }

                    time_predict_last_const = time_current;

                    angvel_avr<<imu_next.angular_velocity.x, imu_next.angular_velocity.y, imu_next.angular_velocity.z;
                    if (nolidar) kf_output.x_.omg = angvel_avr;
                    acc_avr   <<imu_next.linear_acceleration.x, imu_next.linear_acceleration.y, imu_next.linear_acceleration.z;
                    acc_avr_norm = acc_avr * G_m_s2 / acc_norm;
                    kf_output.update_iterated_dyn_share_IMU();
                    imu_deque.pop_front();
                    if (imu_deque.empty()) break;
                    imu_last = imu_next;
                    imu_next = *(imu_deque.front());
                }
                else
                {
                    imu_deque.pop_front();
                    if (imu_deque.empty()) break;
                    imu_last = imu_next;
                    imu_next = *(imu_deque.front());
                }
                }
                }
            }

            // ---- Map maintenance and outputs ----
            // In scan-rate mode, publish odometry after all point corrections.
            if (!publish_odometry_without_downsample)
            {
                publish_odometry(pubOdomAftMapped);
            }

            // Insert the registered scan only after the final state correction.
            if(feats_down_size > 4)
            {
                MapIncremental();
            }

            t5 = omp_get_wtime();
            // Publish optional visualization products from the same final state.
            if (path_en)                         publish_path(pubPath);
            if (scan_pub_en || pcd_save_en)      publish_frame_world(pubLaserCloudFullRes);
            if (scan_pub_en && scan_body_pub_en) publish_frame_body(pubLaserCloudFullRes_body);
            if (p_gnss->rtk_satellite_visualization)
            {
                const Eigen::Vector3d gnss_antenna_local =
                    kf_output.x_.pos + kf_output.x_.rot * p_gnss->Tex_imu_r;
                rtk_satellite_visualizer.publish(
                    rtk_satellite_pub, *p_gnss, gnss_antenna_local,
                    lidar_end_time);
            }

            // Record runtime statistics and trajectory samples when enabled.
            if (runtime_pos_log)
            {
                frame_num ++;
                aver_time_consu = aver_time_consu * (frame_num - 1) / frame_num + (t5 - t0) / frame_num;
                s_plot[time_log_counter] = t5 - t0;
                s_plot3[time_log_counter] = aver_time_consu;
                time_log_counter ++;
                if (!publish_odometry_without_downsample)
                {
                    {
                        {
                            Eigen::Matrix3d R_enu_local_;
                            Eigen::Vector3d pos_r = kf_output.x_.rot * p_gnss->Tex_imu_r + kf_output.x_.pos; // .normalized()
                            time_frame.push_back(lidar_end_time); //(time_predict_last_const);
                            est_poses.push_back(pos_r);
                        }
                        euler_cur = SO3ToEuler(kf_output.x_.rot);
                        fout_out << setw(20) << Measures.lidar_beg_time - first_lidar_time << " " << kf_output.x_.pos.transpose() << " " << euler_cur.transpose() << " " << kf_output.x_.vel.transpose() \
                        <<" "<<state_out.omg.transpose()<<" "<<state_out.acc.transpose()<<" "<<state_out.gravity.transpose()<<" "<<state_out.bg.transpose()<<" "<<state_out.ba.transpose() << endl;
                    }
                }
            }
        }
        status = ros::ok();
        loop_rate.sleep();
    }
    // ---- Orderly shutdown and final result export ----
    fout_out.close();
    // Saving a large accumulated cloud can require substantial memory and may
    // affect real-time performance, so it is guarded by pcd_save_en.
    if (pcl_wait_save->size() > 0 && pcd_save_en)
    {
        string file_name = string("scans.pcd");
        string all_points_dir(string(string(ROOT_DIR) + "PCD/") + file_name);
        pcl::PCDWriter pcd_writer;
        pcd_writer.writeBinary(all_points_dir, *pcl_wait_save);
    }
    {
        Eigen::Matrix3d enu_rot = ecef2rotation(first_pvt_used);
        for (int i = 0; i < time_frame.size(); i++)
        {
            // Eigen::Vector3d euler_ext = SO3ToEuler(local_rots[i]);
            {
                fout_global << setw(20) << time_frame[i] - time_frame[0] << " " << est_poses[i].transpose() << endl; // << " " << local_poses[i].transpose() << " " << euler_ext.transpose() << endl; //"\n"; // p_gnss->pvt_time[0] + 18.0
                // fout_global << setw(20) << time_frame[i] - ppp_ecef[0][0] + 18.0 << " " << est_poses[i].transpose() << endl; //"\n"; // p_gnss->pvt_time[0] + 18.0
            }
            // printf("time: %f, pos: %f %f %f\n", ppp_ecef[0][0] + 18.0, est_poses[i](0), est_poses[i](1), est_poses[i](2));
            // Eigen::Vector3d euler_ext = SO3ToEuler(local_rots[i]);
        }
        fout_global.close();
    }

    for (int i = 0; i < p_gnss->pvt_time.size(); i++)
    {
        fout_rtk << setw(20) << p_gnss->pvt_time[i] - p_gnss->pvt_time[0] << " " << p_gnss->pvt_holder[i].transpose() << " " << p_gnss->diff_holder[i] << " " << p_gnss->float_holder[i] << endl; // "\n";
    }
    fout_rtk.close();
    #ifdef process_ppp
    for (int i = 0; i < ppp_ecef.size(); i++)
    {
        Eigen::Vector3d pos_enu = ecef2enu(p_gnss->first_lla_pvt, ppp_ecef[i].segment<3>(1) - p_gnss->first_xyz_ecef_pvt);
        fout_ppp << setw(20) << ppp_ecef[i][0] - ppp_ecef[0][0] << " " << pos_enu.transpose() << endl;
    }
    fout_ppp.close();
    #endif

    return 0;
}

int main(int argc, char **argv)
{
    // Keep the executable entry point free of mapping policy; the application
    // runner owns initialization, processing, and shutdown.
    return runLaserMappingApplication(argc, argv);
}
