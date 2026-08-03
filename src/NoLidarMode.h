#pragma once

struct NoLidarWindow
{
    double begin_time = 0.0;
    double end_time = 0.0;
};

// Build a fixed-duration processing window driven entirely by IMU timestamps.
// Returns false until the IMU buffer covers the complete window.
bool makeNoLidarWindow(double first_imu_time, double latest_imu_time,
                       double window_duration, NoLidarWindow &window);

// Convert a GNSS timestamp into the local sensor clock and determine whether
// it belongs to the current IMU-driven processing window.
bool gnssBelongsToNoLidarWindow(double gnss_time, double gnss_local_time_offset,
                                const NoLidarWindow &window);

// LiDAR mode publishes at the scan end; GNSS+IMU-only mode publishes at the
// latest state propagation time because no LiDAR timestamp exists.
double mappingOutputTime(bool no_lidar, double state_time,
                         double lidar_scan_end_time);
