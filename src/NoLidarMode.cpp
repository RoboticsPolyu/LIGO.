#include "NoLidarMode.h"

#include <cmath>

bool makeNoLidarWindow(double first_imu_time, double latest_imu_time,
                       double window_duration, NoLidarWindow &window)
{
    if (!std::isfinite(first_imu_time) || !std::isfinite(latest_imu_time) ||
        !std::isfinite(window_duration) || window_duration <= 0.0)
    {
        return false;
    }

    const double end_time = first_imu_time + window_duration;
    if (latest_imu_time < end_time)
    {
        return false;
    }

    window.begin_time = first_imu_time;
    window.end_time = end_time;
    return true;
}

bool gnssBelongsToNoLidarWindow(double gnss_time, double gnss_local_time_offset,
                                const NoLidarWindow &window)
{
    if (!std::isfinite(gnss_time) || !std::isfinite(gnss_local_time_offset))
    {
        return false;
    }
    return gnss_time - gnss_local_time_offset <= window.end_time;
}

double mappingOutputTime(bool no_lidar, double state_time,
                         double lidar_scan_end_time)
{
    return no_lidar ? state_time : lidar_scan_end_time;
}
