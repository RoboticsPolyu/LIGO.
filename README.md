# LIGO 

**LIGO: Tightly Coupled LiDAR-Inertial-GNSS Odometry based on a Hierarchy Fusion Framework for Global Localization with Real-time Mapping**

Code, paper, video are coming soon......
instruction for using will be detailed soon

Our paper is published on [TRO](https://github.com/Joanna-HE/LIGO./blob/main/paper/LIGO_A_Tightly_Coupled_LiDAR-Inertial-GNSS_Odometry_Based_on_a_Hierarchy_Fusion_Framework_for_Global_Localization_With_Real-Time_Mapping.pdf)

Our datasets are uploaded on [Google Drive](https://drive.google.com/drive/folders/1hNwl8u8Pg-SqKh2N808XFixj6PjPf091?usp=sharing)

# Developers:
The codes of this repo are contributed by:
[Dongjiao He (贺东娇)](https://github.com/Joanna-HE)

# Properties

**LIGO is a multi-sensor fusion framework that maximizes the complementary properties of both LiDAR and GNSS systems**. This package achieves the following properties:

1. Competitive accuracy in trajectory estimation across large-scale scenarios.
2. Robustness to malfunctions of either GNSS or LiDAR sensors, enabling seamless handling of added or lost sensor signals during operation.
3. High-output-frequency odometry.
4. Capability of providing globally referenced pose estimations in both indoor and outdoor environments, suitable for ground vehicles and uncrewed aerial vehicles (UAVs).
5. No requirement for GNSS observations to be obtained exactly at the beginning or end time of LiDAR scans.
6. Robustness to large outliers and high noise levels in GNSS observations.

# Hardware setups for self-collected datasets

## Setup

<div align="left">
    <div align="left">
        <img src="https://github.com/Joanna-HE/LIGO./blob/main/image/hardware.jpg" width = 30% >
    </div>
</div>

Platform: DJI Matrice 300   
Onboard computer: DJI Manifold 2-c 256G, CPU: Intel i7-8550U  
LiDAR: Livox Mid360 and Livox Avia  
IMU: Built-in IMU of Livox LiDAR  
GNSS receiver: u-blox C099-F9P-2  
GNSS antenna: B4QA4GGGB  

## Recording rates

LiDAR: 10Hz  
IMU: 200Hz  
GNSS: 10Hz  
RTK: 10Hz  

## Recording software:

Operating system: Ubuntu 20.04  
IMU and LiDAR driver: Livox driver  
GNSS driver: [ublox driver](https://github.com/Joanna-HE/ublox_driver)  

## ROS topics recorded

IMU: /livox/imu  
LiDAR: /livox/lidar  
RAW GNSS: /ublox_driver/range_meas  
GNSS EPHEM: /ublox_driver/ephem and /ublox_driver/glo_ephem  
IONO PARAMETER: /ublox_driver/iono_params  
Onboard pos solution of ublox: /ublox_driver/receiver_pvt and /ublox_driver/receiver_lla  
PPS time info: /ublox_driver/time_pulse_info  

## Time synchronization

PPS: Livox LiDARs can receive pps and gprmc given by the GNSS receiver  
The time difference between LiDAR and IMU is zero, and between LiDAR and GNSS message is 18.0 s  

## RTK solution

Please follow the *Differential GNSS* section shown in [ublox driver](https://github.com/Joanna-HE/ublox_driver) to get the differential GNSS solution online or offline. The self-collected datasets get the online RTK solution which are saved in the topic '/ublox_driver/receiver_pvt', the value of the 'carr_soln' as 1 and 'diff_soln' as 2 indicates the fix RTK solution.

# Build

## Prerequisites

We test LIGO on ubuntu 20.04 with ROS noetic, and C++17 compiler & Eigen 3 & GTSAM 4 & opencv 4.2.0 & pcl 1.10  

### Install Boost using command sudo apt-get install libboost-all-dev

### Install [Livox Driver](https://github.com/Livox-SDK/livox_ros_driver)

### Install [gnss_comm](https://github.com/HKUST-Aerial-Robotics/gnss_comm) with its [instuction](https://github.com/HKUST-Aerial-Robotics/gnss_comm#2-build-gnss_comm-library)

## Make

### clone the code to catkin_ws workspace
```
cd ~/catkin_ws/src/
git clone https://github.com/Joanna-HE/LIGO..git
```
### compile the package
```
cd ~/catkin_ws/
source /PATH/TO/LIVOX_DRIVER/DEVEL/.setup.bash
source /PATH/TO/GNSS_COMM/DEVEL/.setup.bash
catkin_make
source ~/catkin_ws/devel/setup.bash
```

### Base-station carrier data (float/fixed RTK)

Set `gnss/use_double_differences: true` and `gnss/base_station_file` in the
selected YAML file to add base/rover double-difference carrier factors. The
path is resolved under `Data/` and may name a RINEX 3.x observation file, a
ZIP archive containing an observation file, or a directory of such files and
archives. RINEX 3.02 is the tested format. The configured HKSC directory can
be used directly without extracting its hourly ZIP files:

```
base_station_file: "2026-Aug-02_001813"
```

The reader obtains the surveyed ECEF antenna coordinates from `APPROX POSITION
XYZ`, the observation layout from multiline `SYS / # / OBS TYPES` records,
and carrier loss-of-lock from each phase observation's LLI field. Because the
RINEX SSI digit is signal strength rather than uncertainty, configure
`gnss/base_carrier_std_cycles` as the base carrier-phase standard deviation.

Set `gnss/use_secondary: true` to add GPS L2, Galileo E5b, and BeiDou B2.
Set `gnss/use_l5: true` to add GPS L5, Galileo E5a, and BeiDou B2a signals
when they are present at both receivers. The RINEX L1 and L5 observation types
remain separate measurements for the same epoch and satellite. Each band has
its own reference selection, cycle-slip arc, wavelength, and integer
ambiguity; eligible ambiguities from both bands are resolved jointly by
LAMBDA. GLONASS remains L1-only because this signal mapping has no GLONASS L5
equivalent.

When primary and secondary carrier observations are available for the same
satellite, the graph also receives supplementary wide-lane and narrow-lane DD
factors. These combinations have separate ambiguity arcs, propagate the noise
of both source carriers, and are reset when either source carrier loses lock.

Synchronized primary-band code observations also add double-differenced
pseudorange factors. Their standard deviation in metres is configured with
`gnss/double_difference_pseudorange_sigma` (default: `2.0`).

Epoch times must use the same GPST seconds as the rover observations in the
bag. RINEX UTC and BDT headers are converted to GPST. Note that the supplied
HKSC archives contain observations from 21 March 2025 at 08:00--10:00 GPST;
the enclosing 2026 directory name is not the observation date. Ambiguities are persistent
real-valued graph states and reset after a
sampling gap or loss-of-lock. With `gnss/enable_integer_fixing: true`, stable
same-frequency ambiguities are resolved with LAMBDA after the configured lock
period. A covariance gate, ratio test, and partial ambiguity selection guard
the fix. Accepted integers are added as tight fix-and-hold graph constraints;
when validation fails, optimization continues with the float solution.

The RINEX input tests are in `test/test_base_station_data.cpp`, and the RTK
core regression tests are kept separately in `test/test_rtk_core.cpp`. Run
the focused targets with:

```
catkin_make run_tests_ligo_gtest_test_base_station_data
catkin_make run_tests_ligo_gtest_test_base_station_constellation_counts
catkin_make run_tests_ligo_gtest_test_rtk_bag_compatibility
catkin_make run_tests_ligo_gtest_test_rtk_core
```

The constellation-count target reads the supplied `hksc080i` ZIP directly and
prints the number of GPS, GLONASS, Galileo, and BeiDou satellites at its first
epoch, including separate primary-band and L5/E5a counts.
The bag-compatibility target scans every rover GNSS epoch and reports time
overlap, matched carriers, and the number of constructible DD factors.

For a headless integration replay, use
`launch/mapping_avia_rtk_bag_test.launch`. The current bag publishes
`livox_ros_driver2/CustomMsg`, whereas this repository consumes
`livox_ros_driver/CustomMsg`; consequently the launch uses the experimental
GNSS/IMU-only branch. A driver-2 adapter is still required for a complete
LiDAR/GNSS factor-graph replay.

# Demo
**Performance on a sequence with severe LiDAR degeneracy**

<div align="center">
    <div align="center">
        <img src="https://github.com/Joanna-HE/LIGO/blob/main/image/Sample.png" width = 75% >
    </div>
</div>
