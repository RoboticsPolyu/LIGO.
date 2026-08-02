#include "BaseStationData.h"
#include "RtkSignalUtils.h"

#include <gtest/gtest.h>
#include <gnss_comm/GnssMeasMsg.h>
#include <gnss_comm/gnss_constant.hpp>
#include <gnss_comm/gnss_ros.hpp>
#include <gnss_comm/gnss_utility.hpp>
#include <rosbag/bag.h>
#include <rosbag/view.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace
{
const char *constellationName(uint32_t system)
{
  switch (system)
  {
    case SYS_GPS: return "GPS";
    case SYS_GLO: return "GLONASS";
    case SYS_GAL: return "Galileo";
    case SYS_BDS: return "BeiDou";
    default: return "Unknown";
  }
}

double messageTimestamp(const gnss_comm::GnssMeasMsg &message)
{
  if (message.meas.empty()) return 0.0;
  const auto &time = message.meas.front().time;
  return gnss_comm::time2sec(gnss_comm::gpst2time(time.week, time.tow));
}
}  // namespace

TEST(RtkBagCompatibility, ReportsTimeAndSignalCompatibility)
{
  const std::string bag_path =
      std::string(ROOT_DIR) + "Dataset/2025-03-21-16-46-47.bag";
  const std::string base_path =
      std::string(ROOT_DIR) + "Data/2026-Aug-02_001813";

  BaseStationData base;
  std::string error;
  ASSERT_TRUE(base.load(base_path, error)) << error;

  rosbag::Bag bag;
  ASSERT_NO_THROW(bag.open(bag_path, rosbag::bagmode::Read));
  rosbag::View view(
      bag, rosbag::TopicQuery(std::string("/ublox_driver/range_meas")));
  ASSERT_GT(view.size(), 0u);

  double rover_first = 0.0;
  double rover_last = 0.0;
  size_t rover_epoch_count = 0;
  size_t rover_l5_observation_count = 0;
  size_t rover_secondary_observation_count = 0;
  size_t synchronized_epoch_count = 0;
  size_t matched_carrier_count = 0;
  size_t potential_dd_factor_count = 0;
  size_t rover_loss_of_lock_count = 0;
  double maximum_synchronized_offset = 0.0;
  std::map<RtkSignalBand, size_t> matched_by_band;
  std::map<uint32_t, size_t> loss_of_lock_by_satellite;
  std::map<uint8_t, size_t> lli_value_counts;
  std::map<uint32_t, std::set<uint32_t>> first_epoch_satellites;
  for (const rosbag::MessageInstance &instance : view)
  {
    const auto message = instance.instantiate<gnss_comm::GnssMeasMsg>();
    ASSERT_TRUE(static_cast<bool>(message));
    const double timestamp = messageTimestamp(*message);
    if (timestamp <= 0.0) continue;
    if (rover_first == 0.0)
    {
      rover_first = timestamp;
      for (const auto &observation : message->meas)
      {
        const uint32_t system = gnss_comm::satsys(observation.sat, nullptr);
        first_epoch_satellites[system].insert(observation.sat);
      }
    }
    rover_last = timestamp;
    ++rover_epoch_count;
    for (const auto &observation : message->meas)
      for (const double frequency : observation.freqs)
      {
        if (std::abs(frequency - FREQ5) < 1.0)
          ++rover_l5_observation_count;
        if (std::abs(frequency - FREQ2) < 1.0 ||
            std::abs(frequency - FREQ7) < 1.0)
          ++rover_secondary_observation_count;
      }

    double matched_base_timestamp = 0.0;
    const auto *base_epoch = base.epoch(timestamp, 0.05,
                                        &matched_base_timestamp);
    if (!base_epoch) continue;
    ++synchronized_epoch_count;
    maximum_synchronized_offset = std::max(
        maximum_synchronized_offset,
        std::abs(matched_base_timestamp - timestamp));
    std::map<std::pair<uint32_t, RtkSignalBand>, size_t> matched_groups;
    const std::vector<gnss_comm::ObsPtr> rover_observations =
        gnss_comm::msg2meas(message);
    for (const gnss_comm::ObsPtr &observation : rover_observations)
    {
      for (size_t index = 0; index < observation->freqs.size(); ++index)
      {
        if (index >= observation->cp.size() || observation->cp[index] == 0.0)
          continue;
        const double frequency = observation->freqs[index];
        const RtkSignalBand band =
            classifyRtkSignal(observation, frequency, true, true);
        if (band == RtkSignalBand::Unsupported) continue;
        if (index < observation->LLI.size())
        {
          ++lli_value_counts[observation->LLI[index]];
          if (rtkLossOfLock(observation->LLI[index]))
          {
            ++rover_loss_of_lock_count;
            ++loss_of_lock_by_satellite[observation->sat];
          }
        }
        if (!findMatchingBaseCarrier(*base_epoch, observation->sat, frequency))
          continue;
        ++matched_carrier_count;
        ++matched_by_band[band];
        ++matched_groups[{gnss_comm::satsys(observation->sat, nullptr), band}];
      }
    }
    for (const auto &group : matched_groups)
      if (group.second >= 2) potential_dd_factor_count += group.second - 1;
  }
  bag.close();

  ASSERT_GT(rover_first, 0.0);
  ASSERT_GE(rover_last, rover_first);
  const bool overlaps =
      rover_last >= base.firstTimestamp() && rover_first <= base.lastTimestamp();
  const double gap_seconds = base.firstTimestamp() - rover_last;

  std::cout << "\nRTK rosbag/base compatibility\n"
            << "Rover epochs: " << rover_epoch_count << '\n'
            << std::fixed << std::setprecision(3)
            << "Rover GPST range: " << rover_first << " .. " << rover_last << '\n'
            << "Base  GPST range: " << base.firstTimestamp() << " .. "
            << base.lastTimestamp() << '\n'
            << "Base RINEX time system: " << base.timeSystem() << '\n'
            << "Time overlap: " << (overlaps ? "YES" : "NO") << '\n';
  if (!overlaps)
    std::cout << "Gap from rover end to base start: " << gap_seconds
              << " s (" << gap_seconds / 3600.0 << " h)\n";
  std::cout << "Rover L5 observations in bag: "
            << rover_l5_observation_count << '\n'
            << "Rover secondary observations in bag: "
            << rover_secondary_observation_count << '\n'
            << "Synchronized rover/base epochs (0.05 s): "
            << synchronized_epoch_count << '\n'
            << "Maximum synchronized time offset: "
            << maximum_synchronized_offset << " s\n"
            << "Matched supported carrier observations: "
            << matched_carrier_count << '\n'
            << "  primary=" << matched_by_band[RtkSignalBand::Primary]
            << " secondary=" << matched_by_band[RtkSignalBand::Secondary]
            << " l5=" << matched_by_band[RtkSignalBand::L5] << '\n'
            << "Matched rover carriers with LLI loss/half-cycle flag: "
            << rover_loss_of_lock_count
            << " (LLI 1=" << lli_value_counts[1]
            << ", 2=" << lli_value_counts[2]
            << ", 3=" << lli_value_counts[3]
            << ", Galileo BOC bit 4=" << lli_value_counts[4] << ")\n"
            << "Potential DD carrier factors: " << potential_dd_factor_count
            << "\nFirst rover epoch satellites:\n";
  for (const uint32_t system : {SYS_GPS, SYS_GLO, SYS_GAL, SYS_BDS})
    std::cout << "  " << std::setw(8) << std::left << constellationName(system)
              << first_epoch_satellites[system].size() << '\n';
  std::cout << "RTK DD testable with these files: "
            << (overlaps ? "YES" : "NO") << "\n\n";

  std::cout << "Satellites with the most synchronized LLI rejections:\n";
  std::vector<std::pair<size_t, uint32_t>> sorted_lli;
  for (const auto &entry : loss_of_lock_by_satellite)
    sorted_lli.emplace_back(entry.second, entry.first);
  std::sort(sorted_lli.rbegin(), sorted_lli.rend());
  for (size_t index = 0; index < std::min<size_t>(10, sorted_lli.size()); ++index)
    std::cout << "  " << gnss_comm::sat2str(sorted_lli[index].second)
              << ": " << sorted_lli[index].first << '\n';
  std::cout << '\n';

  // The actual integration inputs must overlap and provide at least one
  // same-constellation reference/subject carrier pair.
  EXPECT_EQ(rover_epoch_count, 4243u);
  EXPECT_TRUE(overlaps);
  EXPECT_GT(synchronized_epoch_count, 0u);
  EXPECT_LT(maximum_synchronized_offset, 0.01);
  EXPECT_GT(matched_carrier_count, 0u);
  EXPECT_GT(potential_dd_factor_count, 0u);
  EXPECT_GT(rover_secondary_observation_count, 0u);
}

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
