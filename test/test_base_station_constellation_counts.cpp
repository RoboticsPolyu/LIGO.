#include "BaseStationData.h"

#include <gtest/gtest.h>
#include <gnss_comm/gnss_constant.hpp>
#include <gnss_comm/gnss_utility.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>

namespace
{
double epochSeconds(double year, double month, double day, double hour,
                    double minute, double second)
{
  const double epoch[6] = {year, month, day, hour, minute, second};
  return gnss_comm::time2sec(gnss_comm::epoch2time(epoch));
}

struct ConstellationCounts
{
  std::set<uint32_t> satellites;
  std::set<uint32_t> primary_satellites;
  std::set<uint32_t> secondary_satellites;
  std::set<uint32_t> l5_satellites;
};

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
}  // namespace

TEST(BaseStationConstellationCounts, PrintsFirstEpochCountsFromHkscZip)
{
  const std::string archive =
      std::string(ROOT_DIR) +
      "Data/2026-Aug-02_001813/r3_1s_1h_hksc080i.2025.zip";

  BaseStationData data;
  std::string error;
  ASSERT_TRUE(data.load(archive, error)) << error;

  const double timestamp = epochSeconds(2025, 3, 21, 8, 0, 0.0);
  const auto *observations = data.epoch(timestamp, 0.01);
  ASSERT_NE(observations, nullptr);

  std::map<uint32_t, ConstellationCounts> counts;
  for (const BaseCarrierObservation &observation : *observations)
  {
    const uint32_t system = gnss_comm::satsys(observation.satellite, nullptr);
    ConstellationCounts &system_counts = counts[system];
    system_counts.satellites.insert(observation.satellite);
    if (std::abs(observation.frequency_hz - FREQ5) < 1.0)
      system_counts.l5_satellites.insert(observation.satellite);
    else if (std::abs(observation.frequency_hz - FREQ2) < 1.0 ||
             std::abs(observation.frequency_hz - FREQ7) < 1.0)
      system_counts.secondary_satellites.insert(observation.satellite);
    else
      system_counts.primary_satellites.insert(observation.satellite);
  }

  std::cout << "\nBase observation archive: " << archive << '\n'
            << "Epoch (GPST): 2025-03-21 08:00:00\n"
            << "Base ECEF XYZ (m): " << std::fixed << std::setprecision(4)
            << data.ecefPosition().transpose() << "\n\n"
            << std::left << std::setw(14) << "Constellation"
            << std::right << std::setw(12) << "Satellites"
            << std::setw(12) << "Primary"
            << std::setw(12) << "Secondary"
            << std::setw(12) << "L5/E5a" << '\n';

  size_t total_satellites = 0;
  for (const uint32_t system : {SYS_GPS, SYS_GLO, SYS_GAL, SYS_BDS})
  {
    const ConstellationCounts &system_counts = counts[system];
    std::cout << std::left << std::setw(14) << constellationName(system)
              << std::right << std::setw(12) << system_counts.satellites.size()
              << std::setw(12) << system_counts.primary_satellites.size()
              << std::setw(12) << system_counts.secondary_satellites.size()
              << std::setw(12) << system_counts.l5_satellites.size() << '\n';
    total_satellites += system_counts.satellites.size();
  }
  std::cout << std::left << std::setw(14) << "Total"
            << std::right << std::setw(12) << total_satellites << "\n\n";

  EXPECT_GT(counts[SYS_GPS].satellites.size(), 0u);
  EXPECT_GT(counts[SYS_GLO].satellites.size(), 0u);
  EXPECT_GT(counts[SYS_GAL].satellites.size(), 0u);
  EXPECT_GT(counts[SYS_BDS].satellites.size(), 0u);
  EXPECT_GT(total_satellites, 20u);
}

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
