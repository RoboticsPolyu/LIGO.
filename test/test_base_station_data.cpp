#include "BaseStationData.h"

#include <gtest/gtest.h>
#include <gnss_comm/gnss_constant.hpp>
#include <gnss_comm/gnss_utility.hpp>

#include <algorithm>
#include <cmath>


namespace
{
double epochSeconds(double year, double month, double day, double hour,
                    double minute, double second)
{
  const double epoch[6] = {year, month, day, hour, minute, second};
  return gnss_comm::time2sec(gnss_comm::epoch2time(epoch));
}
}  // namespace

TEST(BaseStationData, LoadsRinex302AndMatchesPrimarySecondaryL5)
{
  BaseStationData data;
  std::string error;
  ASSERT_TRUE(data.load(std::string(ROOT_DIR) +
                            "test/data/base_station_rinex_302.obs",
                        error))
      << error;
  EXPECT_NEAR(data.ecefPosition().x(), -2414266.9200, 1e-6);
  EXPECT_EQ(data.timeSystem(), "GPS");
  EXPECT_EQ(data.epochCount(), 1u);

  const double timestamp = epochSeconds(2025, 3, 21, 16, 0, 0.0);
  const auto *exact = data.epoch(timestamp, 0.01);
  ASSERT_NE(exact, nullptr);
  EXPECT_EQ(exact->size(), 6u);
  EXPECT_EQ(std::count_if(exact->begin(), exact->end(), [](const auto &observation) {
              return std::abs(observation.frequency_hz - FREQ2) < 1.0;
            }), 2);
  EXPECT_EQ(std::count_if(exact->begin(), exact->end(), [](const auto &observation) {
              return std::abs(observation.frequency_hz - FREQ5) < 1.0;
            }), 2);
  const auto l5 = std::find_if(exact->begin(), exact->end(), [](const auto &observation) {
    return observation.satellite == gnss_comm::sat_no(SYS_GPS, 18) &&
           std::abs(observation.frequency_hz - FREQ5) < 1.0;
  });
  ASSERT_NE(l5, exact->end());
  EXPECT_NEAR(l5->carrier_cycles, 80743528.462, 1e-6);
  EXPECT_EQ(l5->loss_of_lock, 0);
  EXPECT_DOUBLE_EQ(l5->carrier_std_cycles, 0.01);

  double matched_timestamp = 0.0;
  EXPECT_NE(data.epoch(timestamp + 0.004, 0.01, &matched_timestamp), nullptr);
  EXPECT_DOUBLE_EQ(matched_timestamp, timestamp);
  EXPECT_EQ(data.epoch(timestamp + 1.0, 0.01, &matched_timestamp), nullptr);
  EXPECT_DOUBLE_EQ(matched_timestamp, timestamp);

}

TEST(BaseStationData, LoadsActualHkscZipDirectory)
{
  BaseStationData data;
  std::string error;
  ASSERT_TRUE(data.load(std::string(ROOT_DIR) + "Data/2026-Aug-02_001813",
                        error))
      << error;
  EXPECT_NEAR(data.ecefPosition().x(), -2414266.9200, 1e-4);
  EXPECT_EQ(data.timeSystem(), "GPS");
  EXPECT_EQ(data.epochCount(), 7200u);
  const auto *epoch =
      data.epoch(epochSeconds(2025, 3, 21, 8, 0, 0.0), 0.01);
  ASSERT_NE(epoch, nullptr);
  EXPECT_GT(epoch->size(), 20u);
  EXPECT_GT(std::count_if(epoch->begin(), epoch->end(), [](const auto &observation) {
              return std::abs(observation.frequency_hz - FREQ5) < 1.0;
            }), 5);
  EXPECT_GT(std::count_if(epoch->begin(), epoch->end(), [](const auto &observation) {
              return std::abs(observation.frequency_hz - FREQ2) < 1.0 ||
                     std::abs(observation.frequency_hz - FREQ7) < 1.0;
            }), 5);
  // Galileo BOC tracking sets RINEX LLI bit 2. It must remain available and
  // must not be interpreted as a cycle slip by the RTK arc validator.
  EXPECT_GT(std::count_if(epoch->begin(), epoch->end(), [](const auto &observation) {
              return observation.loss_of_lock == 4;
            }), 0);

  // The directory loader must aggregate the following hourly archive too.
  const auto *second_hour =
      data.epoch(epochSeconds(2025, 3, 21, 9, 0, 0.0), 0.01);
  ASSERT_NE(second_hour, nullptr);
  EXPECT_GT(second_hour->size(), 20u);
}

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
