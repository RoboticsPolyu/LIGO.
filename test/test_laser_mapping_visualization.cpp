#include "LaserMappingVisualization.h"

#include <gtest/gtest.h>

namespace
{
GNSSProcess::RtkSatelliteDirection direction(uint32_t satellite, double x)
{
  return {satellite, Eigen::Vector3d(x, 0.0, 0.0)};
}
}  // namespace

TEST(LaserMappingVisualization, AppliesHighestCurrentProcessingStage)
{
  const std::vector<GNSSProcess::RtkSatelliteDirection> raw = {
      direction(1, 1.0), direction(2, 2.0), direction(3, 3.0)};
  const std::vector<GNSSProcess::RtkSatelliteDirection> valid = {
      direction(2, 20.0), direction(3, 30.0)};
  const std::vector<GNSSProcess::RtkSatelliteDirection> rtk = {
      direction(3, 300.0)};

  const auto staged = stageRtkSatellites(raw, valid, rtk);
  ASSERT_EQ(staged.size(), 3u);
  EXPECT_EQ(staged.at(1).stage, SatelliteMarkerStage::Raw);
  EXPECT_EQ(staged.at(2).stage, SatelliteMarkerStage::Valid);
  EXPECT_EQ(staged.at(3).stage, SatelliteMarkerStage::RtkCandidate);
  EXPECT_DOUBLE_EQ(staged.at(3).direction.local_unit_direction.x(), 300.0);
}

TEST(LaserMappingVisualization, RejectsStaleOrUnknownRtkCandidates)
{
  const std::vector<GNSSProcess::RtkSatelliteDirection> raw = {
      direction(4, 4.0)};
  const std::vector<GNSSProcess::RtkSatelliteDirection> valid;
  const std::vector<GNSSProcess::RtkSatelliteDirection> rtk = {
      direction(4, 40.0), direction(5, 50.0)};

  const auto staged = stageRtkSatellites(raw, valid, rtk);
  ASSERT_EQ(staged.size(), 1u);
  EXPECT_EQ(staged.at(4).stage, SatelliteMarkerStage::Raw);
  EXPECT_EQ(staged.count(5), 0u);
}

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
