#include <gtest/gtest.h>

#include <limits>

#include "NoLidarMode.h"

TEST(NoLidarMode, WaitsForACompleteImuWindow)
{
    NoLidarWindow window;
    EXPECT_FALSE(makeNoLidarWindow(10.0, 10.09, 0.1, window));
    ASSERT_TRUE(makeNoLidarWindow(10.0, 10.1, 0.1, window));
    EXPECT_DOUBLE_EQ(window.begin_time, 10.0);
    EXPECT_DOUBLE_EQ(window.end_time, 10.1);
}

TEST(NoLidarMode, SelectsGnssEpochsInTheLocalTimeWindow)
{
    const NoLidarWindow window{100.0, 100.1};
    EXPECT_TRUE(gnssBelongsToNoLidarWindow(118.1, 18.0, window));
    EXPECT_FALSE(gnssBelongsToNoLidarWindow(118.1001, 18.0, window));
}

TEST(NoLidarMode, RejectsInvalidWindowConfiguration)
{
    NoLidarWindow window;
    EXPECT_FALSE(makeNoLidarWindow(10.0, 11.0, 0.0, window));
    EXPECT_FALSE(makeNoLidarWindow(
        std::numeric_limits<double>::quiet_NaN(), 11.0, 0.1, window));
}

TEST(NoLidarMode, UsesStateTimeWhenNoLidarTimestampExists)
{
    EXPECT_DOUBLE_EQ(mappingOutputTime(true, 42.5, 0.0), 42.5);
    EXPECT_DOUBLE_EQ(mappingOutputTime(false, 42.5, 42.4), 42.4);
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
