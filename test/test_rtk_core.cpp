#include "LambdaAmbiguityResolver.h"
#include "RtkSignalUtils.h"
#include "RtkSolutionStatus.h"
#include "gnss_factor/gnss_double_diff_carrier_factor.hpp"
#include "gnss_factor/gnss_double_diff_pseudorange_factor.hpp"
#include "gnss_factor/wide_lane_integer_factor.hpp"

#include <gtest/gtest.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/linear/NoiseModel.h>

#include <limits>
#include <random>

namespace
{
// Millimetre-scale translation perturbations avoid subtractive cancellation
// when differencing GNSS ranges of roughly 20,000 km. The same magnitude is
// small enough for the central rotational derivative's second-order error.
constexpr double kDifferenceStep = 1.0e-3;

double doubleDifferenceGeometry(const Eigen::Vector3d &rover,
                                const Eigen::Vector3d &base,
                                const Eigen::Vector3d &satellite,
                                const Eigen::Vector3d &reference)
{
  return (satellite - rover).norm() - (satellite - base).norm() -
         (reference - rover).norm() + (reference - base).norm();
}

template <typename Function>
Eigen::MatrixXd vectorJacobian(const Eigen::VectorXd &value, Function function)
{
  const Eigen::VectorXd nominal = function(value);
  Eigen::MatrixXd jacobian(nominal.size(), value.size());
  for (int column = 0; column < value.size(); ++column)
  {
    Eigen::VectorXd plus = value;
    Eigen::VectorXd minus = value;
    plus[column] += kDifferenceStep;
    minus[column] -= kDifferenceStep;
    jacobian.col(column) =
        (function(plus) - function(minus)) / (2.0 * kDifferenceStep);
  }
  return jacobian;
}

template <typename Function>
Eigen::MatrixXd rotationJacobian(const gtsam::Rot3 &rotation, Function function)
{
  const Eigen::VectorXd nominal = function(rotation);
  Eigen::MatrixXd jacobian(nominal.size(), 3);
  for (int column = 0; column < 3; ++column)
  {
    gtsam::Vector3 delta = gtsam::Vector3::Zero();
    delta[column] = kDifferenceStep;
    jacobian.col(column) =
        (function(rotation.retract(delta)) - function(rotation.retract(-delta))) /
        (2.0 * kDifferenceStep);
  }
  return jacobian;
}

void expectMatrixNear(const Eigen::MatrixXd &actual,
                      const Eigen::MatrixXd &expected, double tolerance)
{
  ASSERT_EQ(actual.rows(), expected.rows());
  ASSERT_EQ(actual.cols(), expected.cols());
  for (int row = 0; row < actual.rows(); ++row)
    for (int column = 0; column < actual.cols(); ++column)
      EXPECT_NEAR(actual(row, column), expected(row, column), tolerance)
          << "at (" << row << ", " << column << ")";
}

struct FactorFixture
{
  Eigen::Vector3d base{120.0, -80.0, 35.0};
  Eigen::Vector3d satellite{20200000.0, 14000000.0, 21700000.0};
  Eigen::Vector3d reference{-15600000.0, 21100000.0, 19300000.0};
  Eigen::Vector3d anchor{135.0, -62.0, 42.0};
  Eigen::Vector3d antenna{0.42, -0.17, 0.83};
  gtsam::Rot3 rotation = gtsam::Rot3::RzRyRx(0.11, -0.07, 0.23);
  gtsam::Vector6 state = (gtsam::Vector6() << 8.0, -4.0, 2.5,
                            1.2, -0.3, 0.1).finished();
  gtsam::Vector1 ambiguity = gtsam::Vector1(3.75);
  gtsam::SharedNoiseModel noise =
      gtsam::noiseModel::Isotropic::Sigma(1, 0.02);
};

Eigen::Vector2d bruteForceBest(const Eigen::Vector2d &floating,
                               const Eigen::Matrix2d &covariance,
                               double *best_norm)
{
  const Eigen::Matrix2d information = covariance.inverse();
  Eigen::Vector2d best = Eigen::Vector2d::Zero();
  *best_norm = std::numeric_limits<double>::infinity();
  for (int first = -10; first <= 10; ++first)
    for (int second = -10; second <= 10; ++second)
    {
      const Eigen::Vector2d integer(first, second);
      const Eigen::Vector2d residual = floating - integer;
      const double norm = residual.transpose() * information * residual;
      if (norm < *best_norm)
      {
        *best_norm = norm;
        best = integer;
      }
    }
  return best;
}
}  // namespace

TEST(DoubleDiffCarrierFactor, LocalFactorHasCorrectResidualAndJacobians)
{
  const FactorFixture data;
  const Eigen::Vector3d rover =
      data.anchor + data.rotation.matrix() * (data.state.head<3>() + data.antenna);
  const double measured = doubleDifferenceGeometry(
                              rover, data.base, data.satellite, data.reference) +
                          data.ambiguity[0];
  const ligo::DoubleDiffCarrierFactor factor(
      gtsam::Symbol('e', 0), gtsam::Symbol('p', 0), gtsam::Symbol('x', 0),
      gtsam::Symbol('n', 0), data.antenna, data.base, data.satellite,
      data.reference, measured, data.noise);

  gtsam::Matrix H1, H2, H3, H4;
  const gtsam::Vector residual = factor.evaluateError(
      data.anchor, data.rotation, data.state, data.ambiguity,
      &H1, &H2, &H3, &H4);
  EXPECT_NEAR(residual[0], 0.0, 1.0e-9);

  const auto evaluate_anchor = [&](const Eigen::VectorXd &anchor) {
    return factor.evaluateError(anchor, data.rotation, data.state,
                                data.ambiguity, nullptr, nullptr, nullptr, nullptr);
  };
  const auto evaluate_rotation = [&](const gtsam::Rot3 &rotation) {
    return factor.evaluateError(data.anchor, rotation, data.state,
                                data.ambiguity, nullptr, nullptr, nullptr, nullptr);
  };
  const auto evaluate_state = [&](const Eigen::VectorXd &state) {
    return factor.evaluateError(data.anchor, data.rotation, state,
                                data.ambiguity, nullptr, nullptr, nullptr, nullptr);
  };
  const auto evaluate_ambiguity = [&](const Eigen::VectorXd &ambiguity) {
    return factor.evaluateError(data.anchor, data.rotation, data.state,
                                ambiguity, nullptr, nullptr, nullptr, nullptr);
  };
  expectMatrixNear(H1, vectorJacobian(data.anchor, evaluate_anchor), 2.0e-5);
  expectMatrixNear(H2, rotationJacobian(data.rotation, evaluate_rotation), 2.0e-5);
  expectMatrixNear(H3, vectorJacobian(data.state, evaluate_state), 2.0e-5);
  expectMatrixNear(H4, vectorJacobian(data.ambiguity, evaluate_ambiguity), 2.0e-5);
}

TEST(DoubleDiffCarrierFactor, EcefFactorHasCorrectResidualAndJacobians)
{
  const FactorFixture data;
  gtsam::Vector12 state = gtsam::Vector12::Zero();
  state.head<3>() = data.anchor;
  const Eigen::Vector3d rover = data.anchor + data.rotation * data.antenna;
  const double measured = doubleDifferenceGeometry(
                              rover, data.base, data.satellite, data.reference) +
                          data.ambiguity[0];
  const ligo::DoubleDiffCarrierFactorNolidar factor(
      gtsam::Symbol('r', 0), gtsam::Symbol('f', 0), gtsam::Symbol('n', 0),
      data.antenna, data.base, data.satellite, data.reference, measured,
      data.noise);

  gtsam::Matrix H1, H2, H3;
  const gtsam::Vector residual = factor.evaluateError(
      data.rotation, state, data.ambiguity, &H1, &H2, &H3);
  EXPECT_NEAR(residual[0], 0.0, 1.0e-9);

  const auto evaluate_rotation = [&](const gtsam::Rot3 &rotation) {
    return factor.evaluateError(rotation, state, data.ambiguity, nullptr,
                                nullptr, nullptr);
  };
  const auto evaluate_state = [&](const Eigen::VectorXd &candidate) {
    return factor.evaluateError(data.rotation, candidate, data.ambiguity,
                                nullptr, nullptr, nullptr);
  };
  const auto evaluate_ambiguity = [&](const Eigen::VectorXd &ambiguity) {
    return factor.evaluateError(data.rotation, state, ambiguity, nullptr,
                                nullptr, nullptr);
  };
  expectMatrixNear(H1, rotationJacobian(data.rotation, evaluate_rotation), 2.0e-5);
  expectMatrixNear(H2, vectorJacobian(state, evaluate_state), 2.0e-5);
  expectMatrixNear(H3, vectorJacobian(data.ambiguity, evaluate_ambiguity), 2.0e-5);
}

TEST(DoubleDiffPseudorangeFactor, LocalFactorHasCorrectResidualAndJacobians)
{
  const FactorFixture data;
  const Eigen::Vector3d rover = data.anchor + data.rotation.matrix() *
      (data.state.head<3>() + data.antenna);
  const double measured = doubleDifferenceGeometry(
      rover, data.base, data.satellite, data.reference);
  const ligo::DoubleDiffPseudorangeFactor factor(
      gtsam::Symbol('e', 0), gtsam::Symbol('p', 0), gtsam::Symbol('x', 0),
      data.antenna, data.base, data.satellite, data.reference, measured,
      data.noise);

  gtsam::Matrix H1, H2, H3;
  EXPECT_NEAR(factor.evaluateError(data.anchor, data.rotation, data.state,
                                   &H1, &H2, &H3)[0], 0.0, 1.0e-9);
  const auto evaluate_anchor = [&](const Eigen::VectorXd &anchor) {
    return factor.evaluateError(anchor, data.rotation, data.state,
                                nullptr, nullptr, nullptr);
  };
  const auto evaluate_rotation = [&](const gtsam::Rot3 &rotation) {
    return factor.evaluateError(data.anchor, rotation, data.state,
                                nullptr, nullptr, nullptr);
  };
  const auto evaluate_state = [&](const Eigen::VectorXd &state) {
    return factor.evaluateError(data.anchor, data.rotation, state,
                                nullptr, nullptr, nullptr);
  };
  expectMatrixNear(H1, vectorJacobian(data.anchor, evaluate_anchor), 2.0e-5);
  expectMatrixNear(H2, rotationJacobian(data.rotation, evaluate_rotation), 2.0e-5);
  expectMatrixNear(H3, vectorJacobian(data.state, evaluate_state), 2.0e-5);
}

TEST(DoubleDiffPseudorangeFactor, EcefFactorHasCorrectResidualAndJacobians)
{
  const FactorFixture data;
  gtsam::Vector12 state = gtsam::Vector12::Zero();
  state.head<3>() = data.anchor;
  const Eigen::Vector3d rover = data.anchor + data.rotation * data.antenna;
  const double measured = doubleDifferenceGeometry(
      rover, data.base, data.satellite, data.reference);
  const ligo::DoubleDiffPseudorangeFactorNolidar factor(
      gtsam::Symbol('r', 0), gtsam::Symbol('f', 0), data.antenna, data.base,
      data.satellite, data.reference, measured, data.noise);

  gtsam::Matrix H1, H2;
  EXPECT_NEAR(factor.evaluateError(data.rotation, state, &H1, &H2)[0],
              0.0, 1.0e-9);
  const auto evaluate_rotation = [&](const gtsam::Rot3 &rotation) {
    return factor.evaluateError(rotation, state, nullptr, nullptr);
  };
  const auto evaluate_state = [&](const Eigen::VectorXd &candidate) {
    return factor.evaluateError(data.rotation, candidate, nullptr, nullptr);
  };
  expectMatrixNear(H1, rotationJacobian(data.rotation, evaluate_rotation), 2.0e-5);
  expectMatrixNear(H2, vectorJacobian(state, evaluate_state), 2.0e-5);
}

TEST(RtkSignalSelection, DistinguishesAndMatchesAllEnabledBands)
{
  auto observation = std::make_shared<gnss_comm::Obs>();
  observation->sat = gnss_comm::sat_no(SYS_GPS, 3);
  observation->freqs = {FREQ1, FREQ5, FREQ2};
  EXPECT_EQ(classifyRtkSignal(observation, FREQ1, true, true),
            RtkSignalBand::Primary);
  EXPECT_EQ(classifyRtkSignal(observation, FREQ5, true, true),
            RtkSignalBand::L5);
  EXPECT_EQ(classifyRtkSignal(observation, FREQ5, false, true),
            RtkSignalBand::Unsupported);
  EXPECT_EQ(classifyRtkSignal(observation, FREQ2, true, true),
            RtkSignalBand::Secondary);
  EXPECT_EQ(classifyRtkSignal(observation, FREQ2, true, false),
            RtkSignalBand::Unsupported);

  const std::vector<BaseCarrierObservation> base_epoch = {
      {observation->sat, FREQ1, 100.0, 0.01, 0},
      {observation->sat, FREQ5, 200.0, 0.02, 0}};
  const BaseCarrierObservation *l1 =
      findMatchingBaseCarrier(base_epoch, observation->sat, FREQ1);
  const BaseCarrierObservation *l5 =
      findMatchingBaseCarrier(base_epoch, observation->sat, FREQ5);
  ASSERT_NE(l1, nullptr);
  ASSERT_NE(l5, nullptr);
  EXPECT_DOUBLE_EQ(l1->carrier_cycles, 100.0);
  EXPECT_DOUBLE_EQ(l5->carrier_cycles, 200.0);
  EXPECT_EQ(findMatchingBaseCarrier(base_epoch,
                                    gnss_comm::sat_no(SYS_GPS, 8), FREQ5),
            nullptr);
}

TEST(RtkSignalSelection, RejectsL5ForUnsupportedConstellation)
{
  auto observation = std::make_shared<gnss_comm::Obs>();
  observation->sat = gnss_comm::sat_no(SYS_GLO, 5);
  observation->freqs = {FREQ1_GLO, FREQ5};
  EXPECT_EQ(classifyRtkSignal(observation, FREQ5, true, true),
            RtkSignalBand::Unsupported);
}

TEST(RtkSignalSelection, FormsWideAndNarrowLaneCombinations)
{
  const double l1_phase_m = 12.5;
  const double l2_phase_m = 11.0;
  const double l1_variance = 0.0004;
  const double l2_variance = 0.0009;
  const RtkLaneCombination wide = rtkLaneCombination(
      l1_phase_m, l1_variance, FREQ1, l2_phase_m, l2_variance, FREQ2, true);
  const RtkLaneCombination narrow = rtkLaneCombination(
      l1_phase_m, l1_variance, FREQ1, l2_phase_m, l2_variance, FREQ2, false);

  EXPECT_NEAR(wide.phase_m,
              (FREQ1 * l1_phase_m - FREQ2 * l2_phase_m) / (FREQ1 - FREQ2),
              1.0e-12);
  EXPECT_NEAR(wide.wavelength_m, LIGHT_SPEED / (FREQ1 - FREQ2), 1.0e-12);
  EXPECT_NEAR(narrow.phase_m,
              (FREQ1 * l1_phase_m + FREQ2 * l2_phase_m) / (FREQ1 + FREQ2),
              1.0e-12);
  EXPECT_NEAR(narrow.wavelength_m, LIGHT_SPEED / (FREQ1 + FREQ2), 1.0e-12);
  EXPECT_GT(wide.variance_m2, narrow.variance_m2);
}

TEST(RtkWideLane, BinaryConstraintUsesRawAmbiguitiesInCycles)
{
  const double lambda1 = 0.190293672798;
  const double lambda2 = 0.244210213425;
  ligo::WideLaneIntegerFactor factor(
      gtsam::Symbol('n', 1), gtsam::Symbol('n', 2), lambda1, lambda2, 7,
      gtsam::noiseModel::Isotropic::Sigma(1, 0.001));
  gtsam::Matrix H1, H2;
  const gtsam::Vector error = factor.evaluateError(
      gtsam::Vector1(12.25 * lambda1), gtsam::Vector1(5.25 * lambda2),
      &H1, &H2);
  ASSERT_EQ(error.size(), 1);
  EXPECT_NEAR(error[0], 0.0, 1.0e-12);
  EXPECT_NEAR(H1(0, 0), 1.0 / lambda1, 1.0e-12);
  EXPECT_NEAR(H2(0, 0), -1.0 / lambda2, 1.0e-12);
}

TEST(RtkWideLane, PropagatesCorrelatedRawCovarianceWithTQt)
{
  const double lambda1 = 0.19;
  const double lambda2 = 0.24;
  Eigen::Matrix4d raw_covariance;
  raw_covariance <<
      0.010, 0.003, 0.002, 0.001,
      0.003, 0.020, 0.004, 0.002,
      0.002, 0.004, 0.030, 0.006,
      0.001, 0.002, 0.006, 0.040;
  Eigen::Matrix<double, 2, 4> transform =
      Eigen::Matrix<double, 2, 4>::Zero();
  transform(0, 0) = 1.0 / lambda1;
  transform(0, 1) = -1.0 / lambda2;
  transform(1, 2) = 1.0 / lambda1;
  transform(1, 3) = -1.0 / lambda2;
  const Eigen::Matrix2d wide_covariance =
      transform * raw_covariance * transform.transpose();
  const double expected_variance0 =
      raw_covariance(0, 0) / (lambda1 * lambda1) +
      raw_covariance(1, 1) / (lambda2 * lambda2) -
      2.0 * raw_covariance(0, 1) / (lambda1 * lambda2);
  EXPECT_NEAR(wide_covariance(0, 0), expected_variance0, 1.0e-12);
  EXPECT_NE(wide_covariance(0, 1), 0.0);
  EXPECT_TRUE(wide_covariance.isApprox(wide_covariance.transpose()));
}

TEST(RtkSignalSelection, MapsGalileoAndBeiDouSecondaryBands)
{
  auto galileo = std::make_shared<gnss_comm::Obs>();
  galileo->sat = gnss_comm::sat_no(SYS_GAL, 15);
  galileo->freqs = {FREQ1, FREQ7, FREQ5};
  EXPECT_EQ(classifyRtkSignal(galileo, FREQ7, true, true),
            RtkSignalBand::Secondary);
  EXPECT_EQ(classifyRtkSignal(galileo, FREQ5, true, true),
            RtkSignalBand::L5);
  EXPECT_EQ(classifyRtkSignal(galileo, FREQ6, true, true),
            RtkSignalBand::Extra);
  EXPECT_EQ(classifyRtkSignal(galileo, FREQ8, true, true),
            RtkSignalBand::Wide);

  auto beidou = std::make_shared<gnss_comm::Obs>();
  beidou->sat = gnss_comm::sat_no(SYS_BDS, 12);
  EXPECT_EQ(classifyRtkSignal(beidou, FREQ2_BDS, true, true),
            RtkSignalBand::Secondary);
  EXPECT_EQ(classifyRtkSignal(beidou, FREQ3_BDS, true, true),
            RtkSignalBand::Extra);
  EXPECT_EQ(classifyRtkSignal(beidou, FREQ8, true, true),
            RtkSignalBand::Wide);
}

TEST(RtkSignalSelection, NormalizesLegacyGalileoGpsL2TagToE5b)
{
  auto galileo = std::make_shared<gnss_comm::Obs>();
  galileo->sat = gnss_comm::sat_no(SYS_GAL, 27);
  galileo->freqs = {FREQ1, FREQ2};

  EXPECT_DOUBLE_EQ(normalizeRtkFrequency(galileo->sat, FREQ2), FREQ7);
  EXPECT_EQ(classifyRtkSignal(galileo, FREQ2, true, true),
            RtkSignalBand::Secondary);
  EXPECT_EQ(classifyRtkSignal(galileo, FREQ2, true, false),
            RtkSignalBand::Unsupported);

  const std::vector<BaseCarrierObservation> base_epoch = {
      {galileo->sat, FREQ7, 123.0, 0.01, 0}};
  const BaseCarrierObservation *matched =
      findMatchingBaseCarrier(base_epoch, galileo->sat, FREQ2);
  ASSERT_NE(matched, nullptr);
  EXPECT_DOUBLE_EQ(matched->frequency_hz, FREQ7);
}

TEST(RtkSignalSelection, DoesNotNormalizeGpsL2)
{
  const uint32_t gps = gnss_comm::sat_no(SYS_GPS, 6);
  EXPECT_DOUBLE_EQ(normalizeRtkFrequency(gps, FREQ2), FREQ2);
}

TEST(RtkSignalSelection, MapsGlonassSecondaryAndExtraBands)
{
  auto glonass = std::make_shared<gnss_comm::Obs>();
  glonass->sat = gnss_comm::sat_no(SYS_GLO, 12);
  const double g2_channel_minus_two = FREQ2_GLO - 2.0 * DFRQ2_GLO;
  EXPECT_EQ(classifyRtkSignal(glonass, g2_channel_minus_two, true, true),
            RtkSignalBand::Secondary);
  EXPECT_EQ(classifyRtkSignal(glonass, FREQ3_GLO, true, true),
            RtkSignalBand::Extra);
}

TEST(RtkSignalSelection, ResolvesEveryGlonassG1FrequencyChannel)
{
  auto glonass = std::make_shared<gnss_comm::Obs>();
  glonass->sat = gnss_comm::sat_no(SYS_GLO, 12);
  for (int channel = -7; channel <= 6; ++channel)
  {
    const double expected = FREQ1_GLO + channel * DFRQ1_GLO;
    // Allow the source metadata to be slightly rounded; L1_freq() should
    // recover the channel centre and its observation-array index.
    glonass->freqs = {FREQ2_GLO, expected + 25000.0};
    int index = -1;
    EXPECT_DOUBLE_EQ(gnss_comm::L1_freq(glonass, &index), expected);
    EXPECT_EQ(index, 1);
    EXPECT_EQ(classifyRtkSignal(glonass, glonass->freqs[1], true, true),
              RtkSignalBand::Primary);
  }
}

TEST(RtkSignalSelection, RejectsNonChannelGlonassG1Metadata)
{
  auto glonass = std::make_shared<gnss_comm::Obs>();
  glonass->sat = gnss_comm::sat_no(SYS_GLO, 6);
  glonass->freqs = {FREQ1_GLO + 0.3 * DFRQ1_GLO};
  int index = 99;
  EXPECT_LT(gnss_comm::L1_freq(glonass, &index), 0.0);
  EXPECT_EQ(index, -1);
}

TEST(RtkArcValidation, InterpretsRinexLliAsBitMask)
{
  EXPECT_FALSE(rtkLossOfLock(0));
  EXPECT_TRUE(rtkLossOfLock(1));
  EXPECT_TRUE(rtkLossOfLock(2));
  EXPECT_TRUE(rtkLossOfLock(3));
  EXPECT_FALSE(rtkLossOfLock(4));  // Galileo BOC tracking, not a cycle slip.
  EXPECT_TRUE(rtkLossOfLock(5));

  EXPECT_FALSE(rtkCycleSlip(0));
  EXPECT_TRUE(rtkCycleSlip(1));
  EXPECT_FALSE(rtkCycleSlip(2));
  EXPECT_TRUE(rtkCycleSlip(3));
  EXPECT_FALSE(rtkHalfCycleInvalid(1));
  EXPECT_TRUE(rtkHalfCycleInvalid(2));
  EXPECT_TRUE(rtkHalfCycleInvalid(3));
  EXPECT_FALSE(rtkHalfCycleInvalid(4));
}

TEST(RtkArcValidation, UsesIndependentBaseEpochGapTolerance)
{
  EXPECT_TRUE(rtkArcIsContinuous(100.0, 101.0, 1.5));
  EXPECT_TRUE(rtkArcIsContinuous(100.0, 101.5, 1.5));
  EXPECT_FALSE(rtkArcIsContinuous(100.0, 101.5001, 1.5));
  EXPECT_FALSE(rtkArcIsContinuous(100.0, 99.0, 1.5));
  EXPECT_FALSE(rtkArcIsContinuous(0.0, 1.0, 1.5));
}

TEST(RtkArcValidation, GeometryFreeSlipHonorsEnableSwitch)
{
  EXPECT_FALSE(rtkGeometryFreeCycleSlip(false, 1.20, 1.40, 0.05));
  EXPECT_FALSE(rtkGeometryFreeCycleSlip(true, 1.20, 1.24, 0.05));
  EXPECT_TRUE(rtkGeometryFreeCycleSlip(true, 1.20, 1.251, 0.05));
  EXPECT_FALSE(rtkGeometryFreeCycleSlip(
      true, std::numeric_limits<double>::quiet_NaN(), 1.3, 0.05));
}

TEST(RtkNoiseModel, PropagatesAllFourCarrierUncertainties)
{
  const double satellite_variance = rtkSingleDifferenceVarianceMeters2(
      0.008, 0.190, 0.010, 0.190);
  const double reference_variance = rtkSingleDifferenceVarianceMeters2(
      0.020, 0.190, 0.010, 0.190);
  const double sigma = rtkDoubleDifferenceSigmaMeters(
      satellite_variance, reference_variance, 0.003);
  const double expected_variance =
      std::pow(0.008 * 0.190, 2) + std::pow(0.010 * 0.190, 2) +
      std::pow(0.020 * 0.190, 2) + std::pow(0.010 * 0.190, 2) +
      std::pow(0.003, 2);

  EXPECT_NEAR(sigma, std::sqrt(expected_variance), 1.0e-15);
  EXPECT_NE(sigma, 0.02);
}

TEST(RtkNoiseModel, RejectsInvalidUncertaintyInputs)
{
  EXPECT_FALSE(std::isfinite(rtkSingleDifferenceVarianceMeters2(
      -0.01, 0.19, 0.01, 0.19)));
  EXPECT_FALSE(std::isfinite(rtkDoubleDifferenceSigmaMeters(
      1.0e-5, -1.0e-5, 0.003)));
}

TEST(RtkResidualScreening, NormalizesAndRejectsCarrierOutliers)
{
  EXPECT_NEAR(rtkStandardizedResidual(0.024, 0.008), 3.0, 1.0e-12);
  EXPECT_TRUE(rtkResidualPassesGate(0.0239, 0.008, 3.0));
  EXPECT_TRUE(rtkResidualPassesGate(-0.024, 0.008, 3.0));
  EXPECT_FALSE(rtkResidualPassesGate(0.0241, 0.008, 3.0));
  EXPECT_FALSE(rtkResidualPassesGate(-0.0241, 0.008, 3.0));
  EXPECT_FALSE(rtkResidualPassesGate(0.01, 0.0, 3.0));
  EXPECT_FALSE(rtkResidualPassesGate(
      std::numeric_limits<double>::quiet_NaN(), 0.01, 3.0));
}

TEST(RtkSolutionStatus, ClassifiesNoAmbiguityFloatAndFixedStates)
{
  EXPECT_EQ(classifyRtkSolution(0, 0), RtkSolutionType::NoAmbiguity);
  EXPECT_STREQ(rtkSolutionTypeName(classifyRtkSolution(0, 0)),
               "NO_AMBIGUITY");
  EXPECT_EQ(classifyRtkSolution(4, 0), RtkSolutionType::Float);
  EXPECT_STREQ(rtkSolutionTypeName(classifyRtkSolution(4, 0)), "FLOAT");
  EXPECT_EQ(classifyRtkSolution(0, 3), RtkSolutionType::Fixed);
  EXPECT_EQ(classifyRtkSolution(5, 3), RtkSolutionType::Fixed);
  EXPECT_STREQ(rtkSolutionTypeName(classifyRtkSolution(5, 3)), "FIXED");
}

TEST(LambdaAmbiguityResolver, FindsKnownCorrelatedSolutionAndOrdersCandidates)
{
  Eigen::Vector3d floating(5.08, -3.04, 12.06);
  Eigen::Matrix3d covariance;
  covariance << 0.040, 0.025, 0.012,
                0.025, 0.035, 0.010,
                0.012, 0.010, 0.030;
  const LambdaResult result =
      LambdaAmbiguityResolver::solve(floating, covariance, 2);
  ASSERT_TRUE(result.valid) << result.error;
  ASSERT_EQ(result.candidates.cols(), 2);
  EXPECT_TRUE(result.candidates.col(0).isApprox(Eigen::Vector3d(5, -3, 12)));
  EXPECT_LT(result.squared_norms[0], result.squared_norms[1]);
}

TEST(LambdaAmbiguityResolver, MatchesDeterministicBruteForceCases)
{
  std::mt19937 generator(20260801);
  std::uniform_real_distribution<double> ambiguity_distribution(-3.0, 3.0);
  std::uniform_real_distribution<double> matrix_distribution(-0.8, 0.8);
  for (int trial = 0; trial < 40; ++trial)
  {
    const Eigen::Vector2d floating(ambiguity_distribution(generator),
                                   ambiguity_distribution(generator));
    Eigen::Matrix2d root;
    root << matrix_distribution(generator), matrix_distribution(generator),
            matrix_distribution(generator), matrix_distribution(generator);
    const Eigen::Matrix2d covariance =
        root * root.transpose() + 0.03 * Eigen::Matrix2d::Identity();
    double brute_force_norm = 0.0;
    const Eigen::Vector2d brute_force =
        bruteForceBest(floating, covariance, &brute_force_norm);
    const LambdaResult result =
        LambdaAmbiguityResolver::solve(floating, covariance, 2);
    ASSERT_TRUE(result.valid) << "trial " << trial << ": " << result.error;
    EXPECT_TRUE(result.candidates.col(0).isApprox(brute_force)) << "trial " << trial;
    EXPECT_NEAR(result.squared_norms[0], brute_force_norm, 1.0e-8)
        << "trial " << trial;
  }
}

TEST(LambdaAmbiguityResolver, IsInvariantToIntegerTranslation)
{
  const Eigen::Vector2d floating(0.23, -0.31);
  Eigen::Matrix2d covariance;
  covariance << 0.08, 0.06, 0.06, 0.09;
  const Eigen::Vector2d translation(17.0, -9.0);
  const LambdaResult original =
      LambdaAmbiguityResolver::solve(floating, covariance, 2);
  const LambdaResult shifted =
      LambdaAmbiguityResolver::solve(floating + translation, covariance, 2);
  ASSERT_TRUE(original.valid);
  ASSERT_TRUE(shifted.valid);
  EXPECT_TRUE(shifted.candidates.col(0).isApprox(
      original.candidates.col(0) + translation));
  EXPECT_NEAR(shifted.squared_norms[0], original.squared_norms[0], 1.0e-10);
}

TEST(LambdaAmbiguityResolver, RejectsInvalidCovariancesAndDimensions)
{
  const Eigen::Vector2d floating(1.1, 2.1);
  Eigen::Matrix2d indefinite;
  indefinite << 1.0, 2.0, 2.0, 1.0;
  EXPECT_FALSE(LambdaAmbiguityResolver::solve(floating, indefinite).valid);

  Eigen::Matrix2d asymmetric;
  asymmetric << 1.0, 0.2, 0.1, 1.0;
  EXPECT_FALSE(LambdaAmbiguityResolver::solve(floating, asymmetric).valid);
  EXPECT_FALSE(LambdaAmbiguityResolver::solve(
      floating, Eigen::Matrix3d::Identity()).valid);
}

TEST(RtkAmbiguityCascade, OrdersWideRawAndNarrowLanes)
{
  EXPECT_GT(rtkAmbiguityCascadeStage(RtkSignalBand::WideLane),
            rtkAmbiguityCascadeStage(RtkSignalBand::Primary));
  EXPECT_GT(rtkAmbiguityCascadeStage(RtkSignalBand::Primary),
            rtkAmbiguityCascadeStage(RtkSignalBand::NarrowLane));
}

TEST(RtkPartialAmbiguityResolution, SupportsAllGiciDeletionCriteria)
{
  Eigen::Vector3d floating(2.05, 4.49, -1.20);
  Eigen::Vector3d variances(0.01, 0.09, 0.04);
  Eigen::Vector3d elevations(0.8, 0.6, 0.2);
  EXPECT_EQ(rtkPartialFixWorstIndex(
                floating, variances, elevations,
                RtkPartialFixDeletionMethod::Elevation),
            2U);
  EXPECT_EQ(rtkPartialFixWorstIndex(
                floating, variances, elevations,
                RtkPartialFixDeletionMethod::Variance),
            1U);
  EXPECT_EQ(rtkPartialFixWorstIndex(
                floating, variances, elevations,
                RtkPartialFixDeletionMethod::Fractional),
            1U);
}

TEST(RtkPartialAmbiguityResolution, EnforcesMinimumFixPercentage)
{
  EXPECT_EQ(rtkMinimumPartialFixSize(10, 4, 0.7), 7U);
  EXPECT_EQ(rtkMinimumPartialFixSize(5, 4, 0.2), 4U);
  EXPECT_EQ(rtkMinimumPartialFixSize(5, 2, 1.5), 5U);
}

TEST(RtkFixConfirmation, RejectsPostFixCostIncrease)
{
  EXPECT_TRUE(rtkPostFixCostPasses(12.0, 12.01, 0.01));
  EXPECT_FALSE(rtkPostFixCostPasses(12.0, 12.011, 0.01));
  EXPECT_FALSE(rtkPostFixCostPasses(
      12.0, std::numeric_limits<double>::quiet_NaN(), 0.01));
}

TEST(RtkHistoricalMatch, PreservesIntegerAcrossReferenceChange)
{
  const std::vector<RtkFixedAmbiguityEdge> history = {
      {1, 3, 17},  // N3-N1
      {1, 2, 5}};  // N2-N1
  long long matched = 0;
  ASSERT_TRUE(rtkFindHistoricalMatch(2, 3, history, &matched));
  EXPECT_EQ(matched, 12);  // (N3-N1) - (N2-N1)
  ASSERT_TRUE(rtkFindHistoricalMatch(3, 2, history, &matched));
  EXPECT_EQ(matched, -12);
  EXPECT_FALSE(rtkFindHistoricalMatch(4, 3, history, &matched));
}

TEST(RtkSatelliteVisualization, ConvertsEcefLineOfSightToLocalFrame)
{
  const Eigen::Vector3d receiver(10.0, 20.0, 30.0);
  const Eigen::Vector3d satellite = receiver + Eigen::Vector3d(0.0, 20000.0, 0.0);
  const Eigen::Matrix3d ecef_from_local =
      Eigen::AngleAxisd(M_PI_2, Eigen::Vector3d::UnitZ()).toRotationMatrix();
  const Eigen::Vector3d local =
      rtkLocalLineOfSight(satellite, receiver, ecef_from_local);
  EXPECT_TRUE(local.isApprox(Eigen::Vector3d::UnitX(), 1.0e-12));
  EXPECT_NEAR(local.norm(), 1.0, 1.0e-12);
}

TEST(RtkSatelliteVisualization, EquivalentPointPreservesDirectionAndRadius)
{
  const Eigen::Vector3d receiver(4.0, -3.0, 2.0);
  const Eigen::Vector3d direction(2.0, 3.0, 6.0);
  const Eigen::Vector3d equivalent =
      rtkEquivalentSatellitePoint(receiver, direction, 50.0);
  EXPECT_NEAR((equivalent - receiver).norm(), 50.0, 1.0e-12);
  EXPECT_TRUE((equivalent - receiver).normalized().isApprox(
      direction.normalized(), 1.0e-12));
  EXPECT_TRUE(rtkEquivalentSatellitePoint(receiver, Eigen::Vector3d::Zero(), 50.0)
                  .isApprox(receiver));
}

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
