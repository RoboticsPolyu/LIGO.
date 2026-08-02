#pragma once

#include <gtsam/base/Vector.h>
#include <gtsam/nonlinear/NoiseModelFactorN.h>

#include <Eigen/Core>

#include <utility>
#include <vector>

namespace ligo
{
// Integer wide-lane relation N_WL = N_primary - N_secondary. The graph stores
// ambiguities in metres, hence the wavelength scaling in the residual.
class WideLaneAmbiguityFactor
    : public gtsam::NoiseModelFactor2<gtsam::Vector1, gtsam::Vector1>
{
 public:
  WideLaneAmbiguityFactor(gtsam::Key primary_key, gtsam::Key secondary_key,
                          double primary_wavelength,
                          double secondary_wavelength,
                          long long wide_lane_integer,
                          const gtsam::SharedNoiseModel &model)
      : gtsam::NoiseModelFactor2<gtsam::Vector1, gtsam::Vector1>(
            model, primary_key, secondary_key),
        primary_wavelength_(primary_wavelength),
        secondary_wavelength_(secondary_wavelength),
        wide_lane_integer_(wide_lane_integer) {}

  gtsam::Vector evaluateError(const gtsam::Vector1 &primary,
                              const gtsam::Vector1 &secondary,
                              gtsam::OptionalMatrixType H1,
                              gtsam::OptionalMatrixType H2) const override
  {
    if (H1) *H1 = (gtsam::Matrix(1, 1) << 1.0 / primary_wavelength_).finished();
    if (H2) *H2 = (gtsam::Matrix(1, 1) << -1.0 / secondary_wavelength_).finished();
    return gtsam::Vector1(primary[0] / primary_wavelength_ -
                          secondary[0] / secondary_wavelength_ -
                          static_cast<double>(wide_lane_integer_));
  }

 private:
  double primary_wavelength_;
  double secondary_wavelength_;
  long long wide_lane_integer_;
};

// T maps raw ambiguities in cycles to correlated wide-lane differences. It is
// deliberately exposed for covariance propagation and deterministic testing.
inline Eigen::MatrixXd wideLaneTransform(
    size_t raw_dimension,
    const std::vector<std::pair<size_t, size_t>> &primary_secondary_indices)
{
  Eigen::MatrixXd transform = Eigen::MatrixXd::Zero(
      primary_secondary_indices.size(), raw_dimension);
  for (size_t row = 0; row < primary_secondary_indices.size(); ++row)
  {
    transform(row, primary_secondary_indices[row].first) = 1.0;
    transform(row, primary_secondary_indices[row].second) = -1.0;
  }
  return transform;
}
}  // namespace ligo
