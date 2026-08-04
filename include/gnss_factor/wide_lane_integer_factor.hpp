#pragma once

#include <gtsam/base/Vector.h>
#include <gtsam/nonlinear/NoiseModelFactorN.h>

namespace ligo
{
// Enforces N_wl = A_primary/lambda_primary - A_secondary/lambda_secondary,
// while both ambiguity states remain in metres in the factor graph.
class WideLaneIntegerFactor
    : public gtsam::NoiseModelFactor2<gtsam::Vector1, gtsam::Vector1>
{
 public:
  WideLaneIntegerFactor(gtsam::Key primary_key, gtsam::Key secondary_key,
                        double primary_wavelength,
                        double secondary_wavelength,
                        long long integer_wide_lane,
                        const gtsam::SharedNoiseModel &model)
      : gtsam::NoiseModelFactor2<gtsam::Vector1, gtsam::Vector1>(
            model, primary_key, secondary_key),
        inverse_primary_wavelength_(1.0 / primary_wavelength),
        inverse_secondary_wavelength_(1.0 / secondary_wavelength),
        integer_wide_lane_(static_cast<double>(integer_wide_lane)) {}

  gtsam::Vector evaluateError(const gtsam::Vector1 &primary,
                              const gtsam::Vector1 &secondary,
                              gtsam::OptionalMatrixType H1,
                              gtsam::OptionalMatrixType H2) const override
  {
    if (H1) *H1 = (gtsam::Matrix(1, 1) << inverse_primary_wavelength_).finished();
    if (H2) *H2 = (gtsam::Matrix(1, 1) << -inverse_secondary_wavelength_).finished();
    return gtsam::Vector1(primary[0] * inverse_primary_wavelength_ -
                          secondary[0] * inverse_secondary_wavelength_ -
                          integer_wide_lane_);
  }

 private:
  double inverse_primary_wavelength_;
  double inverse_secondary_wavelength_;
  double integer_wide_lane_;
};
}  // namespace ligo
