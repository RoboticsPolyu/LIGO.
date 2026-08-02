#pragma once

#include <Eigen/Core>
#include <gtsam/base/Vector.h>
#include <gtsam/geometry/Rot3.h>
#include <gtsam/nonlinear/NoiseModelFactorN.h>

#include <so3_math.h>

namespace ligo
{
class DoubleDiffPseudorangeFactor
    : public gtsam::NoiseModelFactor3<gtsam::Vector3, gtsam::Rot3,
                                      gtsam::Vector6>
{
 public:
  DoubleDiffPseudorangeFactor(
      gtsam::Key anchor_key, gtsam::Key alignment_key, gtsam::Key state_key,
      const Eigen::Vector3d &antenna_offset_local,
      const Eigen::Vector3d &base_ecef,
      const Eigen::Vector3d &satellite_ecef,
      const Eigen::Vector3d &reference_ecef,
      double measured_double_difference,
      const gtsam::SharedNoiseModel &model)
      : gtsam::NoiseModelFactor3<gtsam::Vector3, gtsam::Rot3,
                                  gtsam::Vector6>(
            model, anchor_key, alignment_key, state_key),
        antenna_offset_local_(antenna_offset_local), base_ecef_(base_ecef),
        satellite_ecef_(satellite_ecef), reference_ecef_(reference_ecef),
        measured_(measured_double_difference) {}

  gtsam::Vector evaluateError(const gtsam::Vector3 &anchor,
                              const gtsam::Rot3 &alignment,
                              const gtsam::Vector6 &state,
                              gtsam::OptionalMatrixType H1,
                              gtsam::OptionalMatrixType H2,
                              gtsam::OptionalMatrixType H3) const override
  {
    const Eigen::Vector3d local = state.head<3>() + antenna_offset_local_;
    const Eigen::Vector3d rover = anchor + alignment.matrix() * local;
    const Eigen::Vector3d unit_sat = (satellite_ecef_ - rover).normalized();
    const Eigen::Vector3d unit_ref = (reference_ecef_ - rover).normalized();
    const Eigen::RowVector3d gradient = (unit_ref - unit_sat).transpose();
    Eigen::Matrix3d skew;
    skew << SKEW_SYM_MATRX(local);
    if (H1) *H1 = gradient;
    if (H2) *H2 = -gradient * alignment.matrix() * skew;
    if (H3)
    {
      *H3 = gtsam::Matrix::Zero(1, 6);
      (*H3).block<1, 3>(0, 0) = gradient * alignment.matrix();
    }
    return gtsam::Vector1(geometry(rover) - measured_);
  }

 private:
  double geometry(const Eigen::Vector3d &rover) const
  {
    return (satellite_ecef_ - rover).norm() -
           (satellite_ecef_ - base_ecef_).norm() -
           (reference_ecef_ - rover).norm() +
           (reference_ecef_ - base_ecef_).norm();
  }

  Eigen::Vector3d antenna_offset_local_, base_ecef_;
  Eigen::Vector3d satellite_ecef_, reference_ecef_;
  double measured_;
};

class DoubleDiffPseudorangeFactorNolidar
    : public gtsam::NoiseModelFactor2<gtsam::Rot3, gtsam::Vector12>
{
 public:
  DoubleDiffPseudorangeFactorNolidar(
      gtsam::Key rotation_key, gtsam::Key state_key,
      const Eigen::Vector3d &antenna_offset_imu,
      const Eigen::Vector3d &base_ecef,
      const Eigen::Vector3d &satellite_ecef,
      const Eigen::Vector3d &reference_ecef,
      double measured_double_difference,
      const gtsam::SharedNoiseModel &model)
      : gtsam::NoiseModelFactor2<gtsam::Rot3, gtsam::Vector12>(
            model, rotation_key, state_key),
        antenna_offset_imu_(antenna_offset_imu), base_ecef_(base_ecef),
        satellite_ecef_(satellite_ecef), reference_ecef_(reference_ecef),
        measured_(measured_double_difference) {}

  gtsam::Vector evaluateError(const gtsam::Rot3 &rotation,
                              const gtsam::Vector12 &state,
                              gtsam::OptionalMatrixType H1,
                              gtsam::OptionalMatrixType H2) const override
  {
    const Eigen::Vector3d rover =
        state.head<3>() + rotation * antenna_offset_imu_;
    const Eigen::Vector3d unit_sat = (satellite_ecef_ - rover).normalized();
    const Eigen::Vector3d unit_ref = (reference_ecef_ - rover).normalized();
    const Eigen::RowVector3d gradient = (unit_ref - unit_sat).transpose();
    Eigen::Matrix3d skew;
    skew << SKEW_SYM_MATRX(antenna_offset_imu_);
    if (H1) *H1 = -gradient * rotation.matrix() * skew;
    if (H2)
    {
      *H2 = gtsam::Matrix::Zero(1, 12);
      (*H2).block<1, 3>(0, 0) = gradient;
    }
    const double geometry = (satellite_ecef_ - rover).norm() -
                            (satellite_ecef_ - base_ecef_).norm() -
                            (reference_ecef_ - rover).norm() +
                            (reference_ecef_ - base_ecef_).norm();
    return gtsam::Vector1(geometry - measured_);
  }

 private:
  Eigen::Vector3d antenna_offset_imu_, base_ecef_;
  Eigen::Vector3d satellite_ecef_, reference_ecef_;
  double measured_;
};
}  // namespace ligo
