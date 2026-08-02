#include "LambdaAmbiguityResolver.h"

#include <Eigen/Cholesky>
#include <Eigen/LU>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace
{
// Deterministic nearest-integer convention used by both reduction and search.
double nearestInteger(double value) { return std::floor(value + 0.5); }
// Initial direction for zig-zag enumeration around the nearest integer.
double sign(double value) { return value <= 0.0 ? -1.0 : 1.0; }
}  // namespace

LambdaResult LambdaAmbiguityResolver::solve(
    const Eigen::VectorXd &float_ambiguities, const Eigen::MatrixXd &covariance,
    int candidate_count)
{
  LambdaResult result;
  const int dimension = float_ambiguities.size();

  // Reject malformed inputs before factorization. candidate_count >= 2 is an
  // intentional API requirement: the RTK acceptance test compares the best
  // and second-best candidates.
  if (dimension == 0 || covariance.rows() != dimension ||
      covariance.cols() != dimension || candidate_count < 2)
  {
    result.error = "invalid LAMBDA dimensions";
    return result;
  }
  if (!float_ambiguities.allFinite() || !covariance.allFinite() ||
      !covariance.isApprox(covariance.transpose(), 1e-9))
  {
    result.error = "non-finite or asymmetric ambiguity covariance";
    return result;
  }

  // Q = L^T D L. The reduction stage constructs an integer unimodular matrix
  // Z such that z_float = Z^T a_float has a less-correlated covariance.
  Eigen::MatrixXd lower, transform;
  Eigen::VectorXd diagonal;
  if (!ldFactor(covariance, lower, diagonal))
  {
    result.error = "ambiguity covariance is not positive definite";
    return result;
  }
  transform = Eigen::MatrixXd::Identity(dimension, dimension);
  reduce(lower, diagonal, transform);

  // Search in the decorrelated basis because its conditional variances form a
  // much smaller enumeration tree than the original correlated ambiguities.
  const Eigen::VectorXd decorrelated_float = transform.transpose() * float_ambiguities;
  Eigen::MatrixXd decorrelated_integers;
  if (!search(lower, diagonal, decorrelated_float, candidate_count,
              decorrelated_integers, result.squared_norms))
  {
    result.error = "LAMBDA search did not converge";
    return result;
  }

  // Map integer candidates back with a linear solve instead of explicitly
  // forming (Z^T)^-1. Z is unimodular, so the exact result is integer; final
  // rounding only removes floating-point error from the dense solve.
  result.candidates = transform.transpose().fullPivLu().solve(decorrelated_integers);
  for (int column = 0; column < result.candidates.cols(); ++column)
    for (int row = 0; row < result.candidates.rows(); ++row)
      result.candidates(row, column) = nearestInteger(result.candidates(row, column));
  result.valid = result.candidates.allFinite() && result.squared_norms.allFinite();
  if (!result.valid) result.error = "invalid transformed integer candidates";
  return result;
}

bool LambdaAmbiguityResolver::ldFactor(const Eigen::MatrixXd &covariance,
                                       Eigen::MatrixXd &lower,
                                       Eigen::VectorXd &diagonal)
{
  const int dimension = covariance.rows();
  // Work backward through Q, matching the conventional LAMBDA factorization
  // Q = L^T diag(D) L. Only the lower triangle of workspace is required.
  Eigen::MatrixXd workspace = covariance;
  lower = Eigen::MatrixXd::Zero(dimension, dimension);
  diagonal.resize(dimension);
  for (int row = dimension - 1; row >= 0; --row)
  {
    // A non-positive conditional variance means Q is not positive definite.
    diagonal[row] = workspace(row, row);
    if (diagonal[row] <= 0.0 || !std::isfinite(diagonal[row])) return false;
    const double root = std::sqrt(diagonal[row]);
    // Normalize the active row, then apply its rank-one Schur complement to
    // the unprocessed leading block.
    for (int column = 0; column <= row; ++column)
      lower(row, column) = workspace(row, column) / root;
    for (int column = 0; column < row; ++column)
      for (int inner = 0; inner <= column; ++inner)
        workspace(column, inner) -= lower(row, inner) * lower(row, column);
    for (int column = 0; column <= row; ++column)
      lower(row, column) /= lower(row, row);
  }
  return true;
}

void LambdaAmbiguityResolver::integerGauss(Eigen::MatrixXd &lower,
                                            Eigen::MatrixXd &transform,
                                            int row, int column)
{
  // Subtract the nearest integer multiple of column 'row' from 'column'. This
  // reduces an off-diagonal L entry without changing the integer lattice.
  const double integer = nearestInteger(lower(row, column));
  if (integer == 0.0) return;
  for (int index = row; index < lower.rows(); ++index)
    lower(index, column) -= integer * lower(index, row);
  for (int index = 0; index < transform.rows(); ++index)
    // Accumulate the identical column operation in the unimodular transform Z.
    transform(index, column) -= integer * transform(index, row);
}

void LambdaAmbiguityResolver::permute(Eigen::MatrixXd &lower,
                                      Eigen::VectorXd &diagonal, int column,
                                      double delta, Eigen::MatrixXd &transform)
{
  const int dimension = lower.rows();
  // Reorder adjacent ambiguities when doing so improves the sequence of
  // conditional variances. eta and lambda update the affected 2-by-2 block
  // without reconstructing the full covariance.
  const double eta = diagonal[column] / delta;
  const double lambda = diagonal[column + 1] * lower(column + 1, column) / delta;
  diagonal[column] = eta * diagonal[column + 1];
  diagonal[column + 1] = delta;
  for (int index = 0; index < column; ++index)
  {
    const double first = lower(column, index);
    const double second = lower(column + 1, index);
    lower(column, index) = -lower(column + 1, column) * first + second;
    lower(column + 1, index) = eta * first + lambda * second;
  }
  lower(column + 1, column) = lambda;
  for (int index = column + 2; index < dimension; ++index)
    std::swap(lower(index, column), lower(index, column + 1));
  // The same permutation must be applied to Z to preserve z = Z^T a.
  transform.col(column).swap(transform.col(column + 1));
}

void LambdaAmbiguityResolver::reduce(Eigen::MatrixXd &lower,
                                     Eigen::VectorXd &diagonal,
                                     Eigen::MatrixXd &transform)
{
  const int dimension = lower.rows();
  // Work from right to left. A successful permutation can make later columns
  // reducible again, so restart at the rightmost adjacent pair.
  int column = dimension - 2;
  int reduced_through = column;
  while (column >= 0)
  {
    if (column <= reduced_through)
      for (int row = column + 1; row < dimension; ++row)
        integerGauss(lower, transform, row, column);
    const double delta = diagonal[column] +
                         lower(column + 1, column) * lower(column + 1, column) *
                             diagonal[column + 1];
    if (delta + 1e-6 < diagonal[column + 1])
    {
      // The small tolerance avoids cycling on numerically equal variances.
      permute(lower, diagonal, column, delta, transform);
      reduced_through = column;
      column = dimension - 2;
    }
    else
    {
      --column;
    }
  }
}

bool LambdaAmbiguityResolver::search(
    const Eigen::MatrixXd &lower, const Eigen::VectorXd &diagonal,
    const Eigen::VectorXd &decorrelated_float, int candidate_count,
    Eigen::MatrixXd &integers, Eigen::VectorXd &squared_norms)
{
  // This is the MLAMBDA-style depth-first integer search. The hard iteration
  // bound prevents pathological covariance/input combinations from hanging an
  // online estimator; failure simply leaves the GNSS pipeline in float mode.
  constexpr int maximum_iterations = 100000;
  const int dimension = lower.rows();

  // sums stores recursively accumulated L terms used to form conditional
  // ambiguity centers. distance stores the partial squared norm below a tree
  // level. candidate is the current integer vector, and step generates the
  // nearest, next-nearest, ... zig-zag order at each level.
  Eigen::MatrixXd sums = Eigen::MatrixXd::Zero(dimension, dimension);
  Eigen::VectorXd distance = Eigen::VectorXd::Zero(dimension);
  Eigen::VectorXd conditional = Eigen::VectorXd::Zero(dimension);
  Eigen::VectorXd candidate = Eigen::VectorXd::Zero(dimension);
  Eigen::VectorXd step = Eigen::VectorXd::Zero(dimension);
  integers = Eigen::MatrixXd::Zero(dimension, candidate_count);
  squared_norms = Eigen::VectorXd::Constant(
      candidate_count, std::numeric_limits<double>::infinity());

  // Start at the last variable because L is lower triangular, then descend
  // toward level zero as long as the partial norm remains inside the current
  // search ellipsoid.
  int level = dimension - 1;
  int found = 0;
  int worst = 0;
  double maximum_distance = std::numeric_limits<double>::infinity();
  conditional[level] = decorrelated_float[level];
  candidate[level] = nearestInteger(conditional[level]);
  double offset = conditional[level] - candidate[level];
  step[level] = sign(offset);

  for (int iteration = 0; iteration < maximum_iterations; ++iteration)
  {
    // Incremental objective contribution for the integer selected at 'level'.
    const double new_distance =
        distance[level] + offset * offset / diagonal[level];
    if (new_distance < maximum_distance)
    {
      if (level != 0)
      {
        // Descend one level and calculate its conditional float ambiguity from
        // the already fixed integers at deeper levels.
        --level;
        distance[level] = new_distance;
        sums(level, level) = 0.0;
        for (int index = 0; index <= level; ++index)
          sums(level, index) = sums(level + 1, index) +
                               (candidate[level + 1] - conditional[level + 1]) *
                                   lower(level + 1, index);
        conditional[level] = decorrelated_float[level] + sums(level, level);
        candidate[level] = nearestInteger(conditional[level]);
        offset = conditional[level] - candidate[level];
        step[level] = sign(offset);
      }
      else
      {
        // A leaf is a complete integer vector. Keep only candidate_count best
        // leaves and use the current worst retained norm as the pruning radius.
        if (found < candidate_count)
        {
          integers.col(found) = candidate;
          squared_norms[found] = new_distance;
          if (found == 0 || new_distance > squared_norms[worst]) worst = found;
          ++found;
        }
        else if (new_distance < squared_norms[worst])
        {
          integers.col(worst) = candidate;
          squared_norms[worst] = new_distance;
          worst = 0;
          for (int index = 1; index < candidate_count; ++index)
            if (squared_norms[index] > squared_norms[worst]) worst = index;
        }
        if (found == candidate_count) maximum_distance = squared_norms[worst];
        // Continue enumerating level-zero integers in alternating distance
        // order around its conditional center.
        candidate[0] += step[0];
        offset = conditional[0] - candidate[0];
        step[0] = -step[0] - sign(step[0]);
      }
    }
    else
    {
      // This branch cannot beat the retained candidates. Backtrack and try the
      // next integer at the parent level; termination occurs after exhausting
      // the root level's bounded ellipsoid.
      if (level == dimension - 1) break;
      ++level;
      candidate[level] += step[level];
      offset = conditional[level] - candidate[level];
      step[level] = -step[level] - sign(step[level]);
    }
  }
  if (found < candidate_count) return false;

  // Search storage order depends on discovery/replacement order. Normalize the
  // public result so column zero is always the maximum-likelihood candidate.
  std::vector<int> order(candidate_count);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int left, int right) {
    return squared_norms[left] < squared_norms[right];
  });
  const Eigen::MatrixXd unsorted_integers = integers;
  const Eigen::VectorXd unsorted_norms = squared_norms;
  for (int index = 0; index < candidate_count; ++index)
  {
    integers.col(index) = unsorted_integers.col(order[index]);
    squared_norms[index] = unsorted_norms[order[index]];
  }
  return true;
}
