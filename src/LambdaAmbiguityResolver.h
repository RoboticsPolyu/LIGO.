#pragma once

#include <Eigen/Core>

#include <string>

/** Result of the integer least-squares ambiguity search.
 *
 * candidates has one integer ambiguity vector per column, expressed in the
 * original ambiguity basis. squared_norms contains the corresponding values
 * of (a-z)^T Q^-1 (a-z), sorted from best to worst. An invalid result carries
 * a human-readable error and must not be used to constrain the graph.
 */
struct LambdaResult
{
  bool valid = false;
  Eigen::MatrixXd candidates;
  Eigen::VectorXd squared_norms;
  std::string error;
};

class LambdaAmbiguityResolver
{
 public:
  /**
   * Solve the LAMBDA integer least-squares problem
   *   min_z (a-z)^T Q^-1 (a-z), z in Z^n.
   *
   * The caller supplies float ambiguities and their full covariance in the
   * same units (the GNSS pipeline converts metres to cycles first). At least
   * two candidates are required because ambiguity validation uses their ratio.
   * Candidate columns are returned in increasing squared-norm order.
   */
  static LambdaResult solve(const Eigen::VectorXd &float_ambiguities,
                            const Eigen::MatrixXd &covariance,
                            int candidate_count = 2);

 private:
  /** Factor Q as L^T diag(D) L, where L is unit lower triangular. */
  static bool ldFactor(const Eigen::MatrixXd &covariance,
                       Eigen::MatrixXd &lower,
                       Eigen::VectorXd &diagonal);
  /** Apply an integer Gauss column transform and accumulate its unimodular Z. */
  static void integerGauss(Eigen::MatrixXd &lower, Eigen::MatrixXd &transform,
                           int row, int column);
  /** Swap adjacent ambiguities and update L, D, and the accumulated transform. */
  static void permute(Eigen::MatrixXd &lower, Eigen::VectorXd &diagonal,
                      int column, double delta, Eigen::MatrixXd &transform);
  /** Decorrelate the ambiguities using integer Gauss transforms and permutations. */
  static void reduce(Eigen::MatrixXd &lower, Eigen::VectorXd &diagonal,
                     Eigen::MatrixXd &transform);
  /** Enumerate the best integers in the decorrelated basis using a bounded
   * depth-first zig-zag search (nearest integer, then alternating neighbors).
   */
  static bool search(const Eigen::MatrixXd &lower,
                     const Eigen::VectorXd &diagonal,
                     const Eigen::VectorXd &decorrelated_float,
                     int candidate_count, Eigen::MatrixXd &integers,
                     Eigen::VectorXd &squared_norms);
};
