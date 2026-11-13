/** @file local_curvature.hpp
 * Validated inverse joint-likelihood curvature, without ridge regularization.
 */
#pragma once
#include "delta_method.hpp"
namespace swordfish {
/// Full optimizer-coordinate covariance plus two-step numerical diagnostics.
struct LocalCurvature {
  Eigen::MatrixXd covariance, coarse_covariance;
  double asymmetry, step_difference, minimum_eigenvalue, maximum_eigenvalue,
      maximum_se_step_change;
};
/// Validate two finite-difference Hessians and invert positive local curvature.
/// Both steps must be symmetric/stable to 1e-5 relative norm, positive definite,
/// and give marginal SEs stable to 0.1%. Singular/indefinite curvature is rejected;
/// no eigenvalue clipping, ridge, diagonal-only approximation or saved covariance.
inline LocalCurvature validated_local_curvature(const Eigen::MatrixXd &coarse,
                                                const Eigen::MatrixXd &fine) {
  if(coarse.rows()==0||coarse.rows()!=coarse.cols()||fine.rows()!=coarse.rows()||fine.cols()!=coarse.cols()||!coarse.allFinite()||!fine.allFinite()||!(fine.norm()>0)||!(coarse.norm()>0))
    throw std::invalid_argument("Invalid local Hessian dimensions/values");
  const double asym=std::max((fine-fine.transpose()).norm()/fine.norm(),(coarse-coarse.transpose()).norm()/coarse.norm());
  const double difference=(coarse-fine).norm()/fine.norm();
  if(asym>1e-5||difference>1e-5)throw std::runtime_error("Local Hessian symmetry/step validation failed");
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> a((coarse+coarse.transpose())*.5),b((fine+fine.transpose())*.5);
  if(a.info()!=Eigen::Success||b.info()!=Eigen::Success||a.eigenvalues().minCoeff()<=0||b.eigenvalues().minCoeff()<=0)
    throw std::runtime_error("Local Hessian is not positive definite");
  const Eigen::MatrixXd vc=a.eigenvectors()*a.eigenvalues().cwiseInverse().asDiagonal()*a.eigenvectors().transpose();
  const Eigen::MatrixXd vf=b.eigenvectors()*b.eigenvalues().cwiseInverse().asDiagonal()*b.eigenvectors().transpose();
  validate_parameter_covariance(vc);validate_parameter_covariance(vf);
  double change=0.;for(int i=0;i<vf.rows();++i)change=std::max(change,std::abs(std::sqrt(vc(i,i)/vf(i,i))-1.));
  if(change>1e-3)throw std::runtime_error("Local parameter SEs sensitive to Hessian step");
  return {vf,vc,asym,difference,b.eigenvalues().minCoeff(),b.eigenvalues().maxCoeff(),change};
}
} // namespace swordfish
