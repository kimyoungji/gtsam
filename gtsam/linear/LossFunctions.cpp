/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file NoiseModel.cpp
 * @date Jan 13, 2010
 * @author Richard Roberts
 * @author Frank Dellaert
 */

#include <gtsam/linear/LossFunctions.h>

#include <iostream>

using namespace std;

namespace gtsam {
namespace noiseModel {

/* ************************************************************************* */
// M-Estimator
/* ************************************************************************* */

namespace mEstimator {

Vector Base::weight(const Vector& error) const {
  const size_t n = error.rows();
  Vector w(n);
  for (size_t i = 0; i < n; ++i)
    w(i) = weight(error(i));
  return w;
}

// The following three functions re-weight block matrices and a vector
// according to their weight implementation

void Base::reweight(Vector& error) const {
  if (reweight_ == Scalar) {
    const double w = sqrtWeight(error.norm());// customweight(error);//
    error *= w;
  } else if(reweight_== Block){
    error.array() *= weight(error).cwiseSqrt().array();
  }
  else {
    //Vector W = error.cwiseSqrt();
    error.array() = pairweight(error).cwiseProduct(error).array();
    //error = sqrtWeight(pairweight(W)).cwiseProduct(error);//.cwiseSqrt().array();
  }
}

// Reweight n block matrices with one error vector
void Base::reweight(vector<Matrix> &A, Vector &error) const {
  if ( reweight_ == Scalar ) {
    const double w = sqrtWeight(error.norm());
    for(Matrix& Aj: A) {
      Aj *= w;
    }
    error *= w;
  }
  else if(reweight_== Block){
    const Vector W = sqrtWeight(error);
    for(Matrix& Aj: A) {
      vector_scale_inplace(W,Aj);
    }
    error.array() = W.cwiseProduct(error).array();
  }
  else{
    //Vector W_pre = error.cwiseSqrt();
    const Vector W = pairweight(error);
    for(Matrix& Aj: A) {
      vector_scale_inplace(W,Aj);
    }
    error.array() = W.cwiseProduct(error).array();
  }
}

// Reweight one block matrix with one error vector
void Base::reweight(Matrix &A, Vector &error) const {
  if ( reweight_ == Scalar ) {
    const double w = sqrtWeight(error.norm());
    A *= w;
    error *= w;
  }
  else if(reweight_== Block){
    Vector W = sqrtWeight(error);
    vector_scale_inplace(W,A);
    error.array() = W.cwiseProduct(error).array();
  }
  else{
    //Vector W_pre = error.cwiseSqrt();
    const Vector W = pairweight(error);
    vector_scale_inplace(W,A);
    error.array() = W.cwiseProduct(error).array();
  }
}

// Reweight two block matrix with one error vector
void Base::reweight(Matrix &A1, Matrix &A2, Vector &error) const {
  if ( reweight_ == Scalar ) {
    const double w = sqrtWeight(error.norm());
    A1 *= w;
    A2 *= w;
    error *= w;
  }
  else if(reweight_== Block){
    const Vector W = sqrtWeight(error);
    vector_scale_inplace(W,A1);
    vector_scale_inplace(W,A2);
    error.array() = W.cwiseProduct(error).array();
  }
  else{
    //Vector W_pre = error.cwiseSqrt();
    const Vector W = pairweight(error);
    vector_scale_inplace(W,A1);
    vector_scale_inplace(W,A2);
    error.array() = W.cwiseProduct(error).array();
  }
}

// Reweight three block matrix with one error vector
void Base::reweight(Matrix &A1, Matrix &A2, Matrix &A3, Vector &error) const {
  if ( reweight_ == Scalar) {
    const double w = sqrtWeight(error.norm());
    A1 *= w;
    A2 *= w;
    A3 *= w;
    error *= w;
  }
  else if(reweight_== Block ){
    const Vector W = sqrtWeight(error);
    vector_scale_inplace(W,A1);
    vector_scale_inplace(W,A2);
    vector_scale_inplace(W,A3);
    error.array() = W.cwiseProduct(error).array();
  }
  else{
    //Vector W_pre = error.cwiseSqrt();
    const Vector W = pairweight(error);
    vector_scale_inplace(W,A1);
    vector_scale_inplace(W,A2);
    vector_scale_inplace(W,A3);
    error.array() = W.cwiseProduct(error).array();
  }
}

/* ************************************************************************* */
// Null model
/* ************************************************************************* */

void Null::print(const std::string &s="") const
{ cout << s << "null ()" << endl; }

Null::shared_ptr Null::Create()
{ return shared_ptr(new Null()); }

/* ************************************************************************* */
// Fair
/* ************************************************************************* */

Fair::Fair(double c, const ReweightScheme reweight) : Base(reweight), c_(c) {
  if (c_ <= 0) {
    throw runtime_error("mEstimator Fair takes only positive double in constructor.");
  }
}

double Fair::weight(double error) const {
  return 1.0 / (1.0 + std::abs(error) / c_);
}

double Fair::residual(const Vector error) const {
  const double absError = std::abs(error.norm());
  const double normalizedError = absError / c_;
  const double c_2 = c_ * c_;
  return c_2 * (normalizedError - std::log1p(normalizedError));
}

void Fair::print(const std::string &s="") const
{ cout << s << "fair (" << c_ << ")" << endl; }

bool Fair::equals(const Base &expected, double tol) const {
  const Fair* p = dynamic_cast<const Fair*> (&expected);
  if (p == NULL) return false;
  return std::abs(c_ - p->c_ ) < tol;
}

Fair::shared_ptr Fair::Create(double c, const ReweightScheme reweight)
{ return shared_ptr(new Fair(c, reweight)); }

/* ************************************************************************* */
// Huber
/* ************************************************************************* */

Huber::Huber(double k, const ReweightScheme reweight) : Base(reweight), k_(k) {
  if (k_ <= 0) {
    throw runtime_error("mEstimator Huber takes only positive double in constructor.");
  }
}

double Huber::weight(double error) const {
  const double absError = std::abs(error);
  return (absError <= k_) ? (1.0) : (k_ / absError);
}

double Huber::residual(const Vector error) const {
  const double absError = std::abs(error.norm());
  if (absError <= k_) {  // |x| <= k
    return error.norm()*error.norm() / 2;
  } else { // |x| > k
    return k_ * (absError - (k_/2));
  }
}

void Huber::print(const std::string &s="") const {
  cout << s << "huber (" << k_ << ")" << endl;
}

bool Huber::equals(const Base &expected, double tol) const {
  const Huber* p = dynamic_cast<const Huber*>(&expected);
  if (p == NULL) return false;
  return std::abs(k_ - p->k_) < tol;
}

Huber::shared_ptr Huber::Create(double c, const ReweightScheme reweight) {
  return shared_ptr(new Huber(c, reweight));
}

/* ************************************************************************* */
// Cauchy
/* ************************************************************************* */

Cauchy::Cauchy(double k, const ReweightScheme reweight) : Base(reweight), k_(k), ksquared_(k * k) {
  if (k <= 0) {
    throw runtime_error("mEstimator Cauchy takes only positive double in constructor.");
  }
}

double Cauchy::weight(double error) const {
  return ksquared_ / (ksquared_ + error*error);
}

double Cauchy::residual(const Vector error) const {
  const double val = std::log1p(error.norm() * error.norm() / ksquared_);
  return ksquared_ * val * 0.5;
}

void Cauchy::print(const std::string &s="") const {
  cout << s << "cauchy (" << k_ << ")" << endl;
}

bool Cauchy::equals(const Base &expected, double tol) const {
  const Cauchy* p = dynamic_cast<const Cauchy*>(&expected);
  if (p == NULL) return false;
  return std::abs(ksquared_ - p->ksquared_) < tol;
}

Cauchy::shared_ptr Cauchy::Create(double c, const ReweightScheme reweight) {
  return shared_ptr(new Cauchy(c, reweight));
}

/* ************************************************************************* */
// Tukey
/* ************************************************************************* */

Tukey::Tukey(double c, const ReweightScheme reweight) : Base(reweight), c_(c), csquared_(c * c) {
  if (c <= 0) {
    throw runtime_error("mEstimator Tukey takes only positive double in constructor.");
  }
}

double Tukey::weight(double error) const {
  if (std::abs(error) <= c_) {
    const double one_minus_xc2 = 1.0 - error*error/csquared_;
    return one_minus_xc2 * one_minus_xc2;
  }
  return 0.0;
}

double Tukey::residual(const Vector error) const {
  double absError = std::abs(error.norm());
  if (absError <= c_) {
    const double one_minus_xc2 = 1.0 - error.norm()*error.norm()/csquared_;
    const double t = one_minus_xc2*one_minus_xc2*one_minus_xc2;
    return csquared_ * (1 - t) / 6.0;
  } else {
    return csquared_ / 6.0;
  }
}

void Tukey::print(const std::string &s="") const {
  std::cout << s << ": Tukey (" << c_ << ")" << std::endl;
}

bool Tukey::equals(const Base &expected, double tol) const {
  const Tukey* p = dynamic_cast<const Tukey*>(&expected);
  if (p == NULL) return false;
  return std::abs(c_ - p->c_) < tol;
}

Tukey::shared_ptr Tukey::Create(double c, const ReweightScheme reweight) {
  return shared_ptr(new Tukey(c, reweight));
}

/* ************************************************************************* */
// Welsch
/* ************************************************************************* */

Welsch::Welsch(double c, const ReweightScheme reweight) : Base(reweight), c_(c), csquared_(c * c) {}

double Welsch::weight(double error) const {
  const double xc2 = (error*error)/csquared_;
  return std::exp(-xc2);
}

double Welsch::residual(const Vector error) const {
  const double xc2 = (error.norm()*error.norm())/csquared_;
  return csquared_ * 0.5 * -std::expm1(-xc2);
}

void Welsch::print(const std::string &s="") const {
  std::cout << s << ": Welsch (" << c_ << ")" << std::endl;
}

bool Welsch::equals(const Base &expected, double tol) const {
  const Welsch* p = dynamic_cast<const Welsch*>(&expected);
  if (p == NULL) return false;
  return std::abs(c_ - p->c_) < tol;
}

Welsch::shared_ptr Welsch::Create(double c, const ReweightScheme reweight) {
  return shared_ptr(new Welsch(c, reweight));
}

/* ************************************************************************* */
// GemanMcClure
/* ************************************************************************* */
GemanMcClure::GemanMcClure(double c, const ReweightScheme reweight)
  : Base(reweight), c_(c) {
}

double GemanMcClure::weight(double error) const {
  const double c2 = c_*c_;
  const double c4 = c2*c2;
  const double c2error = c2 + error*error;
  return c4/(c2error*c2error);
}

double GemanMcClure::residual(const Vector error) const {
  const double c2 = c_*c_;
  const double error2 = error.norm()*error.norm();
  return 0.5 * (c2 * error2) / (c2 + error2);
}

void GemanMcClure::print(const std::string &s="") const {
  std::cout << s << ": Geman-McClure (" << c_ << ")" << std::endl;
}

bool GemanMcClure::equals(const Base &expected, double tol) const {
  const GemanMcClure* p = dynamic_cast<const GemanMcClure*>(&expected);
  if (p == NULL) return false;
  return std::abs(c_ - p->c_) < tol;
}

GemanMcClure::shared_ptr GemanMcClure::Create(double c, const ReweightScheme reweight) {
  return shared_ptr(new GemanMcClure(c, reweight));
}

/* ************************************************************************* */
// DCS
/* ************************************************************************* */
DCS::DCS(double c, const ReweightScheme reweight)
  : Base(reweight), c_(c) {
}

double DCS::weight(double error) const {
  const double e2 = error*error;
  if (e2 > c_)
  {
    const double w = 2.0*c_/(c_ + e2);
    return w*w;
  }

  return 1.0;
}

double DCS::residual(const Vector error) const {
  // This is the simplified version of Eq 9 from (Agarwal13icra)
  // after you simplify and cancel terms.
  const double e2 = error.norm()*error.norm();
  const double e4 = e2*e2;
  const double c2 = c_*c_;

  return (c2*e2 + c_*e4) / ((e2 + c_)*(e2 + c_));
}

void DCS::print(const std::string &s="") const {
  std::cout << s << ": DCS (" << c_ << ")" << std::endl;
}

bool DCS::equals(const Base &expected, double tol) const {
  const DCS* p = dynamic_cast<const DCS*>(&expected);
  if (p == NULL) return false;
  return std::abs(c_ - p->c_) < tol;
}

DCS::shared_ptr DCS::Create(double c, const ReweightScheme reweight) {
  return shared_ptr(new DCS(c, reweight));
}

/* ************************************************************************* */
// PairDCS
/* ************************************************************************* */
PairDCS::PairDCS(double c, const ReweightScheme reweight)
  : Base(reweight), c_(c) {
}

double PairDCS::weight(double error) const {
  const double e2 = error*error;
  //return e2;
  if (e2 > c_)
  {
    const double w = 2.0*c_/(c_ + e2);
    return w*w;
  }

  return 1.0;
}

//double PairDCS::weight(Vector& error) const {
//    int dim;
//  if (error.size()==4)dim = 2;
//  else dim = 3;
//  const double e2_pair1 = error.head(dim).transpose()*error.head(dim);//error.head(dim).norm();//
//  const double e2_pair2 = error.tail(dim).transpose()*error.tail(dim);//error.tail(dim).norm();//
//  const double h_mean = e2_pair1*e2_pair2/(e2_pair1 + e2_pair2);
//  return c_/(c_+h_mean);
//}
double PairDCS::customweight(Vector& error) const {
  int dim;
  if (error.size()==4)dim = 2;
  else dim = 3;
  const double e2_pair1 = error.head(dim).transpose()*error.head(dim);//error.head(dim).norm();//
  const double e2_pair2 = error.tail(dim).transpose()*error.tail(dim);//error.tail(dim).norm();//
  const double h_mean = e2_pair1*e2_pair2/(e2_pair1 + e2_pair2);
  double s = c_/(c_+h_mean);
//  if (h_mean>c_) s = 2.0*c_/(c_+h_mean);
//  else s = 1.0;
  return s;
}
Vector PairDCS::pairweight(Vector& error) const {

  int dim;
  if (error.size()==4)dim = 2;
  else dim = 3;
  const double e2_pair1 = error.head(dim).transpose()*error.head(dim);//error.head(dim).norm();//
  const double e2_pair2 = error.tail(dim).transpose()*error.tail(dim);//error.tail(dim).norm();//
  const double h_mean = e2_pair1*e2_pair2/(e2_pair1 + e2_pair2);
  double s = c_/(c_+h_mean);
  if (h_mean>c_) s = 2.0*c_/(c_+h_mean);
  else s = 1.0;
  double c1 = 2.0*e2_pair2/(e2_pair1 + e2_pair2);
  if (c1>1.0) c1 = 1.0;
  double c2 = 2.0*e2_pair1/(e2_pair1 + e2_pair2);
  if (c2>1.0) c2 = 1.0;

  if (dim ==2){
  return (Vector(4) << c1*s, c1*s, c2*s, c2*s).finished();//s, s, s, s).finished();//c1, c1, c2, c2).finished();//
  }
  else{
  return (Vector(6) << c1*s, c1*s, c1*s, c2*s, c2*s, c2*s).finished();//s, s, s, s, s, s).finished();//c1, c1, c1, c2, c2, c2).finished();//
  }
}

double PairDCS::residual(const Vector error) const {

  int dim;
  if (error.size()==4)dim = 2;
  else dim = 3;
  const double e2_pair1 = error.head(dim).transpose()*error.head(dim);//error.head(dim).norm();//
  const double e2_pair2 = error.tail(dim).transpose()*error.tail(dim);//error.tail(dim).norm();//
  const double h_mean = e2_pair1*e2_pair2/(e2_pair1 + e2_pair2);

  return h_mean*c_/(h_mean+c_);//error.transpose()*error;//
}

void PairDCS::print(const std::string &s="") const {
  std::cout << s << ": PairDCS (" << c_ << ")" << std::endl;
}

bool PairDCS::equals(const Base &expected, double tol) const {
  const PairDCS* p = dynamic_cast<const PairDCS*>(&expected);
  if (p == NULL) return false;
  return std::abs(c_ - p->c_) < tol;
}

PairDCS::shared_ptr PairDCS::Create(double c, const ReweightScheme reweight) {
  return shared_ptr(new PairDCS(c, reweight));
}

/* ************************************************************************* */
// MM for multi-modal
/* ************************************************************************* */
MM_multi::MM_multi( double c, const ReweightScheme reweight)
  : Base(reweight), c_(c){
}

double MM_multi::weight(double error) const {
  const double e2 = error*error;
  if (e2 > c_)
  {
    const double w = 2.0*c_/(c_ + e2);
    return w*w;
  }

  return 1.0;
}

Vector MM_multi::pairweight(Vector& error) const {

  int dim;
  if (error.size()==4)dim = 2;
  else dim = 3;
  Vector err1 = error.head(dim);
  Vector err2 = error.tail(dim);
  const double e2_pair1 = err1.transpose()*err1;//error.head(dim).norm();//
  const double e2_pair2 = err2.transpose()*err2;//error.tail(dim).norm();//

  double e2 = e2_pair2;
  Vector c1 = err2.cwiseProduct(err1.cwiseInverse());//Vector::Zero(dim);// 
  Vector c2 = err2.cwiseProduct(err2.cwiseInverse());//Vector::Ones(dim);// 
  if(e2_pair2>e2_pair1){
    e2 = e2_pair1;
//    err1 = c1;
//    c1 = c2;
//    c2 = err1;
    c1 = err1.cwiseProduct(err1.cwiseInverse());//Vector::Zero(dim);// 
    c2 = err1.cwiseProduct(err2.cwiseInverse());//Vector::Ones(dim);// 
  }
  //const double e2 = min(e2_pair1, e2_pair2);
  double s = 1.0;
  if(e2>c_) s = 2.0*c_/(c_ + e2);
  if (dim ==2)return (Vector(4) << s*c1(0), s*c1(1), s*c2(0), s*c2(1)).finished();//s*c1,s*c1,s*c2,s*c2).finished();//
  else return (Vector(6) << s*c1(0), s*c1(1), s*c1(2), s*c2(0), s*c2(1), s*c2(2)).finished();//s*c1,s*c1,s*c1,s*c2,s*c2,s*c2).finished();//

//  if (e2 > log(2*3.141592/(b_*b_))/(a_*a_-b_*b_))
//  {
//    const double w = b_*b_/(a_*a_) + log(2*3.141592/(b_*b_))/e2;
//    if (dim ==2)return (Vector(4) << w, w, w, w).finished();
//    else return (Vector(6) << w, w, w, w, w, w).finished();
//  }

//  if (dim ==2)return (Vector(4) << 1.0,1.0,1.0,1.0).finished();
//  else return (Vector(6) << 1.0,1.0,1.0,1.0,1.0,1.0).finished();
  
}

double MM_multi::residual(const Vector error) const {

  int dim = 3;
  if (error.size()==4)dim = 2;
  const double e2_pair1 = error.head(dim).transpose()*error.head(dim);//error.head(dim).norm();//
  const double e2_pair2 = error.tail(dim).transpose()*error.tail(dim);//error.tail(dim).norm();//
  const double e2 = min(e2_pair1, e2_pair2);

  const double e4 = e2*e2;
  const double c2 = c_*c_;

  return (c2*e2 + c_*e4) / ((e2 + c_)*(e2 + c_));

//  if (e2 > log(2*3.141592/(b_*b_))/(a_*a_-b_*b_))
//  {
//    const double w = b_*b_/(a_*a_)*e2 + log(2*3.141592/(b_*b_));
//    return w;
//  }

//  return e2;
}

void MM_multi::print(const std::string &s="") const {
  std::cout << s << ": PairDCS (" << c_ << ")" << std::endl;
}

bool MM_multi::equals(const Base &expected, double tol) const {
  const MM_multi* p = dynamic_cast<const MM_multi*>(&expected);
  if (p == NULL) return false;
  return std::abs(c_ - p->c_) < tol;
}

MM_multi::shared_ptr MM_multi::Create(double c, const ReweightScheme reweight) {
  return shared_ptr(new MM_multi(c, reweight));
}

/* ************************************************************************* */
// MM
/* ************************************************************************* */
MM::MM(double a, double b, const ReweightScheme reweight)
  : Base(reweight), a_(a), b_(b) {
}

double MM::weight(double error) const {
  const double e2 = error*error;
  if (e2 > log(2*3.141592/(b_*b_))/(a_*a_-b_*b_))
  {
    const double w = b_*b_/(a_*a_) + log(2*3.141592/(b_*b_))/e2;
    return w;
  }

  return 1.0;
}

double MM::residual(const Vector error) const {
  // This is the simplified version of Eq 9 from (Agarwal13icra)
  // after you simplify and cancel terms.
  const double e2 = error.norm()*error.norm();
  if (e2 > log(2*3.141592/(b_*b_))/(a_*a_-b_*b_))
  {
    const double w = b_*b_/(a_*a_)*e2 + log(2*3.141592/(b_*b_));
    return w;
  }

  return e2;
}

void MM::print(const std::string &s="") const {
  std::cout << s << ": MM (" << a_ << ")" << std::endl;
}

bool MM::equals(const Base &expected, double tol) const {
  const MM* p = dynamic_cast<const MM*>(&expected);
  if (p == NULL) return false;
  return std::abs(a_ - p->a_) < tol;
}

MM::shared_ptr MM::Create(double a, double b, const ReweightScheme reweight) {
  return shared_ptr(new MM(a, b, reweight));
}

/* ************************************************************************* */
// DCE
/* ************************************************************************* */
DCE::DCE(double a, double b, const ReweightScheme reweight)
  : Base(reweight), a_(a), b_(b) {
}

double DCE::weight(double error) const {
  const double e2 = error*error;
  if (e2 > b_*b_/(a_*a_))
  {
    const double w = log(b_*b_/(a_*a_)*e2)/e2 +1.0/e2;
    return w;
  }

  return 1.0;
}

double DCE::residual(const Vector error) const {
  // This is the simplified version of Eq 9 from (Agarwal13icra)
  // after you simplify and cancel terms.
  const double e2 = error.norm()*error.norm();
  if (e2 > b_*b_/(a_*a_))
  {
    const double w = log(b_*b_/(a_*a_)*e2) + 1.0;
    return w;
  }

  return e2;
}

void DCE::print(const std::string &s="") const {
  std::cout << s << ": MM (" << a_ << ")" << std::endl;
}

bool DCE::equals(const Base &expected, double tol) const {
  const DCE* p = dynamic_cast<const DCE*>(&expected);
  if (p == NULL) return false;
  return std::abs(a_ - p->a_) < tol;
}

DCE::shared_ptr DCE::Create(double a, double b, const ReweightScheme reweight) {
  return shared_ptr(new DCE(a, b, reweight));
}
/* ************************************************************************* */
// L2WithDeadZone
/* ************************************************************************* */

L2WithDeadZone::L2WithDeadZone(double k, const ReweightScheme reweight)
 : Base(reweight), k_(k) {
  if (k_ <= 0) {
    throw runtime_error("mEstimator L2WithDeadZone takes only positive double in constructor.");
  }
}

double L2WithDeadZone::weight(double error) const {
  // note that this code is slightly uglier than residual, because there are three distinct
  // cases to handle (left of deadzone, deadzone, right of deadzone) instead of the two
  // cases (deadzone, non-deadzone) in residual.
  if (std::abs(error) <= k_) return 0.0;
  else if (error > k_) return (-k_+error)/error;
  else return (k_+error)/error;
}

double L2WithDeadZone::residual(const Vector error) const {
  const double abs_error = std::abs(error.norm());
  return (abs_error < k_) ? 0.0 : 0.5*(k_-abs_error)*(k_-abs_error);
}

void L2WithDeadZone::print(const std::string &s="") const {
  std::cout << s << ": L2WithDeadZone (" << k_ << ")" << std::endl;
}

bool L2WithDeadZone::equals(const Base &expected, double tol) const {
  const L2WithDeadZone* p = dynamic_cast<const L2WithDeadZone*>(&expected);
  if (p == NULL) return false;
  return std::abs(k_ - p->k_) < tol;
}

L2WithDeadZone::shared_ptr L2WithDeadZone::Create(double k, const ReweightScheme reweight) {
  return shared_ptr(new L2WithDeadZone(k, reweight));
}

} // namespace mEstimator
} // namespace noiseModel
} // gtsam
