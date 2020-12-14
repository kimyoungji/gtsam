/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 *  @file   GPSFactor.h
 *  @author Frank Dellaert
 *  @brief  Header file for GPS factor
 *  @date   January 22, 2014
 **/
#pragma once

#include <gtsam/nonlinear/NonlinearFactor.h>
#include <gtsam/navigation/NavState.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam/geometry/Pose2.h>

namespace gtsam {

/**
 * Prior on position in a Cartesian frame.
 * Possibilities include:
 *   ENU: East-North-Up navigation frame at some local origin
 *   NED: North-East-Down navigation frame at some local origin
 *   ECEF: Earth-centered Earth-fixed, origin at Earth's center
 * See Farrell08book or e.g. http://www.dirsig.org/docs/new/coordinates.html
 * @addtogroup Navigation
 */
class GTSAM_EXPORT GPSFactor: public NoiseModelFactor1<Pose3> {

private:

  typedef NoiseModelFactor1<Pose3> Base;

  Point3 nT_; ///< Position measurement in cartesian coordinates

public:

  /// shorthand for a smart pointer to a factor
  typedef boost::shared_ptr<GPSFactor> shared_ptr;

  /// Typedef to this class
  typedef GPSFactor This;

  /** default constructor - only use for serialization */
  GPSFactor(): nT_(0, 0, 0) {}

  virtual ~GPSFactor() {}

  /**
   * @brief Constructor from a measurement in a Cartesian frame.
   * Use GeographicLib to convert from geographic (latitude and longitude) coordinates
   * @param key of the Pose3 variable that will be constrained
   * @param gpsIn measurement already in correct coordinates
   * @param model Gaussian noise model
   */
  GPSFactor(Key key, const Point3& gpsIn, const SharedNoiseModel& model) :
      Base(model, key), nT_(gpsIn) {
  }

  /// @return a deep copy of this factor
  virtual gtsam::NonlinearFactor::shared_ptr clone() const {
    return boost::static_pointer_cast<gtsam::NonlinearFactor>(
        gtsam::NonlinearFactor::shared_ptr(new This(*this)));
  }

  /// print
  virtual void print(const std::string& s, const KeyFormatter& keyFormatter =
      DefaultKeyFormatter) const;

  /// equals
  virtual bool equals(const NonlinearFactor& expected, double tol = 1e-9) const;

  /// vector of errors
  Vector evaluateError(const Pose3& p,
      boost::optional<Matrix&> H = boost::none) const;

  inline const Point3 & measurementIn() const {
    return nT_;
  }

  /**
   *  Convenience function to estimate state at time t, given two GPS
   *  readings (in local NED Cartesian frame) bracketing t
   *  Assumes roll is zero, calculates yaw and pitch from NED1->NED2 vector.
   */
  static std::pair<Pose3, Vector3> EstimateState(double t1, const Point3& NED1,
      double t2, const Point3& NED2, double timestamp);

private:

  /// Serialization function
  friend class boost::serialization::access;
  template<class ARCHIVE>
  void serialize(ARCHIVE & ar, const unsigned int /*version*/) {
    ar
        & boost::serialization::make_nvp("NoiseModelFactor1",
            boost::serialization::base_object<Base>(*this));
    ar & BOOST_SERIALIZATION_NVP(nT_);
  }
};

/**
 * Version of GPSFactor for NavState
 * @addtogroup Navigation
 */
class GTSAM_EXPORT GPSFactor2: public NoiseModelFactor1<NavState> {

private:

  typedef NoiseModelFactor1<NavState> Base;

  Point3 nT_; ///< Position measurement in cartesian coordinates

public:

  /// shorthand for a smart pointer to a factor
  typedef boost::shared_ptr<GPSFactor2> shared_ptr;

  /// Typedef to this class
  typedef GPSFactor2 This;

  /// default constructor - only use for serialization
  GPSFactor2():nT_(0, 0, 0) {}

  virtual ~GPSFactor2() {}

  /// Constructor from a measurement in a Cartesian frame.
  GPSFactor2(Key key, const Point3& gpsIn, const SharedNoiseModel& model) :
      Base(model, key), nT_(gpsIn) {
  }

  /// @return a deep copy of this factor
  virtual gtsam::NonlinearFactor::shared_ptr clone() const {
    return boost::static_pointer_cast<gtsam::NonlinearFactor>(
        gtsam::NonlinearFactor::shared_ptr(new This(*this)));
  }

  /// print
  virtual void print(const std::string& s, const KeyFormatter& keyFormatter =
      DefaultKeyFormatter) const;

  /// equals
  virtual bool equals(const NonlinearFactor& expected, double tol = 1e-9) const;

  /// vector of errors
  Vector evaluateError(const NavState& p,
      boost::optional<Matrix&> H = boost::none) const;

  inline const Point3 & measurementIn() const {
    return nT_;
  }

private:

  /// Serialization function
  friend class boost::serialization::access;
  template<class ARCHIVE>
  void serialize(ARCHIVE & ar, const unsigned int /*version*/) {
    ar
        & boost::serialization::make_nvp("NoiseModelFactor1",
            boost::serialization::base_object<Base>(*this));
    ar & BOOST_SERIALIZATION_NVP(nT_);
  }
};

class GPSFactorPair3D: public NoiseModelFactor1<Pose3> {

  // The factor will hold a measurement consisting of an (X,Y) location
  // We could this with a Point2 but here we just use two doubles
  Point3 nT1_, nT2_;

public:
  /// shorthand for a smart pointer to a factor
  typedef boost::shared_ptr<GPSFactorPair3D> shared_ptr;

  // The constructor requires the variable key, the (X, Y) measurement value, and the noise model
  GPSFactorPair3D(Key j, const Point3& gpsIn, const Point3& vpsIn, const SharedNoiseModel& model):
    NoiseModelFactor1<Pose3>(model, j), nT1_(gpsIn), nT2_(vpsIn){}

  virtual ~GPSFactorPair3D() {}

  /// vector of errors
  Vector evaluateError(const Pose3& p, boost::optional<Matrix&> H = boost::none) const {
  Vector res = (Vector(6) << p.translation(H) -nT1_, p.translation(H) -nT2_).finished();
//  std::cout<<res<<std::endl;
//  std::cout<<*H<<std::endl;
  if (H){
    (*H) = (Matrix(6,6) << *H, *H).finished();
//    (*H) = (Matrix(6,6) << 1.0,0.0,0.0,0.0,0.0,0.0, 0.0,1.0,0.0,0.0,0.0,0.0, 0.0,0.0,1.0,0.0,0.0,0.0, 1.0,0.0,0.0,0.0,0.0,0.0, 0.0,1.0,0.0,0.0,0.0,0.0, 0.0,0.0,1.0,0.0,0.0,0.0).finished();
  }
  return res;
  }

  // The second is a 'clone' function that allows the factor to be copied. Under most
  // circumstances, the following code that employs the default copy constructor should
  // work fine.
  virtual gtsam::NonlinearFactor::shared_ptr clone() const {
    return boost::static_pointer_cast<gtsam::NonlinearFactor>(
        gtsam::NonlinearFactor::shared_ptr(new GPSFactorPair3D(*this))); }

  // Additionally, we encourage you the use of unit testing your custom factors,
  // (as all GTSAM factors are), in which you would need an equals and print, to satisfy the
  // GTSAM_CONCEPT_TESTABLE_INST(T) defined in Testable.h, but these are not needed below.

}; // GPSFactorPair3D

//class GPSFactor3D: public NoiseModelFactor1<Pose3> {

//  // The factor will hold a measurement consisting of an (X,Y) location
//  // We could this with a Point2 but here we just use two doubles
//  double mx_, my_, mz_;

//public:
//  /// shorthand for a smart pointer to a factor
//  typedef boost::shared_ptr<GPSFactor3D> shared_ptr;

//  // The constructor requires the variable key, the (X, Y) measurement value, and the noise model
//  GPSFactor3D(Key j, double x, double y, double z, const SharedNoiseModel& model):
//    NoiseModelFactor1<Pose3>(model, j), mx_(x), my_(y), mz_(z) {}

//  virtual ~GPSFactor3D() {}

//  // Using the NoiseModelFactor1 base class there are two functions that must be overridden.
//  // The first is the 'evaluateError' function. This function implements the desired measurement
//  // function, returning a vector of errors when evaluated at the provided variable value. It
//  // must also calculate the Jacobians for this measurement function, if requested.
//  Vector evaluateError(const Pose3& q, boost::optional<Matrix&> H = boost::none) const
//  {
//    // The measurement function for a GPS-like measurement is simple:
//    // error_x = pose.x - measurement.x
//    // error_y = pose.y - measurement.y
//    // Consequently, the Jacobians are:
//    // [ derror_x/dx  derror_x/dy  derror_x/dtheta ] = [1 0 0]
//    // [ derror_y/dx  derror_y/dy  derror_y/dtheta ] = [0 1 0]
//    if (H){
//        (*H) = (Matrix(3,6) << 0.0,0.0,0.0,1.0,0.0,0.0, 0.0,0.0,0.0,0.0,1.0,0.0, 0.0,0.0,0.0,0.0,0.0,1.0).finished();
////        (*H1) = (Matrix(2,3) << 1.0,0.0,0.0, 0.0,1.0,0.0).finished();
////        (*H2) = (Matrix(2,3) << 1.0,0.0,0.0, 0.0,1.0,0.0).finished();
//    }
//    return (Vector(3) << q.x() - mx_, q.y() - my_, q.z() - mz_).finished();
//  }

//  // The second is a 'clone' function that allows the factor to be copied. Under most
//  // circumstances, the following code that employs the default copy constructor should
//  // work fine.
//  virtual gtsam::NonlinearFactor::shared_ptr clone() const {
//    return boost::static_pointer_cast<gtsam::NonlinearFactor>(
//        gtsam::NonlinearFactor::shared_ptr(new GPSFactor3D(*this))); }

//  // Additionally, we encourage you the use of unit testing your custom factors,
//  // (as all GTSAM factors are), in which you would need an equals and print, to satisfy the
//  // GTSAM_CONCEPT_TESTABLE_INST(T) defined in Testable.h, but these are not needed below.

//}; // GPSFactor3D

class GPSFactorPair: public NoiseModelFactor1<Pose2> {

  // The factor will hold a measurement consisting of an (X,Y) location
  // We could this with a Point2 but here we just use two doubles
  double mx1_, my1_, mx2_, my2_;

public:
  /// shorthand for a smart pointer to a factor
  typedef boost::shared_ptr<GPSFactorPair> shared_ptr;

  // The constructor requires the variable key, the (X, Y) measurement value, and the noise model
  GPSFactorPair(Key j, double x1, double y1, double x2, double y2, const SharedNoiseModel& model):
    NoiseModelFactor1<Pose2>(model, j), mx1_(x1), my1_(y1), mx2_(x2), my2_(y2){}

  virtual ~GPSFactorPair() {}

  // Using the NoiseModelFactor1 base class there are two functions that must be overridden.
  // The first is the 'evaluateError' function. This function implements the desired measurement
  // function, returning a vector of errors when evaluated at the provided variable value. It
  // must also calculate the Jacobians for this measurement function, if requested.
  Vector evaluateError(const Pose2& q, boost::optional<Matrix&> H = boost::none) const
  {
    // The measurement function for a GPS-like measurement is simple:
    // error_x = pose.x - measurement.x
    // error_y = pose.y - measurement.y
    // Consequently, the Jacobians are:
    // [ derror_x/dx  derror_x/dy  derror_x/dtheta ] = [1 0 0]
    // [ derror_y/dx  derror_y/dy  derror_y/dtheta ] = [0 1 0]
    if (H){
        (*H) = (Matrix(4,3) << 1.0,0.0,0.0, 0.0,1.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0).finished();
//        (*H1) = (Matrix(2,3) << 1.0,0.0,0.0, 0.0,1.0,0.0).finished();
//        (*H2) = (Matrix(2,3) << 1.0,0.0,0.0, 0.0,1.0,0.0).finished();
    }
    return (Vector(4) << q.x() - mx1_, q.y() - my1_, q.x() - mx2_, q.y() - my2_).finished();
  }

  // The second is a 'clone' function that allows the factor to be copied. Under most
  // circumstances, the following code that employs the default copy constructor should
  // work fine.
  virtual gtsam::NonlinearFactor::shared_ptr clone() const {
    return boost::static_pointer_cast<gtsam::NonlinearFactor>(
        gtsam::NonlinearFactor::shared_ptr(new GPSFactorPair(*this))); }

  // Additionally, we encourage you the use of unit testing your custom factors,
  // (as all GTSAM factors are), in which you would need an equals and print, to satisfy the
  // GTSAM_CONCEPT_TESTABLE_INST(T) defined in Testable.h, but these are not needed below.

}; // GPSFactorPair

class GPSFactor2D: public NoiseModelFactor1<Pose2> {

  // The factor will hold a measurement consisting of an (X,Y) location
  // We could this with a Point2 but here we just use two doubles
  double mx_, my_;

public:
  /// shorthand for a smart pointer to a factor
  typedef boost::shared_ptr<GPSFactor2D> shared_ptr;

  // The constructor requires the variable key, the (X, Y) measurement value, and the noise model
  GPSFactor2D(Key j, double x, double y, const SharedNoiseModel& model):
    NoiseModelFactor1<Pose2>(model, j), mx_(x), my_(y) {}

  virtual ~GPSFactor2D() {}

  // Using the NoiseModelFactor1 base class there are two functions that must be overridden.
  // The first is the 'evaluateError' function. This function implements the desired measurement
  // function, returning a vector of errors when evaluated at the provided variable value. It
  // must also calculate the Jacobians for this measurement function, if requested.
  Vector evaluateError(const Pose2& q, boost::optional<Matrix&> H = boost::none) const
  {
    // The measurement function for a GPS-like measurement is simple:
    // error_x = pose.x - measurement.x
    // error_y = pose.y - measurement.y
    // Consequently, the Jacobians are:
    // [ derror_x/dx  derror_x/dy  derror_x/dtheta ] = [1 0 0]
    // [ derror_y/dx  derror_y/dy  derror_y/dtheta ] = [0 1 0]
    if (H){
        (*H) = (Matrix(2,3) << 1.0,0.0,0.0, 0.0,1.0,0.0).finished();
//        (*H1) = (Matrix(2,3) << 1.0,0.0,0.0, 0.0,1.0,0.0).finished();
//        (*H2) = (Matrix(2,3) << 1.0,0.0,0.0, 0.0,1.0,0.0).finished();
    }
    return (Vector(2) << q.x() - mx_, q.y() - my_).finished();
  }

  // The second is a 'clone' function that allows the factor to be copied. Under most
  // circumstances, the following code that employs the default copy constructor should
  // work fine.
  virtual gtsam::NonlinearFactor::shared_ptr clone() const {
    return boost::static_pointer_cast<gtsam::NonlinearFactor>(
        gtsam::NonlinearFactor::shared_ptr(new GPSFactor2D(*this))); }

  // Additionally, we encourage you the use of unit testing your custom factors,
  // (as all GTSAM factors are), in which you would need an equals and print, to satisfy the
  // GTSAM_CONCEPT_TESTABLE_INST(T) defined in Testable.h, but these are not needed below.

}; // GPSFactor2D



} /// namespace gtsam
