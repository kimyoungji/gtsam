/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation, 
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file Expression.h
 * @date September 18, 2014
 * @author Frank Dellaert
 * @author Paul Furgale
 * @brief Expressions for Block Automatic Differentiation
 */

#include <gtsam/nonlinear/NonlinearFactor.h>
#include <gtsam/inference/Key.h>

#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>
#include <boost/bind.hpp>

namespace gtsam {

///-----------------------------------------------------------------------------
/// Expression node. The superclass for objects that do the heavy lifting
/// An Expression<T> has a pointer to an ExpressionNode<T> underneath
/// allowing Expressions to have polymorphic behaviour even though they
/// are passed by value. This is the same way boost::function works.
/// http://loki-lib.sourceforge.net/html/a00652.html
template<class T>
class ExpressionNode {
protected:
  ExpressionNode() {
  }
public:
  virtual ~ExpressionNode() {
  }

  /// Return keys that play in this expression as a set
  virtual std::set<Key> keys() const = 0;

  /// Return value and optional derivatives
  virtual T value(const Values& values,
      boost::optional<std::map<Key, Matrix>&> = boost::none) const = 0;
};

template<typename T>
class Expression;

/// Constant Expression
template<class T>
class ConstantExpression: public ExpressionNode<T> {

  T value_;

  /// Constructor with a value, yielding a constant
  ConstantExpression(const T& value) :
      value_(value) {
  }

  friend class Expression<T> ;

public:

  virtual ~ConstantExpression() {
  }

  /// Return keys that play in this expression, i.e., the empty set
  virtual std::set<Key> keys() const {
    std::set<Key> keys;
    return keys;
  }

  /// Return value and optional derivatives
  virtual T value(const Values& values,
      boost::optional<std::map<Key, Matrix>&> jacobians = boost::none) const {
    return value_;
  }
};

//-----------------------------------------------------------------------------
/// Leaf Expression
template<class T>
class LeafExpression: public ExpressionNode<T> {

  Key key_;

  /// Constructor with a single key
  LeafExpression(Key key) :
      key_(key) {
  }

  friend class Expression<T> ;

public:

  virtual ~LeafExpression() {
  }

  /// Return keys that play in this expression
  virtual std::set<Key> keys() const {
    std::set<Key> keys;
    keys.insert(key_);
    return keys;
  }

  /// Return value and optional derivatives
  virtual T value(const Values& values,
      boost::optional<std::map<Key, Matrix>&> jacobians = boost::none) const {
    const T& value = values.at<T>(key_);
    if (jacobians) {
      std::map<Key, Matrix>::iterator it = jacobians->find(key_);
      if (it != jacobians->end()) {
        it->second += Eigen::MatrixXd::Identity(value.dim(), value.dim());
      } else {
        (*jacobians)[key_] = Eigen::MatrixXd::Identity(value.dim(),
            value.dim());
      }
    }
    return value;
  }

};

//-----------------------------------------------------------------------------
/// Unary Expression
template<class T, class E>
class UnaryExpression: public ExpressionNode<T> {

public:

  typedef boost::function<T(const E&, boost::optional<Matrix&>)> function;

private:

  boost::shared_ptr<ExpressionNode<E> > expression_;
  function f_;

  /// Constructor with a unary function f, and input argument e
  UnaryExpression(function f, const Expression<E>& e) :
      expression_(e.root()), f_(f) {
  }

  friend class Expression<T> ;

public:

  virtual ~UnaryExpression() {
  }

  /// Return keys that play in this expression
  virtual std::set<Key> keys() const {
    return expression_->keys();
  }

  /// Return value and optional derivatives
  virtual T value(const Values& values,
      boost::optional<std::map<Key, Matrix>&> jacobians = boost::none) const {

    T value;
    if (jacobians) {
      Eigen::MatrixXd H;
      value = f_(expression_->value(values, jacobians), H);
      std::map<Key, Matrix>::iterator it = jacobians->begin();
      for (; it != jacobians->end(); ++it) {
        it->second = H * it->second;
      }
    } else {
      value = f_(expression_->value(values), boost::none);
    }
    return value;
  }

};

//-----------------------------------------------------------------------------
/// Binary Expression

template<class T, class E1, class E2>
class BinaryExpression: public ExpressionNode<T> {

public:

  typedef boost::function<
      T(const E1&, const E2&, boost::optional<Matrix&>,
          boost::optional<Matrix&>)> function;
private:

  boost::shared_ptr<ExpressionNode<E1> > expression1_;
  boost::shared_ptr<ExpressionNode<E2> > expression2_;
  function f_;

  /// Constructor with a binary function f, and two input arguments
  BinaryExpression(function f, //
      const Expression<E1>& e1, const Expression<E2>& e2) :
      expression1_(e1.root()), expression2_(e2.root()), f_(f) {
  }

  friend class Expression<T> ;

public:

  virtual ~BinaryExpression() {
  }

  /// Return keys that play in this expression
  virtual std::set<Key> keys() const {
    std::set<Key> keys1 = expression1_->keys();
    std::set<Key> keys2 = expression2_->keys();
    keys1.insert(keys2.begin(), keys2.end());
    return keys1;
  }

  /// Return value and optional derivatives
  virtual T value(const Values& values,
      boost::optional<std::map<Key, Matrix>&> jacobians = boost::none) const {
    T val;
    if (jacobians) {
      std::map<Key, Matrix> terms1;
      std::map<Key, Matrix> terms2;
      Matrix H1, H2;
      val = f_(expression1_->value(values, terms1),
          expression2_->value(values, terms2), H1, H2);
      // TODO: both Jacobians and terms are sorted. There must be a simple
      //       but fast algorithm that does this.
      typedef std::pair<Key, Matrix> Pair;
      BOOST_FOREACH(const Pair& term, terms1) {
        std::map<Key, Matrix>::iterator it = jacobians->find(term.first);
        if (it != jacobians->end()) {
          it->second += H1 * term.second;
        } else {
          (*jacobians)[term.first] = H1 * term.second;
        }
      }
      BOOST_FOREACH(const Pair& term, terms2) {
        std::map<Key, Matrix>::iterator it = jacobians->find(term.first);
        if (it != jacobians->end()) {
          it->second += H2 * term.second;
        } else {
          (*jacobians)[term.first] = H2 * term.second;
        }
      }
    } else {
      val = f_(expression1_->value(values), expression2_->value(values),
          boost::none, boost::none);
    }
    return val;
  }

};

/**
 * Expression class that supports automatic differentiation
 */
template<typename T>
class Expression {
public:

  // Construct a constant expression
  Expression(const T& value) :
      root_(new ConstantExpression<T>(value)) {
  }

  // Construct a leaf expression
  Expression(const Key& key) :
      root_(new LeafExpression<T>(key)) {
  }

  /// Construct a unary expression
  template<typename E>
  Expression(typename UnaryExpression<T, E>::function f,
      const Expression<E>& expression) {
    // TODO Assert that root of expression is not null.
    root_.reset(new UnaryExpression<T, E>(f, expression));
  }

  /// Construct a binary expression
  template<typename E1, typename E2>
  Expression(typename BinaryExpression<T, E1, E2>::function f,
      const Expression<E1>& expression1, const Expression<E2>& expression2) {
    // TODO Assert that root of expressions 1 and 2 are not null.
    root_.reset(new BinaryExpression<T, E1, E2>(f, expression1, expression2));
  }

  /// Return keys that play in this expression
  std::set<Key> keys() const {
    return root_->keys();
  }

  /// Return value and optional derivatives
  T value(const Values& values,
      boost::optional<std::map<Key, Matrix>&> jacobians = boost::none) const {
    return root_->value(values, jacobians);
  }

  const boost::shared_ptr<ExpressionNode<T> >& root() const {
    return root_;
  }
private:
  boost::shared_ptr<ExpressionNode<T> > root_;
};

// http://stackoverflow.com/questions/16260445/boost-bind-to-operator
template<class T>
struct apply_compose {
  typedef T result_type;
  T operator()(const T& x, const T& y, boost::optional<Matrix&> H1,
      boost::optional<Matrix&> H2) const {
    return x.compose(y, H1, H2);
  }
};

/// Construct a product expression, assumes T::compose(T) -> T
template<typename T>
Expression<T> operator*(const Expression<T>& expression1,
    const Expression<T>& expression2) {
  return Expression<T>(boost::bind(apply_compose<T>(), _1, _2, _3, _4),
      expression1, expression2);
}

//-----------------------------------------------------------------------------
/// AD Factor
template<class T>
class BADFactor: NonlinearFactor {

  const T measurement_;
  const Expression<T> expression_;

  /// get value from expression and calculate error with respect to measurement
  Vector unwhitenedError(const Values& values) const {
    const T& value = expression_.value(values);
    return value.localCoordinates(measurement_);
  }

public:

  /// Constructor
  BADFactor(const T& measurement, const Expression<T>& expression) :
      measurement_(measurement), expression_(expression) {
  }
  /// Constructor
  BADFactor(const T& measurement, const ExpressionNode<T>& expression) :
      measurement_(measurement), expression_(expression) {
  }
  /**
   * Calculate the error of the factor.
   * This is the log-likelihood, e.g. \f$ 0.5(h(x)-z)^2/\sigma^2 \f$ in case of Gaussian.
   * In this class, we take the raw prediction error \f$ h(x)-z \f$, ask the noise model
   * to transform it to \f$ (h(x)-z)^2/\sigma^2 \f$, and then multiply by 0.5.
   */
  virtual double error(const Values& values) const {
    if (this->active(values)) {
      const Vector e = unwhitenedError(values);
      return 0.5 * e.squaredNorm();
    } else {
      return 0.0;
    }
  }

  /// get the dimension of the factor (number of rows on linearization)
  size_t dim() const {
    return 0;
  }

  /// linearize to a GaussianFactor
  boost::shared_ptr<GaussianFactor> linearize(const Values& values) const {
    // We will construct an n-ary factor below, where  terms is a container whose
    // value type is std::pair<Key, Matrix>, specifying the
    // collection of keys and matrices making up the factor.
    std::map<Key, Matrix> terms;
    expression_.value(values, terms);
    Vector b = unwhitenedError(values);
    SharedDiagonal model = SharedDiagonal();
    return boost::shared_ptr<JacobianFactor>(
        new JacobianFactor(terms, b, model));
  }

};
}
