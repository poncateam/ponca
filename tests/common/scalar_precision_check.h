/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.

 \author Nicolas Mellado nmellado0@gmail.com
*/


#ifndef SCALAR_PRECISION_CHECK_H
#define SCALAR_PRECISION_CHECK_H

#include <limits>

template <typename _lowPrecisionScalarT,
          typename _highPrecisionScalarT,
          bool autoCheck>
class ScalarPrecisionCheck;

namespace internal{

template <bool check>
struct _conditionnal_precision_assert{
    template <typename lst, typename hst, bool acheck>
    static inline const ScalarPrecisionCheck<lst, hst, acheck>&
    run(const ScalarPrecisionCheck<lst, hst, acheck>& s) { return s; }
};
template <>
struct _conditionnal_precision_assert<true>{
    template <typename lst, typename hst, bool acheck>
    static inline const ScalarPrecisionCheck<lst, hst, acheck>&
    run(const ScalarPrecisionCheck<lst, hst, acheck>& s)
    { assert(s.checkPrecision()); return s; }
};
}

/*!
 * \brief Scalar class running arithmetic operations on two scalar types to
 *  compare their precision.
 *
 * \tparam _lowPrecisionScalarT is aliased as LPScalar
 * \tparam _highPrecisionScalarT is aliased as HPScalar
 * \tparam autocheck enables precision check assertion
 *
 * This class support some std::math functions, right now pow, abs.
 *
 * \code
    typedef ScalarPrecisionCheck<float, double, false> ScalarPC;
    typedef ScalarPrecisionCheck<float, double, true> ScalarPCAuto;

    double epsilon = 10e-6;

    // epsilon value will be propagated from a and b
    ScalarPC a = 0.0000000000003; a.setEpsilon(epsilon);
    ScalarPC b = 17;              b.setEpsilon(epsilon);
    ScalarPC c = (a+b) * 0.4;
    c /= 0.4;
    c = abs(17 - c);
    cout << c.epsilon() << " >? " << c.precisionDelta() << endl;

    c = pow(c,3);
    cout << c.epsilon() << " >? " << c.precisionDelta() << endl;
    ScalarPCAuto ccheck = c;

    (void)c;
    \endcode
 */
template <typename _lowPrecisionScalarT  = float,
          typename _highPrecisionScalarT = double,
          bool autoCheck = true>
class ScalarPrecisionCheck {
public:
    typedef _lowPrecisionScalarT  LPScalar;
    typedef _highPrecisionScalarT HPScalar;
    typedef ScalarPrecisionCheck<LPScalar, HPScalar, autoCheck> ScalarPC;

    /*!
     * \fixme Not sure about how to properly set this epsilon value
     */
    template <typename RealScalar>
    inline ScalarPrecisionCheck(RealScalar v,
                                HPScalar epsilon = std::numeric_limits<HPScalar>::epsilon())
        : _lp(v), _hp(v), _epsilon(epsilon) {
        conditionnal_precision_assert::run(*this);
    }

    template <bool otherCheck> inline ScalarPrecisionCheck
    (const ScalarPrecisionCheck<LPScalar, HPScalar, otherCheck>& __rhs)
        : _lp(__rhs._lp), _hp(__rhs._hp), _epsilon(__rhs._epsilon) {
        conditionnal_precision_assert::run(*this);
    }

    inline ScalarPrecisionCheck(LPScalar lp,
                                HPScalar hp,
                                HPScalar epsilon = std::numeric_limits<HPScalar>::epsilon())
        : _lp(lp), _hp(hp), _epsilon(epsilon) {
        conditionnal_precision_assert::run(*this);
    }

    inline void setEpsilon(HPScalar epsilon) { _epsilon = epsilon; }
    inline HPScalar epsilon() const { return _epsilon; }

    inline HPScalar precisionDelta() const { return std::abs(HPScalar(_lp) - _hp); }
    inline bool checkPrecision() const { return ( precisionDelta() < _epsilon); }

    inline LPScalar lp() const { return _lp; }
    inline HPScalar hp() const { return _hp; }

    ////////////////////////////////////////////////////////////////////////////
    /// Pairwise operators with scalars
    ///
    template <typename RealScalar> inline ScalarPC operator* (RealScalar __rhs)
    { return ScalarPC(__rhs, this->_epsilon) * (*this); }

    template <typename RealScalar> inline ScalarPC operator/ (RealScalar __rhs)
    { return ScalarPC(__rhs, this->_epsilon) / (*this); }

    template <typename RealScalar> inline ScalarPC operator+ (RealScalar __rhs)
    { return ScalarPC(__rhs, this->_epsilon) + (*this); }

    template <typename RealScalar> inline ScalarPC operator- (RealScalar __rhs)
    { return ScalarPC(__rhs, this->_epsilon) - (*this); }

    template <typename RealScalar> inline ScalarPC& operator*= (RealScalar __rhs)
    { return *this = ScalarPC(__rhs, this->_epsilon) * (*this); }

    template <typename RealScalar> inline ScalarPC& operator/= (RealScalar __rhs)
    { return *this = ScalarPC(__rhs, this->_epsilon) / (*this); }

    template <typename RealScalar> inline ScalarPC& operator+= (RealScalar __rhs)
    { return *this = ScalarPC(__rhs, this->_epsilon) + (*this); }

    template <typename RealScalar> inline ScalarPC& operator-= (RealScalar __rhs)
    { return *this = ScalarPC(__rhs, this->_epsilon) - (*this); }

    ////////////////////////////////////////////////////////////////////////////
    /// Pairwise operators between ScalarPCs
    ///
    template <bool otherCheck> inline ScalarPC operator*
    (const ScalarPrecisionCheck<LPScalar, HPScalar, otherCheck>& __rhs)
    { return ScalarPC ( _lp * __rhs._lp, _hp * __rhs._hp, _epsilon ); }

    template <bool otherCheck> inline ScalarPC operator/
    (const ScalarPrecisionCheck<LPScalar, HPScalar, otherCheck>& __rhs)
    { return ScalarPC ( _lp / __rhs._lp, _hp / __rhs._hp, _epsilon ); }

    template <bool otherCheck> inline ScalarPC operator+
    (const ScalarPrecisionCheck<LPScalar, HPScalar, otherCheck>& __rhs)
    { return ScalarPC ( _lp + __rhs._lp, _hp + __rhs._hp, _epsilon ); }

    template <bool otherCheck> inline ScalarPC operator-
    (const ScalarPrecisionCheck<LPScalar, HPScalar, otherCheck>& __rhs)
    { return ScalarPC ( _lp - __rhs._lp, _hp - __rhs._hp, _epsilon ); }

    template <bool otherCheck>  inline ScalarPC& operator*=
    (const ScalarPrecisionCheck<LPScalar, HPScalar, otherCheck>& __rhs)
    { return *this = __rhs * (*this); }

    template <bool otherCheck>  inline ScalarPC& operator/=
    (const ScalarPrecisionCheck<LPScalar, HPScalar, otherCheck>& __rhs)
    { return *this = __rhs / (*this); }

    template <bool otherCheck>  inline ScalarPC& operator+=
    (const ScalarPrecisionCheck<LPScalar, HPScalar, otherCheck>& __rhs)
    { return *this = __rhs + (*this); }

    template <bool otherCheck>  inline ScalarPC& operator-=
    (const ScalarPrecisionCheck<LPScalar, HPScalar, otherCheck>& __rhs)
    { return *this = __rhs - (*this); }

private:
    LPScalar _lp;
    HPScalar _hp;
    HPScalar _epsilon;
    typedef internal::_conditionnal_precision_assert<autoCheck>
    conditionnal_precision_assert;

    // used to access to member fields easily in operators
    friend class ScalarPrecisionCheck<LPScalar, HPScalar, ! autoCheck>;
}; // class ScalarPrecisionCheck


template <typename RealScalar, typename lps, typename hps, bool acheck>
inline ScalarPrecisionCheck<lps,hps,acheck> operator*  (RealScalar __lhs,
                            const ScalarPrecisionCheck<lps,hps,acheck>& __rhs)
{ return ScalarPrecisionCheck<lps,hps,acheck>(__lhs, __rhs.epsilon()) * __rhs; }

template <typename RealScalar, typename lps, typename hps, bool acheck>
inline ScalarPrecisionCheck<lps,hps,acheck> operator/  (RealScalar __lhs,
                            const ScalarPrecisionCheck<lps,hps,acheck>& __rhs)
{ return ScalarPrecisionCheck<lps,hps,acheck>(__lhs, __rhs.epsilon()) / __rhs; }

template <typename RealScalar, typename lps, typename hps, bool acheck>
inline ScalarPrecisionCheck<lps,hps,acheck> operator+  (RealScalar __lhs,
                            const ScalarPrecisionCheck<lps,hps,acheck>& __rhs)
{ return ScalarPrecisionCheck<lps,hps,acheck>(__lhs, __rhs.epsilon()) + __rhs; }

template <typename RealScalar, typename lps, typename hps, bool acheck>
inline ScalarPrecisionCheck<lps,hps,acheck> operator-  (RealScalar __lhs,
                            const ScalarPrecisionCheck<lps,hps,acheck>& __rhs)
{ return ScalarPrecisionCheck<lps,hps,acheck>(__lhs, __rhs.epsilon()) - __rhs; }

// create definitions for common std math functions
namespace std{
template <typename lps, typename hps, bool acheck>
inline ScalarPrecisionCheck<lps,hps,acheck> abs(
        const ScalarPrecisionCheck<lps,hps,acheck>& spc){
    return ScalarPrecisionCheck<lps,hps,acheck>(std::abs(spc.lp()),
                                                std::abs(spc.hp()),
                                                spc.epsilon());
}

template <typename lps, typename hps, bool acheck, typename PowScalar>
inline ScalarPrecisionCheck<lps,hps,acheck> pow(
        const ScalarPrecisionCheck<lps,hps,acheck>& spc, PowScalar p){
    return ScalarPrecisionCheck<lps,hps,acheck>(std::pow(spc.lp(), p),
                                                std::pow(spc.hp(), p),
                                                spc.epsilon());
}

} // namespace std

#endif // SCALAR_PRECISION_CHECK_H
