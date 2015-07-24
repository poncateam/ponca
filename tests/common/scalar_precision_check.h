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

template <bool check>
struct _conditionnal_AND_assert{
    static inline bool run (bool b0, bool b1){ return b0 && b1; }
};

template <>
struct _conditionnal_AND_assert<true>{
    static inline bool run(bool lp, bool hp){ assert (lp && hp); return lp && hp; }
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
    typedef ScalarPrecisionCheck<float, double> ScalarPC;
    ScalarPC::epsilon = 10e-4;

    ScalarPC a = 2;
    ScalarPC b = 10e-6;
    ScalarPC c = pow(a, 3) * b;

    std::cout << c.lp() << " " << c.hp() << std::endl;
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
    inline ScalarPrecisionCheck(RealScalar v)
        : _lp(v), _hp(v) {
        precision_assert::run(*this);
    }

    template <bool otherCheck> inline ScalarPrecisionCheck
    (const ScalarPrecisionCheck<LPScalar, HPScalar, otherCheck>& __rhs)
        : _lp(__rhs._lp), _hp(__rhs._hp) {
        precision_assert::run(*this);
    }

    inline ScalarPrecisionCheck(LPScalar lp,
                                HPScalar hp)
        : _lp(lp), _hp(hp) {
        precision_assert::run(*this);
    }

    inline HPScalar precisionDelta() const { return std::abs(HPScalar(_lp) - _hp); }
    inline bool checkPrecision() const { return ( precisionDelta() < epsilon); }

    inline LPScalar lp() const { return _lp; }
    inline HPScalar hp() const { return _hp; }

    ////////////////////////////////////////////////////////////////////////////
    /// Pairwise operators with scalars
    ///
    template <typename RealScalar> inline ScalarPC operator* (RealScalar __rhs) const
    { return ScalarPC(__rhs) * (*this); }

    template <typename RealScalar> inline ScalarPC operator/ (RealScalar __rhs) const
    { return ScalarPC(__rhs) / (*this); }

    template <typename RealScalar> inline ScalarPC operator+ (RealScalar __rhs) const
    { return ScalarPC(__rhs) + (*this); }

    template <typename RealScalar> inline ScalarPC operator- (RealScalar __rhs) const
    { return ScalarPC(__rhs) - (*this); }

    template <typename RealScalar> inline ScalarPC& operator*= (RealScalar __rhs)
    { return *this = ScalarPC(__rhs) * (*this); }

    template <typename RealScalar> inline ScalarPC& operator/= (RealScalar __rhs)
    { return *this = ScalarPC(__rhs) / (*this); }

    template <typename RealScalar> inline ScalarPC& operator+= (RealScalar __rhs)
    { return *this = ScalarPC(__rhs) + (*this); }

    template <typename RealScalar> inline ScalarPC& operator-= (RealScalar __rhs)
    { return *this = ScalarPC(__rhs) - (*this); }

    template <typename RealScalar> inline bool operator< (RealScalar __rhs) const
    { return AND_assert::run( _lp < __rhs, _hp < __rhs ); }

    template <typename RealScalar> inline bool operator> (RealScalar __rhs) const
    { return AND_assert::run( _lp > __rhs, _hp > __rhs ); }

    template <typename RealScalar> inline bool operator== (RealScalar __rhs) const
    { return AND_assert::run( _lp == __rhs, _hp == __rhs ); }

    template <typename RealScalar> inline bool operator!= (RealScalar __rhs) const
    { return AND_assert::run( _lp != __rhs, _hp != __rhs ); }

    ////////////////////////////////////////////////////////////////////////////
    /// Pairwise operators between ScalarPCs
    ///
    template <bool otherCheck> inline ScalarPC operator*
    (const ScalarPrecisionCheck<LPScalar, HPScalar, otherCheck>& __rhs) const
    { return ScalarPC ( _lp * __rhs._lp, _hp * __rhs._hp ); }

    template <bool otherCheck> inline ScalarPC operator/
    (const ScalarPrecisionCheck<LPScalar, HPScalar, otherCheck>& __rhs) const
    { return ScalarPC ( _lp / __rhs._lp, _hp / __rhs._hp ); }

    template <bool otherCheck> inline ScalarPC operator+
    (const ScalarPrecisionCheck<LPScalar, HPScalar, otherCheck>& __rhs) const
    { return ScalarPC ( _lp + __rhs._lp, _hp + __rhs._hp ); }

    template <bool otherCheck> inline ScalarPC operator-
    (const ScalarPrecisionCheck<LPScalar, HPScalar, otherCheck>& __rhs) const
    { return ScalarPC ( _lp - __rhs._lp, _hp - __rhs._hp ); }

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

    template <bool otherCheck> inline bool operator<
    (const ScalarPrecisionCheck<LPScalar, HPScalar, otherCheck>& __rhs) const
    { return AND_assert::run( _lp < __rhs._lp, _hp < __rhs._hp ); }


    template <bool otherCheck> inline bool operator>
    (const ScalarPrecisionCheck<LPScalar, HPScalar, otherCheck>& __rhs) const
    { return AND_assert::run( _lp > __rhs._lp, _hp > __rhs._hp ); }

    template <bool otherCheck> inline bool operator==
    (const ScalarPrecisionCheck<LPScalar, HPScalar, otherCheck>& __rhs) const
    { return AND_assert::run( _lp == __rhs._lp, _hp == __rhs._hp ); }


    template <bool otherCheck> inline bool operator!=
    (const ScalarPrecisionCheck<LPScalar, HPScalar, otherCheck>& __rhs) const
    { return AND_assert::run( _lp != __rhs._lp, _hp != __rhs._hp ); }

public:
    static HPScalar epsilon;

private:
    LPScalar _lp;
    HPScalar _hp;
    typedef ::internal::_conditionnal_precision_assert<autoCheck> precision_assert;
    typedef ::internal::_conditionnal_AND_assert<autoCheck> AND_assert;

    // used to access to member fields easily in operators
    friend class ScalarPrecisionCheck<LPScalar, HPScalar, ! autoCheck>;
}; // class ScalarPrecisionCheck

template <typename lps, typename hps, bool c>
typename ScalarPrecisionCheck<lps,hps,c>::HPScalar
ScalarPrecisionCheck<lps,hps,c>::epsilon =
         std::numeric_limits<typename ScalarPrecisionCheck<lps,hps,c>::HPScalar>::epsilon();


template <typename RealScalar, typename lps, typename hps, bool acheck>
inline ScalarPrecisionCheck<lps,hps,acheck> operator*  (RealScalar __lhs,
                            const ScalarPrecisionCheck<lps,hps,acheck>& __rhs)
{ return ScalarPrecisionCheck<lps,hps,acheck>(__lhs) * __rhs; }

template <typename RealScalar, typename lps, typename hps, bool acheck>
inline ScalarPrecisionCheck<lps,hps,acheck> operator/  (RealScalar __lhs,
                            const ScalarPrecisionCheck<lps,hps,acheck>& __rhs)
{ return ScalarPrecisionCheck<lps,hps,acheck>(__lhs) / __rhs; }

template <typename RealScalar, typename lps, typename hps, bool acheck>
inline ScalarPrecisionCheck<lps,hps,acheck> operator+  (RealScalar __lhs,
                            const ScalarPrecisionCheck<lps,hps,acheck>& __rhs)
{ return ScalarPrecisionCheck<lps,hps,acheck>(__lhs) + __rhs; }

template <typename RealScalar, typename lps, typename hps, bool acheck>
inline ScalarPrecisionCheck<lps,hps,acheck> operator-  (RealScalar __lhs,
                            const ScalarPrecisionCheck<lps,hps,acheck>& __rhs)
{ return ScalarPrecisionCheck<lps,hps,acheck>(__lhs) - __rhs; }

template <typename RealScalar, typename lps, typename hps, bool acheck>
inline bool operator<  (RealScalar __lhs,
                        const ScalarPrecisionCheck<lps,hps,acheck>& __rhs)
{ return __rhs > __lhs; }

template <typename RealScalar, typename lps, typename hps, bool acheck>
inline bool operator>  (RealScalar __lhs,
                        const ScalarPrecisionCheck<lps,hps,acheck>& __rhs)
{ return __rhs.operator< (__lhs); }

template <typename RealScalar, typename lps, typename hps, bool acheck>
inline bool operator==  (RealScalar __lhs,
                        const ScalarPrecisionCheck<lps,hps,acheck>& __rhs)
{ return __rhs == __lhs; }

template <typename RealScalar, typename lps, typename hps, bool acheck>
inline bool operator!=  (RealScalar __lhs,
                        const ScalarPrecisionCheck<lps,hps,acheck>& __rhs)
{ return __rhs != __lhs; }

// create definitions for common std math functions
namespace std{
template <typename lps, typename hps, bool acheck>
inline ScalarPrecisionCheck<lps,hps,acheck> abs(
        const ScalarPrecisionCheck<lps,hps,acheck>& spc){
    return ScalarPrecisionCheck<lps,hps,acheck>(std::abs(spc.lp()),
                                                std::abs(spc.hp()));
}

template <typename lps, typename hps, bool acheck, typename PowScalar>
inline ScalarPrecisionCheck<lps,hps,acheck> pow(
        const ScalarPrecisionCheck<lps,hps,acheck>& spc, PowScalar p){
    return ScalarPrecisionCheck<lps,hps,acheck>(std::pow(spc.lp(), p),
                                                std::pow(spc.hp(), p));
}

template <typename lps, typename hps, bool acheck>
inline ScalarPrecisionCheck<lps,hps,acheck> sqrt(
        const ScalarPrecisionCheck<lps,hps,acheck>& spc){
    return ScalarPrecisionCheck<lps,hps,acheck>(std::sqrt(spc.lp()),
                                                std::sqrt(spc.hp()));
}

template <typename lps, typename hps, bool acheck>
inline ScalarPrecisionCheck<lps,hps,acheck> cos(
        const ScalarPrecisionCheck<lps,hps,acheck>& spc){
    return ScalarPrecisionCheck<lps,hps,acheck>(std::cos(spc.lp()),
                                                std::cos(spc.hp()));
}

template <typename lps, typename hps, bool acheck>
inline ScalarPrecisionCheck<lps,hps,acheck> sin(
        const ScalarPrecisionCheck<lps,hps,acheck>& spc){
    return ScalarPrecisionCheck<lps,hps,acheck>(std::sin(spc.lp()),
                                                std::sin(spc.hp()));
}

template <typename lps, typename hps, bool acheck>
inline ScalarPrecisionCheck<lps,hps,acheck> tan(
        const ScalarPrecisionCheck<lps,hps,acheck>& spc){
    return ScalarPrecisionCheck<lps,hps,acheck>(std::tan(spc.lp()),
                                                std::tan(spc.hp()));
}

template <typename lps, typename hps, bool acheck>
inline ScalarPrecisionCheck<lps,hps,acheck> atan(
        const ScalarPrecisionCheck<lps,hps,acheck>& spc){
    return ScalarPrecisionCheck<lps,hps,acheck>(std::atan(spc.lp()),
                                                std::atan(spc.hp()));
}

} // namespace std

#endif // SCALAR_PRECISION_CHECK_H
