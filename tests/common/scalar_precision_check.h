/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.

 \author Nicolas Mellado nmellado0@gmail.com
 \author Gael Guennebaud gael.guennebaud@inria.com
*/


#ifndef SCALAR_PRECISION_CHECK_H
#define SCALAR_PRECISION_CHECK_H

#include <limits>
#include <cassert>
#include <cmath>
#include <algorithm>

template <typename _lowPrecisionScalarT,typename _highPrecisionScalarT>
class ScalarPrecisionCheck;

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
    ScalarPC::check_enabled = true;

    ScalarPC a = 2;
    ScalarPC b = 10e-6;
    ScalarPC c = pow(a, 3) * b;

    std::cout << c.lp() << " " << c.hp() << std::endl;
    \endcode
 */
template <typename _lowPrecisionScalarT  = float,
          typename _highPrecisionScalarT = double>
class ScalarPrecisionCheck {
public:
    typedef _lowPrecisionScalarT  LPScalar;
    typedef _highPrecisionScalarT HPScalar;

    inline ScalarPrecisionCheck() {}
    
    /*!
     */
    template <typename RealScalar>
    inline ScalarPrecisionCheck(RealScalar v)
        : _lp(v), _hp(v) {
        precision_assert(*this);
    }

    inline ScalarPrecisionCheck(const ScalarPrecisionCheck& __rhs)
        : _lp(__rhs._lp), _hp(__rhs._hp) {
        precision_assert(*this);
        
    }
    
    inline ScalarPrecisionCheck(LPScalar lp,
                                HPScalar hp)
        : _lp(lp), _hp(hp) {
        //if(CheckEnabled && relativePrecisionDelta()>0)
        //  std::cerr << "    " << _lp << " " << _hp << " -> " << relativePrecisionDelta() << "\n";
        precision_assert(*this);
    }

    inline HPScalar precisionDelta() const { return std::abs(HPScalar(_lp) - _hp); }
    inline HPScalar relativePrecisionDelta() const {
      using std::abs;
      using std::max;
      using std::min;
      HPScalar diff = abs(HPScalar(_lp) - _hp);
      if(min(abs(HPScalar(_lp)),abs(_hp))<7e-5)
        return 0;
      return diff<std::numeric_limits<LPScalar>::min() ? 0 : diff/abs(HPScalar(_hp));
    }
    inline bool checkPrecision() const {
      return ( relativePrecisionDelta() < epsilon);
    }

    operator const LPScalar () const { return _lp; }
    inline LPScalar lp() const { return _lp; }
    inline HPScalar hp() const { return _hp; }

    ////////////////////////////////////////////////////////////////////////////
    /// Pairwise operators between ScalarPrecisionChecks
    ///
    inline ScalarPrecisionCheck operator*(const ScalarPrecisionCheck& __rhs) const
    { return ScalarPrecisionCheck ( _lp * __rhs._lp, _hp * __rhs._hp ); }

    inline ScalarPrecisionCheck operator/(const ScalarPrecisionCheck& __rhs) const
    { return ScalarPrecisionCheck ( _lp / __rhs._lp, _hp / __rhs._hp ); }

    inline ScalarPrecisionCheck operator+(const ScalarPrecisionCheck& __rhs) const
    { return ScalarPrecisionCheck ( _lp + __rhs._lp, _hp + __rhs._hp ); }

    inline ScalarPrecisionCheck operator-(const ScalarPrecisionCheck& __rhs) const
    { return ScalarPrecisionCheck ( _lp - __rhs._lp, _hp - __rhs._hp ); }

    inline ScalarPrecisionCheck& operator*=(const ScalarPrecisionCheck& __rhs)
    { return *this = __rhs * (*this); }

    inline ScalarPrecisionCheck& operator/=(const ScalarPrecisionCheck& __rhs)
    { return *this = __rhs / (*this); }

     inline ScalarPrecisionCheck& operator+=(const ScalarPrecisionCheck& __rhs)
    { return *this = __rhs + (*this); }

     inline ScalarPrecisionCheck& operator-=(const ScalarPrecisionCheck& __rhs)
    { return *this = __rhs - (*this); }

    inline bool operator<(const ScalarPrecisionCheck& __rhs) const
    { return AND_assert( _lp < __rhs._lp, _hp < __rhs._hp ); }


    inline bool operator>(const ScalarPrecisionCheck& __rhs) const
    { return AND_assert( _lp > __rhs._lp, _hp > __rhs._hp ); }

    inline bool operator==(const ScalarPrecisionCheck& __rhs) const
    { return AND_assert( _lp == __rhs._lp, _hp == __rhs._hp ); }


    inline bool operator!=(const ScalarPrecisionCheck& __rhs) const
    { return AND_assert( _lp != __rhs._lp, _hp != __rhs._hp ); }

public:
    static HPScalar epsilon;
    static bool check_enabled;

private:
    LPScalar _lp;
    HPScalar _hp;
    static inline const ScalarPrecisionCheck& precision_assert(const ScalarPrecisionCheck& s)
    {
      if(check_enabled) assert(s.checkPrecision());
      return s;
    }
    
    static inline bool AND_assert(bool lp, bool hp)
    {
      if(check_enabled) assert (lp==hp);
      return lp && hp;
    }
    
}; // class ScalarPrecisionCheck

template <typename lps, typename hps>
hps ScalarPrecisionCheck<lps,hps>::epsilon = 2*std::numeric_limits<lps>::epsilon();

template <typename lps, typename hps>
bool ScalarPrecisionCheck<lps,hps>::check_enabled = true;

// create definitions for common std math functions
namespace std{
template <typename lps, typename hps>
inline ScalarPrecisionCheck<lps,hps> abs(const ScalarPrecisionCheck<lps,hps>& spc){
    return ScalarPrecisionCheck<lps,hps>(std::abs(spc.lp()), std::abs(spc.hp()));
}

template <typename lps, typename hps, typename PowScalar>
inline ScalarPrecisionCheck<lps,hps> pow(const ScalarPrecisionCheck<lps,hps>& spc, PowScalar p){
    return ScalarPrecisionCheck<lps,hps>(std::pow(spc.lp(), p), std::pow(spc.hp(), p));
}

template <typename lps, typename hps>
inline ScalarPrecisionCheck<lps,hps> sqrt(const ScalarPrecisionCheck<lps,hps>& spc){
    return ScalarPrecisionCheck<lps,hps>(std::sqrt(spc.lp()), std::sqrt(spc.hp()));
}

template <typename lps, typename hps>
inline ScalarPrecisionCheck<lps,hps> cos(const ScalarPrecisionCheck<lps,hps>& spc){
    return ScalarPrecisionCheck<lps,hps>(std::cos(spc.lp()), std::cos(spc.hp()));
}

template <typename lps, typename hps>
inline ScalarPrecisionCheck<lps,hps> sin(
        const ScalarPrecisionCheck<lps,hps>& spc){
    return ScalarPrecisionCheck<lps,hps>(std::sin(spc.lp()),
                                                std::sin(spc.hp()));
}

template <typename lps, typename hps>
inline ScalarPrecisionCheck<lps,hps> tan(const ScalarPrecisionCheck<lps,hps>& spc){
    return ScalarPrecisionCheck<lps,hps>(std::tan(spc.lp()), std::tan(spc.hp()));
}

template <typename lps, typename hps>
inline ScalarPrecisionCheck<lps,hps> atan(const ScalarPrecisionCheck<lps,hps>& spc){
    return ScalarPrecisionCheck<lps,hps>(std::atan(spc.lp()), std::atan(spc.hp()));
}

} // namespace std

#endif // SCALAR_PRECISION_CHECK_H
