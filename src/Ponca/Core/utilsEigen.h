/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
 \author Nicolas Mellado
*/

#pragma once

#include <Eigen/Dense>

namespace Ponca
{

  template<class ArgType>
  struct _symmetric_cov_helper {
      using MatrixType = Eigen::Matrix<typename ArgType::Scalar,
      ArgType::SizeAtCompileTime,
      ArgType::SizeAtCompileTime,
      Eigen::ColMajor,
      ArgType::MaxSizeAtCompileTime,
      ArgType::MaxSizeAtCompileTime>;
  };

  template<class ArgType>
  class _symmetric_cov_functor {
    using Scalar = typename ArgType::Scalar;
  public:
      inline _symmetric_cov_functor(Scalar w, const ArgType& arg) : m_w{w}, m_vec{arg} {}
    inline Scalar operator() (Eigen::Index row, Eigen::Index col) const {
      return (col <= row) ? m_w * m_vec(col) * m_vec(row) : Scalar(0);
    }
  private:
    Scalar m_w;
    const ArgType &m_vec;
  };

  template <class ArgType>
  Eigen::CwiseNullaryOp<_symmetric_cov_functor<ArgType>,
                        typename _symmetric_cov_helper<ArgType>::MatrixType>
  covarianceAccumulate(const Eigen::MatrixBase<ArgType>& arg,
                       typename ArgType::Scalar w = typename ArgType::Scalar(1.))
  {
    using MatrixType = typename _symmetric_cov_helper<ArgType>::MatrixType;
    return MatrixType::NullaryExpr(arg.size(), arg.size(),
                                   _symmetric_cov_functor<ArgType>(w, arg.derived()));
  }

}
