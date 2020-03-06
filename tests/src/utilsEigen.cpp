/*
 Copyright (C) 2020 Nicolas Mellado <nmellado0@gmail.com>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
 \file tests/src/utilsEigen.cpp
 \brief Test validity of the utility functions for Eigen
 \authors Nicolas Mellado
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <Ponca/Core/utilsEigen.h>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <vector>

int main(int argc, char **argv)
{
  if(!init_testing(argc, argv))
  {
      return EXIT_FAILURE;
  }

  Eigen::VectorXd vec(4);
  vec << 1, 2, 4, 8;
  Eigen::Matrix4d mat = Eigen::Matrix4d::Zero();

  double epsilon = testEpsilon<double>();


  int n = 1000000;
  std::vector<Eigen::Vector4d> vvs (n);
  std::generate(vvs.begin(), vvs.end(), []() { return Eigen::Vector4d::Random();});

  auto start = std::chrono::system_clock::now();
  for(const auto&v : vvs)
    mat += Ponca::covarianceAccumulate(v);
  Eigen::MatrixXd tmp = mat.template selfadjointView<Eigen::Lower>();
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout<<"Sparse total time : "<< elapsed_seconds.count() << std::endl;
  std::cout << tmp << std::endl;

  mat = Eigen::Matrix4d::Zero();

  start = std::chrono::system_clock::now();
  for(const auto&v : vvs)
      mat += v * v.transpose();
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  std::cout<<"Full total time : "<< elapsed_seconds.count() << std::endl;
  std::cout << mat << std::endl;


  VERIFY( (mat.array() - tmp.array()).sum() < epsilon );


  exit(0);
}
