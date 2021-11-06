//
// Created by yunfan on 2021/11/6.
//

#ifndef POLY_DATA_STRUCTURE_COMMON_INCLUDE_HPP
#define POLY_DATA_STRUCTURE_COMMON_INCLUDE_HPP

#include "root_finder.hpp"
#include <vector>
#include <list>
#include <Eigen/Eigen>
#include <iostream>
#include "fmt/color.h"
#include "scope_timer.hpp"

typedef Eigen::Matrix<double, 3, 1> Vec3;
typedef Eigen::Matrix<double, 3, 2> StatePV;
typedef Eigen::Matrix<double, 3, 3> StatePVA;
typedef Eigen::Matrix<double, 3, 4> StatePVAJ;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> DynamicMat;

#endif //POLY_DATA_STRUCTURE_COMMON_INCLUDE_HPP
