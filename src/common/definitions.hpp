/*  
*   Copyright 2019 Simon Raschke
*
*   Licensed under the Apache License, Version 2.0 (the "License");
*   you may not use this file except in compliance with the License.
*   You may obtain a copy of the License at
*
*       http://www.apache.org/licenses/LICENSE-2.0
*
*   Unless required by applicable law or agreed to in writing, software
*   distributed under the License is distributed on an "AS IS" BASIS,
*   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*   See the License for the specific language governing permissions and
*   limitations under the License.
*/

#pragma once

#define REAL float

constexpr const double PI_4  = 0.785398163397448309616;
constexpr const double PI_2  = 1.57079632679489661923;
constexpr const double PI    = 3.14159265358979323846;
constexpr const double TWOPI = 3.14159265358979323846*2;

#define likely(x)      __builtin_expect(!!(x), 1)
#define unlikely(x)    __builtin_expect(!!(x), 0)


#include <csignal>
#include <iostream>
#ifndef NDEBUG
    #define vesDEBUG(x) {std::cerr << "[DEBUG] "; do { std::cerr << x; } while (0); std::cerr << '\n';}
#else
    #define vesDEBUG(x)
#endif
#define vesLOG(x) {std::clog << "[LOG] "; do { std::clog << x; } while (0); std::clog << '\n';}
#define vesWARNING(x) {std::clog << "[WARNING] "; do { std::clog << x; } while (0); std::clog << '\n';}
#define vesCRITICAL(x) {std::cerr << "[ERROR] "<< __FILE__ <<":" << __LINE__ << "  "; do { std::cerr << x; } while (0); std::cerr <<" raising SIGABRT\n"; std::exit(SIGABRT);}


#if __has_include(<Eigen/Core>)
#include <Eigen/Core>
#elif __has_include(<eigen3/Eigen/Core>)
#include <eigen3/Eigen/Core>
#endif
//-------------------Eigen::IOFormat( prec, flag,                 coeffSep, rowSep, rowPre, rowSuf, matPre, matSuf )
#define ROWFORMAT    Eigen::IOFormat( 3,    Eigen::DontAlignCols, ", ",     " ",    " ",    "",     " ",    " " )
#define PYTHONFORMAT Eigen::IOFormat( 4,    0,                    ", ",     "\n",   "[",    "]",    "[",    "]" )