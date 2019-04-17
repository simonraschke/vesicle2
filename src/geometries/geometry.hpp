/*  
*   Copyright 2017-2019 Simon Raschke
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



#include "common/definitions.hpp"
#include <cmath>
#include <vector>
#if __has_include(<Eigen/Core>)
#include <Eigen/Core>
#include <Eigen/Geometry>
#elif __has_include(<eigen3/Eigen/Core>)
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#endif



namespace ves 
{
    struct Geometry;
    struct GridGeometry;
    struct PlaneGeometry;
    struct SphereGeometry;
    struct SphereGridGeometry;
    struct HelixGeometry;
    struct CircleGeometry;
    struct TubeGeometry;
}



struct ves::Geometry
{
    using cartesian = Eigen::Matrix<REAL,3,1>;

    virtual ~Geometry() = default;
    virtual void generate() = 0;
    virtual void scale(const cartesian&) = 0;
    virtual void shift(const cartesian&) = 0;

    std::vector<cartesian> points {};

protected:
    Geometry() = default;
};