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

#include "common/definitions.hpp"
#include "common/component.hpp"

// #include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <cmath>
#include <memory>
#include <exception>



namespace ves 
{
    struct Bounding; 
    struct CoordinatesBounding; 
    struct OrientationBounding; 
}



struct ves::Bounding
{
protected:
    Bounding() = default;
    std::unique_ptr<Eigen::Matrix<REAL,3,1>> origin {nullptr};

public:
    virtual ~Bounding() = default;
    virtual bool isAllowed(const Eigen::Matrix<REAL,3,1>&) = 0;
};



struct ves::OrientationBounding
    : public Bounding
{
    float value = TWOPI;

    inline virtual bool isAllowed(const Eigen::Matrix<REAL,3,1>& compare) override
    {
        return std::acos(compare.normalized().dot(origin->normalized())) < value; 
    };
};



struct ves::CoordinatesBounding
    : public Bounding
{
    inline virtual bool isAllowed(const Eigen::Matrix<REAL,3,1> & compare) override
    {
        if(!origin)
        {
            origin = std::make_unique<Eigen::Matrix<REAL,3,1>>(compare);
        }

        if(bounding_box)
        {
            return bounding_box->contains(compare);
        }
        else if(sphere_bounds)
        {
            const auto offset = (compare - (*origin)).norm();
            return offset >= sphere_bounds->first && offset <= sphere_bounds->second;
        }
        else if(!bounding_box && !sphere_bounds)
        {
            return true;
        }
        else
        {
            vesCRITICAL(__PRETTY_FUNCTION__ << " bad decision");
            return false;
        }
    };



    inline void setBoundingSphere(const REAL from, const REAL to)
    {
        sphere_bounds = std::make_unique<std::pair<REAL,REAL>>(std::make_pair(from,to));
        bounding_box.reset(nullptr);
    };



    inline void setBoundingBox(Eigen::AlignedBox<REAL,3>&& _box)
    {
        bounding_box = std::make_unique<Eigen::AlignedBox<REAL,3>>(std::move(_box));
        sphere_bounds.reset(nullptr);
    };



protected:
    std::unique_ptr<Eigen::AlignedBox<REAL,3>> bounding_box = {nullptr};
    std::unique_ptr<std::pair<REAL,REAL>> sphere_bounds = {nullptr};
};