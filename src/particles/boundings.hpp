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

public:
    std::unique_ptr<Eigen::Matrix<REAL,3,1>> origin {nullptr};

    virtual ~Bounding() = default;
    virtual bool isAllowed(const decltype(origin)::element_type&) = 0;
    virtual void forcefullyShift(const Eigen::Matrix<REAL,3,1>&) = 0;
};



struct ves::OrientationBounding
    : public Bounding
{
    REAL value = TWOPI;

    inline virtual bool isAllowed(const decltype(origin)::element_type& compare) override
    {
        if(!origin)
        {
            origin = std::make_unique<decltype(origin)::element_type>(compare);
            return true;
        }
        else
        {
            return (std::acos(compare.normalized().dot(origin->normalized())) < value); 
            const bool allowed = (std::acos(compare.normalized().dot(origin->normalized())) < value); 
            if(allowed)
            {
                origin = std::make_unique<decltype(origin)::element_type>(compare);
            }
            return allowed;
        }
    };

    virtual void forcefullyShift(const Eigen::Matrix<REAL,3,1>& s) override
    {
        if(origin)
            origin->operator+=(s);
    }
};



struct ves::CoordinatesBounding
    : public Bounding
{
    inline virtual bool isAllowed(const decltype(origin)::element_type& compare) override
    {
        if(!origin)
        {
            bool allowed = false;
            if(bounding_box)
            {
                allowed = bounding_box->contains(compare);
            }
            else if(sphere_bounds)
            {
                const auto offset = (compare - (*origin)).norm();
                allowed = offset >= sphere_bounds->first && offset <= sphere_bounds->second;
            }
            else if(!bounding_box && !sphere_bounds)
            {
                allowed = true;
            }
            else
            {
                vesCRITICAL(__PRETTY_FUNCTION__ << " bad decision");
                allowed = false;
            }
            
            if(allowed)
                origin = std::make_unique<decltype(origin)::element_type>(compare);
            
            return allowed;
        }
        else
        {
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



    virtual void forcefullyShift(const Eigen::Matrix<REAL,3,1>& s) override
    {
        if(origin)
            origin->operator+=(s);
        if(bounding_box)
            bounding_box->translate(s);
    }



    inline bool isBoxBound() const { return static_cast<bool>(bounding_box); };
    inline bool isSphereBound() const { return static_cast<bool>(sphere_bounds); };

    inline auto getBoundingBox() const { return *bounding_box; };
    inline auto getSphereBounds() const { return *sphere_bounds; };


protected:
    std::unique_ptr<Eigen::AlignedBox<REAL,3>> bounding_box = {nullptr};
    std::unique_ptr<std::pair<REAL,REAL>> sphere_bounds = {nullptr};
};