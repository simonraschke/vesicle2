/*  
*   Copyright 2017-2018 Simon Raschke
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

#include <cmath>
#include <cstdint>
#include <sstream>
#include <exception>
#if __has_include(<Eigen/Core>)
#include <Eigen/Core>
#elif __has_include(<eigen3/Eigen/Core>)
#include <eigen3/Eigen/Core>
#endif
    


namespace enhance
{
    // number of digits
    template <typename T>
    std::uint16_t numDigits(T number)
    {
        std::uint16_t digits = 0;
        if (number < 0) digits = 1; // remove this line if '-' counts as a digit
        std::uint32_t helper = static_cast<int>(number);
        while (helper) 
        {
            helper /= 10;
            digits++;
        }
        return digits;
    }
    


    // angle between two eigen vectors
    // vectors dont have to be normalized
    // will give signed result
    template<typename DERIVED1, typename DERIVED2>
    constexpr float directed_angle(const DERIVED1& v1, const DERIVED2& v2)
    {
        assert(std::isfinite(std::atan2(v2.normalized()(1), v2.normalized()(0)) - std::atan2(v1.normalized()(1), v1.normalized()(0))));
        return std::atan2(v2.normalized()(1), v2.normalized()(0)) - std::atan2(v1.normalized()(1), v1.normalized()(0));
    }
    


    // will give unsigned result
    template<typename DERIVED1, typename DERIVED2>
    constexpr float absolute_angle(const DERIVED1& v1, const DERIVED2& v2)
    {
        return std::abs(directed_angle(v1,v2));
    }
    

    
    //directed angle normalized to 0 to 360°
    template<typename DERIVED1, typename DERIVED2>
    constexpr float normalized_angle(const DERIVED1& v1, const DERIVED2& v2)
    {
        const float angle = directed_angle(v1,v2);
        return angle < 0.f ? angle+M_PI : angle;
    }
        
    
    
    // calculate rad from given deg
    // expects floating pioint type
    template<typename T, typename ENABLER = typename std::enable_if<std::is_floating_point<T>::value>::type>
    constexpr T deg_to_rad(const T& __deg) noexcept
    {
        return __deg*M_PI/180;
    }
    
    
    
    // calculate deg from given rad
    // expects floating pioint type
    template<typename T, typename ENABLER = typename std::enable_if<std::is_floating_point<T>::value>::type>
    constexpr T rad_to_deg(const T& __rad) noexcept
    {
        return __rad/M_PI*180;
    }
    
    
    
    // calculates the volume of a sphere given its radius
    // may throw if radius is negative
    template<typename T>
    constexpr T sphere_volume(const T& __rad) 
    {
        if(__rad < 0) throw std::logic_error("radius must not be negative");
        return M_PI*__rad*__rad*__rad*4.f/3.f;
    }
    
    
    
    // calculates the surface of a sphere given its radius
    // may throw if radius is negative
    template<typename T>
    constexpr T sphere_surface(const T& __rad) 
    {
        if(__rad < 0) throw std::logic_error("radius must not be negative");
        return M_PI*__rad*__rad*4.f;
    }
    
    
    
    // calculates the area of a circle given its radius
    // may throw if radius is negative
    template<typename T>
    constexpr T circle_area(const T& __rad) 
    {
        if(__rad < 0) throw std::logic_error("radius must not be negative");
        return M_PI*__rad*__rad;
    }
    
    
    
    // calculates the radius of a circle given its area
    // may throw if area is negative
    template<typename T>
    constexpr T circle_area_to_radius(const T& __area)
    {
        if(__area < 0) throw std::logic_error("area must not be negative");
        return std::sqrt( __area/M_PI );
    }
    
    
    
    // calculates the volume of a cone given its base radius and heigt
    // may throw if area is negative
    // may throw if height is negative
    template<typename T>
    constexpr T cone_volume(const T& __rad, const T& __height) 
    {
        if(__rad < 0) throw std::logic_error("radius must not be negative");
        if(__height < 0) throw std::logic_error("height must not be negative");
        return enhance::circle_area(__rad)*__height/(static_cast<T>(3));
    }



    template<std::size_t N, typename T>
    constexpr typename std::enable_if<std::is_floating_point<T>::value,T>::type nth_root(const T& __val )
    {
        return std::pow(__val, 1.0/N);
    }
}