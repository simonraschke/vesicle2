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

#include "io/parameters.hpp"
#include "particles/base.hpp"
#include "enhance/random.hpp"
#include <memory>
#if __has_include(<Eigen/Core>)
#include <Eigen/Geometry>
#elif __has_include(<eigen3/Eigen/Core>)
#include <eigen3/Eigen/Geometry>
#endif



namespace ves {
    enum class PERIODIC : bool 
    { 
        ON=true, 
        OFF=false}
    ;

    template<PERIODIC P> class Box;
};



// a simulation box implementation
// represents a box from 0 to x,x,z
// is ParameterDependentComponent
//
// may eiter be PERIODIC::ON to calculate with periodic boundary conditions
// or may be PERIODIC::OFF to caclulate without periodic boundary conditions
namespace ves{
template<PERIODIC P>
class Box
{
public:
    typedef Eigen::Matrix<REAL,3,1> cartesian;

    // set x by mutableAccess to Parameter base class
    void setLengthX(REAL);

    // set y by mutableAccess to Parameter base class
    void setLengthY(REAL);

    // set z by mutableAccess to Parameter base class
    void setLengthZ(REAL);

    REAL getLengthX() const;
    REAL getLengthY() const;
    REAL getLengthZ() const;

    cartesian getCenter() const;

    // calcaulates the distance vector of two ves::Particle::Bases
    // depending on PERIODIC ON or OFF
    // called from anywhere else
    cartesian distanceVector(const cartesian&, const cartesian&) const;
    // implementation for ves::Particle::Base base class. calls cartesian version
    cartesian distanceVector(const ves::Particle::Base&, const ves::Particle::Base&) const;


    // distance
    // calls squared_distance and calculates std::sqrt
    REAL distance(const cartesian&, const cartesian&) const;
    // implementation for ves::Particle::Base base class. calls cartesian version
    REAL distance(const ves::Particle::Base&, const ves::Particle::Base&) const;


    // squared distance
    REAL squared_distance(const cartesian&, const cartesian&) const;
    // implementation for ves::Particle::Base base class. calls cartesian version
    REAL squared_distance(const ves::Particle::Base&, const ves::Particle::Base&) const;


    // scales down any give coordinates into the simulation box
    cartesian scaleDown(cartesian) const ;
    // implementation for ves::Particle::Base base class. calls cartesian version
    cartesian scaleDown(const ves::Particle::Base&) const;


    // scales down any give coordinates into the simulation box
    // VMD simulation box is from -x/2 to x/2
    // scales accordingly
    cartesian scaleDownForVMD(cartesian) const ;
    // implementation for ves::Particle::Base base class. calls cartesian version
    cartesian scaleDownForVMD(const ves::Particle::Base&) const;

    // checks if bounding_box contains coordinates
    // if PERIODIC::ON calls scaleDown before
    bool contains(const cartesian&) const;
    // implementation for ves::Particle::Base base class. calls cartesian version
    bool contains(const ves::Particle::Base&) const;


    // destroy if derived is destroyed
    virtual ~Box() = default;

    // check if all parameters are set to make bounding_box
    // necessary for contains(const ves::Particle::Base&)
    void check_for_aligned_box_setup();

    cartesian randomPointInside() const;

protected:

private:
    REAL x {0};
    REAL y {0};
    REAL z {0};

    std::unique_ptr<Eigen::AlignedBox<REAL,3>> bounding_box {nullptr};
};



template<PERIODIC P>
void Box<P>::setLengthX(REAL l)
{
    x = l;
    check_for_aligned_box_setup();
}



template<PERIODIC P>
void Box<P>::setLengthY(REAL l)
{
    y = l;
    check_for_aligned_box_setup();
}



template<PERIODIC P>
void Box<P>::setLengthZ(REAL l)
{
    z = l;
    check_for_aligned_box_setup();
}



template<PERIODIC P>
REAL Box<P>::getLengthX() const
{
    return x;
}



template<PERIODIC P>
REAL Box<P>::getLengthY() const
{
    return x;
}



template<PERIODIC P>
REAL Box<P>::getLengthZ() const
{
    return z;
}



template<PERIODIC P>
typename Box<P>::cartesian Box<P>::getCenter() const
{
    return bounding_box->center();
}



template<PERIODIC P>
void Box<P>::check_for_aligned_box_setup()
{
    bounding_box.reset(nullptr);
    cartesian vec (x,y,z);
    bounding_box = std::make_unique<Eigen::AlignedBox<REAL,3>>( cartesian::Zero(), vec );
}



template<>
EIGEN_STRONG_INLINE Box<PERIODIC::ON>::cartesian Box<PERIODIC::ON>::distanceVector(const cartesian& c1, const cartesian& c2) const
{
    cartesian distance_cartesian = c2-c1;
    distance_cartesian(0) = distance_cartesian(0) - x * std::round(static_cast<REAL>(distance_cartesian(0)/(x)));
    distance_cartesian(1) = distance_cartesian(1) - y * std::round(static_cast<REAL>(distance_cartesian(1)/(y)));
    distance_cartesian(2) = distance_cartesian(2) - z * std::round(static_cast<REAL>(distance_cartesian(2)/(z)));
    return distance_cartesian;
}



template<>
EIGEN_STRONG_INLINE Box<PERIODIC::OFF>::cartesian Box<PERIODIC::OFF>::distanceVector(const cartesian& c1, const cartesian& c2) const
{
    return (c2-c1);
}



template<PERIODIC P>
EIGEN_STRONG_INLINE typename Box<P>::cartesian Box<P>::distanceVector(const ves::Particle::Base& p1, const ves::Particle::Base& p2) const
{
    return distanceVector(p1.getCoordinates(),p2.getCoordinates());
}



template<PERIODIC P>
EIGEN_STRONG_INLINE REAL Box<P>::squared_distance(const cartesian& c1, const cartesian& c2) const 
{
    return distanceVector(c1,c2).squaredNorm();
}



template<PERIODIC P>
EIGEN_STRONG_INLINE REAL Box<P>::squared_distance(const ves::Particle::Base& p1, const ves::Particle::Base& p2) const 
{
    return squared_distance(p1.getCoordinates(),p2.getCoordinates());
}



template<PERIODIC P>
EIGEN_STRONG_INLINE REAL Box<P>::distance(const cartesian& c1, const cartesian& c2) const 
{
    return std::sqrt(squared_distance(c1,c2));
}



template<PERIODIC P>
EIGEN_STRONG_INLINE REAL Box<P>::distance(const ves::Particle::Base& p1, const ves::Particle::Base& p2) const 
{
    return distance(p1.getCoordinates(),p2.getCoordinates());
}



template<PERIODIC P>
typename Box<P>::cartesian Box<P>::scaleDown(cartesian c) const 
{
    while( c(0) > x ) c(0) -= x;
    while( c(1) > y ) c(1) -= y;
    while( c(2) > z ) c(2) -= z;

    while( c(0) < 0.f ) c(0) += x;
    while( c(1) < 0.f ) c(1) += y;
    while( c(2) < 0.f ) c(2) += z;
    return c;
}



template<PERIODIC P>
typename Box<P>::cartesian Box<P>::scaleDown(const ves::Particle::Base& p) const 
{
    return scaleDown(p.getCoordinates());
}



template<PERIODIC P>
typename Box<P>::cartesian Box<P>::scaleDownForVMD(cartesian c) const 
{
    c(0) = c(0) - x * std::round(c(0)/(x));
    c(1) = c(1) - y * std::round(c(1)/(y));
    c(2) = c(2) - z * std::round(c(2)/(z));
    return c;
}



template<PERIODIC P>
typename Box<P>::cartesian Box<P>::scaleDownForVMD(const ves::Particle::Base& p) const 
{
    return scaleDownForVMD(p.getCoordinates());
}



template<>
inline bool Box<PERIODIC::ON>::contains(const cartesian& c) const 
{
    assert(bounding_box);
    return bounding_box->contains(scaleDown(c));
}



template<>
inline bool Box<PERIODIC::OFF>::contains(const cartesian& c) const 
{
    assert(bounding_box);
    return bounding_box->contains(c);
}



template<PERIODIC P>
inline bool Box<P>::contains(const ves::Particle::Base& p) const 
{
    return contains(p.getCoordinates());
}



template<PERIODIC P>
EIGEN_STRONG_INLINE typename Box<P>::cartesian Box<P>::randomPointInside() const 
{
    return cartesian
    (
        enhance::random<cartesian::Scalar>(0.f,getLengthX()),
        enhance::random<cartesian::Scalar>(0.f,getLengthY()),
        enhance::random<cartesian::Scalar>(0.f,getLengthZ())
    );
}

}; // namespace ves