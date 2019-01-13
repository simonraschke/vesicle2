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

#if __has_include(<Eigen/Core>)
#include <Eigen/Geometry>
#elif __has_include(<eigen3/Eigen/Core>)
#include <eigen3/Eigen/Geometry>
#endif

#include <memory>
#include "particles/particle.hpp"
#include "vesicleIO/parameters.hpp"



enum class PERIODIC : bool { ON=true, OFF=false};



// a simulation box implementation
// represents a box from 0 to x,x,z
// is ParameterDependentComponent
//
// may eiter be PERIODIC::ON to calculate with periodic boundary conditions
// or may be PERIODIC::OFF to caclulate without periodic boundary conditions
template<PERIODIC P>
class Box
    : virtual public ParameterDependentComponent
{
public:
    typedef Particle::real real;
    typedef Particle::cartesian cartesian;

    // set x by mutableAccess to Parameter base class
    void setLengthX(real);

    // set y by mutableAccess to Parameter base class
    void setLengthY(real);

    // set z by mutableAccess to Parameter base class
    void setLengthZ(real);

    real getLengthX() const;
    real getLengthY() const;
    real getLengthZ() const;

    cartesian getCenter() const;

    // calcaulates the distance vector of two particles
    // depending on PERIODIC ON or OFF
    // called from anywhere else
    cartesian distanceVector(const cartesian&, const cartesian&) const;
    // implementation for Particle base class. calls cartesian version
    cartesian distanceVector(const Particle&, const Particle&) const;


    // distance
    // calls squared_distance and calculates std::sqrt
    real distance(const cartesian&, const cartesian&) const;
    // implementation for Particle base class. calls cartesian version
    real distance(const Particle&, const Particle&) const;


    // squared distance
    real squared_distance(const cartesian&, const cartesian&) const;
    // implementation for Particle base class. calls cartesian version
    real squared_distance(const Particle&, const Particle&) const;


    // scales down any give coordinates into the simulation box
    cartesian scaleDown(cartesian) const ;
    // implementation for Particle base class. calls cartesian version
    cartesian scaleDown(const Particle&) const;


    // scales down any give coordinates into the simulation box
    // VMD simulation box is from -x/2 to x/2
    // scales accordingly
    cartesian scaleDownForVMD(cartesian) const ;
    // implementation for Particle base class. calls cartesian version
    cartesian scaleDownForVMD(const Particle&) const;

    // checks if bounding_box contains coordinates
    // if PERIODIC::ON calls scaleDown before
    bool contains(const cartesian&) const;
    // implementation for Particle base class. calls cartesian version
    bool contains(const Particle&) const;


    // destroy if derived is destroyed
    virtual ~Box() = default;

    // check if all parameters are set to make bounding_box
    // necessary for contains(const Particle&)
    void check_for_aligned_box_setup();

protected:
    using ParameterDependentComponent::mutableAccess;

private:

    std::unique_ptr<Eigen::AlignedBox<real,3>> bounding_box {nullptr};
};



template<PERIODIC P>
void Box<P>::setLengthX(real l)
{
    mutableAccess().x = l;
    check_for_aligned_box_setup();
}



template<PERIODIC P>
void Box<P>::setLengthY(real l)
{
    mutableAccess().y = l;
    check_for_aligned_box_setup();
}



template<PERIODIC P>
void Box<P>::setLengthZ(real l)
{
    mutableAccess().z = l;
    check_for_aligned_box_setup();
}



template<PERIODIC P>
typename Box<P>::real Box<P>::getLengthX() const
{
    return getParameters().x;
}



template<PERIODIC P>
typename Box<P>::real Box<P>::getLengthY() const
{
    return getParameters().x;
}



template<PERIODIC P>
typename Box<P>::real Box<P>::getLengthZ() const
{
    return getParameters().z;
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
    cartesian vec;
    vec(0) = getParameters().x;
    vec(1) = getParameters().y;
    vec(2) = getParameters().z;
    bounding_box = std::make_unique<Eigen::AlignedBox<real,3>>( cartesian::Zero(), vec );
}



template<>
EIGEN_STRONG_INLINE Box<PERIODIC::ON>::cartesian Box<PERIODIC::ON>::distanceVector(const cartesian& c1, const cartesian& c2) const
{
    cartesian distance_cartesian = c2-c1;
    distance_cartesian(0) = distance_cartesian(0) - getParameters().x * std::round(distance_cartesian(0)/(getParameters().x));
    distance_cartesian(1) = distance_cartesian(1) - getParameters().y * std::round(distance_cartesian(1)/(getParameters().y));
    distance_cartesian(2) = distance_cartesian(2) - getParameters().z * std::round(distance_cartesian(2)/(getParameters().z));
    return distance_cartesian;
}



template<>
EIGEN_STRONG_INLINE Box<PERIODIC::OFF>::cartesian Box<PERIODIC::OFF>::distanceVector(const cartesian& c1, const cartesian& c2) const
{
    return (c2-c1);
}



template<PERIODIC P>
EIGEN_STRONG_INLINE typename Box<P>::cartesian Box<P>::distanceVector(const Particle&p1, const Particle& p2) const
{
    return distanceVector(p1.coords(),p2.coords());
}



template<PERIODIC P>
EIGEN_STRONG_INLINE typename Box<P>::real Box<P>::squared_distance(const cartesian& c1, const cartesian& c2) const 
{
    return distanceVector(c1,c2).squaredNorm();
}



template<PERIODIC P>
EIGEN_STRONG_INLINE typename Box<P>::real Box<P>::squared_distance(const Particle& p1, const Particle& p2) const 
{
    return squared_distance(p1.coords(),p2.coords());
}



template<PERIODIC P>
EIGEN_STRONG_INLINE typename Box<P>::real Box<P>::distance(const cartesian& c1, const cartesian& c2) const 
{
    return std::sqrt(squared_distance(c1,c2));
}



template<PERIODIC P>
EIGEN_STRONG_INLINE typename Box<P>::real Box<P>::distance(const Particle& p1, const Particle& p2) const 
{
    return distance(p1.coords(),p2.coords());
}



template<PERIODIC P>
typename Box<P>::cartesian Box<P>::scaleDown(cartesian c) const 
{
    while( c(0) > getParameters().x ) c(0) -= getParameters().x;
    while( c(1) > getParameters().y ) c(1) -= getParameters().y;
    while( c(2) > getParameters().z ) c(2) -= getParameters().z;

    while( c(0) < 0.f ) c(0) += getParameters().x;
    while( c(1) < 0.f ) c(1) += getParameters().y;
    while( c(2) < 0.f ) c(2) += getParameters().z;
    return c;
}



template<PERIODIC P>
typename Box<P>::cartesian Box<P>::scaleDown(const Particle& p) const 
{
    return scaleDown(p.coords());
}



template<PERIODIC P>
typename Box<P>::cartesian Box<P>::scaleDownForVMD(cartesian c) const 
{
    c(0) = c(0) - getParameters().x * std::round(c(0)/(getParameters().x));
    c(1) = c(1) - getParameters().y * std::round(c(1)/(getParameters().y));
    c(2) = c(2) - getParameters().z * std::round(c(2)/(getParameters().z));
    return c;
}



template<PERIODIC P>
typename Box<P>::cartesian Box<P>::scaleDownForVMD(const Particle& p) const 
{
    return scaleDownForVMD(p.coords());
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
inline bool Box<P>::contains(const Particle& p) const 
{
    return contains(p.coords());
}