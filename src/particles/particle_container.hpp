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
#include "common/box.hpp"
#include "particles/mobile.hpp"
#include "particles/frame.hpp"
#include "particles/osmotic.hpp"
#include "geometries/sphere_grid.hpp"
#include "geometries/plane.hpp"
#include "io/parameters.hpp"
#include "enhance/container_class_base.hpp"
#include <deque>
#include <vector>
#include <algorithm>
#include <memory>
#include <type_traits>
#include <filesystem>
#include <tbb/mutex.h>
#include <boost/multi_array.hpp>

#ifdef __clang_major__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wexceptions"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wnon-virtual-dtor"
#include "h5xx/h5xx.hpp"
#pragma clang diagnostic pop
#elif  __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wexceptions"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#include "h5xx/h5xx.hpp"
#pragma GCC diagnostic pop
#else
    #error no valid compiler
#endif



namespace ves { struct ParticleContainer; };



struct ves::ParticleContainer
    : public enhance::ContainerBase<std::deque<std::unique_ptr<Particle::Base>>>
{
    using particle_ptr_t = element_t;
    using particle_t = element_t::element_type;
    using PATH = Parameters::PATH;
    using FSTREAM = std::ofstream;
    using array1d_t = boost::multi_array<std::uint32_t,1>;
    using array2d_t = boost::multi_array<REAL,2>;
    using cartesian = Particle::Base::cartesian;


    template<typename P, typename ENABLER = typename std::enable_if<std::is_base_of<Particle::Base,P>::value>::type>
    auto addParticle();
    void removeParticle(const particle_t&);
    void removeParticle(const particle_ptr_t&);

    bool placement_conflict(const Particle::Base&, REAL) const;
    void setup();

    template<ves::Particle::TYPE T>
    std::size_t numType() const;

    ves::Box<PERIODIC::ON> box;

    std::size_t getTime() const;

protected:
    void setupFromNew();
    void setupFromH5();

    h5xx::file h5file;
    std::vector<std::string> getGroupNames();
    std::size_t time_from_input = 0;

    tbb::mutex mutex;
};



template<typename P, typename ENABLER>
auto ves::ParticleContainer::addParticle()
{
    tbb::mutex::scoped_lock lock(mutex);
    data.emplace_back(std::make_unique<P>());
}



template<ves::Particle::TYPE T>
std::size_t ves::ParticleContainer::numType() const
{
    return std::accumulate(begin(), end(), std::size_t(0), [](auto i, const particle_ptr_t& p){ return p->getType() == T ? i+1 : i; });
}