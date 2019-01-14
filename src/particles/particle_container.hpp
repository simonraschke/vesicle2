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
#include "io/parameters.hpp"
#include "enhance/container_class_base.hpp"
#include <deque>
#include <memory>
#include <tbb/mutex.h>
// #include <tbb/cache_aligned_allocator.h>




namespace ves { struct ParticleContainer; };



struct ves::ParticleContainer
    : public enhance::ContainerBase<std::deque<std::unique_ptr<Particle::Base>>>
{
    using particle_ptr_t = element_t;
    using particle_t = element_t::element_type;

    template<typename P, typename ENABLER = typename std::enable_if<std::is_base_of<Particle::Base,P>::value>::type>
    auto addParticle();
    void removeParticle(const particle_t&);
    void removeParticle(const particle_ptr_t&);

    bool placement_conflict(const Particle::Base&, REAL) const;
    void setup();

    template<ves::Particle::TYPE T>
    std::size_t numType() const;

    ves::Box<PERIODIC::ON> box;

protected:
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