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



#include "cell.hpp"



ves::Cell::Cell()
{
    box.setLengthX(Parameters::getInstance().getOption("system.box.x").as<REAL>());
    box.setLengthY(Parameters::getInstance().getOption("system.box.y").as<REAL>());
    box.setLengthZ(Parameters::getInstance().getOption("system.box.z").as<REAL>());
}




bool ves::Cell::areNeighbours(const Cell& a, const Cell& b)
{
    static REAL x = Parameters::getInstance().getOption("system.box.x").as<REAL>();
    static REAL y = Parameters::getInstance().getOption("system.box.y").as<REAL>();
    static REAL z = Parameters::getInstance().getOption("system.box.z").as<REAL>();
    ves::Box<PERIODIC::ON> periodic;
    periodic.setLengthX(x);
    periodic.setLengthY(y);
    periodic.setLengthZ(z);

    const cartesian connection_vector = periodic.distanceVector(a.getBoundaries().center(), b.getBoundaries().center()).cwiseAbs();

    if( a == b )
    {
        // vesDEBUG(__func__ << " " << connection_vector.format(ROWFORMAT) << " return:  false  (Reason: same Cell)")
        return false;
    }
    else if((connection_vector(0) < a.getBoundaries().sizes()(0) + 0.01) 
        &&  (connection_vector(1) < a.getBoundaries().sizes()(1) + 0.01) 
        &&  (connection_vector(2) < a.getBoundaries().sizes()(2) + 0.01) )
    { 
        // vesDEBUG(__func__ << " " << connection_vector.format(ROWFORMAT) << " return:  true  (Reason: in reach)")
        return true;
    }
    else 
    { 
        // vesDEBUG(__func__ << " " << connection_vector.format(ROWFORMAT) << " return:  false (Reason: not in reach)")
        return false;
    }
}



void ves::Cell::setBoundaries(const cartesian& from, const cartesian& to)
{
    bounding_box.setEmpty();
    bounding_box.extend(from);
    bounding_box.extend(to);
    vesDEBUG("cell boundaries were set from " << from.format(ROWFORMAT) << " to " << to.format(ROWFORMAT));
}



auto ves::Cell::getBoundaries() const -> const decltype(bounding_box)&
{
    return bounding_box;
}



void ves::Cell::clear()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    tbb::spin_rw_mutex::scoped_lock lock(particles_access_mutex, true);
    data.clear();
}



void ves::Cell::removeParticle(const particle_t& to_remove)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    tbb::spin_rw_mutex::scoped_lock lock(particles_access_mutex, true);
    data.erase( std::remove_if(begin(), end(), [&](const particle_ptr_t& to_compare)
    { 
        return to_remove == *to_compare;
    ;}), end() );

}



bool ves::Cell::try_add(particle_t* particle)
{
    if(!contains(*particle) && contains(box.scaleDown(particle->getCoordinates())))
    {   
        tbb::spin_rw_mutex::scoped_lock lock(particles_access_mutex, true);
        data.emplace_back(nonstd::make_observer<particle_t>(particle));
        lock.release();
        // vesDEBUG(__func__ << "  Cell contains particle after insertion " << std::boolalpha << contains(&particle))
        assert(contains(*particle));
        return true;
    }
    else
    {
        // vesDEBUG(__func__ << "  Cell contains particle after NO insertion " << std::boolalpha << contains(&particle))
        assert(!contains(*particle));
        return false;
    }
}



bool ves::Cell::contains(const cartesian& c) const
{
    return bounding_box.contains(box.scaleDown(c));
}



bool ves::Cell::contains(const particle_t& other)
{
    tbb::spin_rw_mutex::scoped_lock lock(particles_access_mutex, false);
    return std::any_of(cbegin(), cend(), [&](const auto& p ){ return other == *p; });
}



auto ves::Cell::getLeavers() -> decltype(data)
{
    tbb::spin_rw_mutex::scoped_lock lock(particles_access_mutex, false);
    decltype(data) leavers;
    for(const auto& particle : data)
    {
        if(!contains(particle->getCoordinates()))
        {
            leavers.emplace_back(nonstd::make_observer<particle_t>(particle.get()));
            assert(!contains(leavers.back()->getCoordinates()));
        }
    }
    return leavers;
}



auto ves::Cell::getListFromParticlesInRegion() const
{
    std::deque<particle_ptr_t> particles;
    for(const Cell& cell : region)
        for(const particle_ptr_t& particle : cell.data)
            particles.emplace_back(particle);
    return particles;
}



REAL ves::Cell::potentialOfSingleParticle(const particle_t& particle) const 
{
    return std::accumulate(std::begin(region), std::end(region), REAL(0), [&](REAL __val, const auto& cell)
    {
        return __val + std::accumulate(std::begin(cell.get()), std::end(cell.get()), REAL(0), [&](REAL _val, const auto& compare)
        {
            return particle == *compare ? _val : _val + interaction.calculate(particle, *compare);
        });
    });
    // const auto particles = getListFromParticlesInRegion();
    // REAL sum = 0;
    // for(std::size_t i = 0; i < particles.size(); ++i)
    // {
    //     if(particle != *particles[i] )
    //         sum += interaction.calculate(*particles[i], particle);
    // }
    // return sum;
}



REAL ves::Cell::potentialIgnoreParticle(particle_t& particle) const 
{
    auto particles = getListFromParticlesInRegion();
    const auto size_before = particles.size();
    particles.erase(std::remove_if(std::begin(particles), std::end(particles), [&](const particle_ptr_t& comp){ return particle == *comp; }), std::end(particles));
    if( size_before == particles.size())
        vesCRITICAL("could not erase particle");

    REAL sum = 0;
    for(std::size_t i = 0; i < particles.size(); ++i)
    {
        for(std::size_t j = 0; j<i; ++j)
        {
            // if(*particles[i] != particle && *particles[j] != particle)
                sum += interaction.calculate(*particles[i], *particles[j]);
        }
    }
    return sum;
}



//The chemical potential from computer simulation Test particle method with umbrella sampling 
REAL ves::Cell::chemicalPotentialIgnoreParticle(particle_t& particle) const 
{
    auto particles = getListFromParticlesInRegion();
    const auto size_before = particles.size();
    particles.erase(std::remove_if(std::begin(particles), std::end(particles), [&](const particle_ptr_t& comp){ return particle == *comp; }), std::end(particles));
    if( size_before == particles.size())
        vesCRITICAL("could not erase particle");
    
    REAL sum = 0;
    for(std::size_t i = 0; i < particles.size(); ++i)
    {
        for(std::size_t j = 0; j<i; ++j)
        {
            // if(*particles[i] != particle && *particles[j] != particle)
                sum += interaction.calculate(*particles[i], *particles[j]);
        }
    }
    
    // return complete_sum/(particles.size() - 1);
    return 0;



    // return std::exp(-std::accumulate(std::begin(particles), std::end(particles), REAL(0), [&](REAL _val, const auto& compare)
    // {
    //     return particle == *compare ? _val : _val + interaction.calculate(particle, *compare);
    // })/temperature) / (particles.size() - 1);

    
}



//The chemical potential from computer simulation Test particle method with umbrella sampling 
REAL ves::Cell::chemicalPotentialWithPhantomParticle(particle_t& particle) const 
{
    auto particles = getListFromParticlesInRegion();
    if(contains(particle.getCoordinates()))
    {
        particles.emplace_back(nonstd::make_observer<particle_t>(&particle));
    }

    REAL complete_sum = 0;
    for(const auto& origin : particles)
    {
        complete_sum += std::exp( -std::accumulate(std::begin(particles), std::end(particles), REAL(0), [&](REAL _val, const auto& compare)
        {
            return origin == compare ? _val : _val + interaction.calculate(*origin, *compare);
        }) / temperature);
    }
    
    return complete_sum/(particles.size() - 1);



    // return std::exp(-std::accumulate(std::begin(particles), std::end(particles), REAL(0), [&](REAL _val, const auto& compare)
    // {
    //     return particle == *compare ? _val : _val + interaction.calculate(particle, *compare);
    // })/temperature) / (particles.size() - 1);

    
}



REAL ves::Cell::potentialWithPhantomParticle(particle_t& particle) const 
{
    auto particles = getListFromParticlesInRegion();
    particles.emplace_back(nonstd::make_observer<particle_t>(&particle));
    REAL sum = 0;
    for(std::size_t i = 0; i < particles.size(); ++i)
    {
        for(std::size_t j = 0; j<i; ++j)
        {
            sum += interaction.calculate(*particles[i], *particles[j]);
        }
    }
    return sum;
}



REAL ves::Cell::potential() const 
{
    const auto particles = getListFromParticlesInRegion();
    REAL sum = 0;
    for(std::size_t i = 0; i < particles.size(); ++i)
    {
        for(std::size_t j = 0; j<i; ++j)
        {
            sum += interaction.calculate(*particles[i], *particles[j]);
        }
    }
    return sum;
}