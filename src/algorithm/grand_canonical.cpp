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



#include "grand_canonical.hpp"



bool ves::GrandCanonicalInsertion::isValid(const MonteCarloSystem& sys, REAL delta_U) const
{
    static const REAL thermal_wavelength_cubic = std::pow(std::sqrt(REAL(1) / (PI*2*temperature)), 3);

    const REAL prob = sys.getBox().get().getVolume() * std::exp(-(delta_U-mu)/temperature) / thermal_wavelength_cubic / (sys.getParticles().get().data.size() + 1);

    // vesLOG("");
    // vesLOG("INSERTION");
    // vesLOG("delta_U " << delta_U);
    // vesLOG("mu " << mu);
    // vesLOG("boltzmann " << std::exp(-(delta_U-mu)/temperature));
    // // vesLOG("boltzmann " << boltzmann);
    // // vesLOG("cavity_prob " << cavity_prob);
    // vesLOG("prob "<< prob);
    // vesLOG("return " << std::boolalpha << (std::min(REAL(1),prob) > enhance::random<REAL>(0.0,1.0)));
    // vesLOG("");
    return std::min(REAL(1),prob) > enhance::random<REAL>(0.0,1.0);

}



bool ves::GrandCanonicalDeletion::isValid(const MonteCarloSystem& sys, REAL delta_U) const
{
    static const REAL thermal_wavelength_cubic = std::pow(std::sqrt(REAL(1) / (PI*2*temperature)), 3);

    const REAL prob = thermal_wavelength_cubic * sys.getParticles().get().data.size() * std::exp(-(delta_U+mu)/temperature) / sys.getBox().get().getVolume();

    // vesLOG("");
    // vesLOG("DELETION");
    // vesLOG("delta_U " << delta_U);
    // vesLOG("mu " << mu);
    // vesLOG("boltzmann " << std::exp(-(delta_U+mu)/temperature));
    // // vesLOG("boltzmann " << boltzmann);
    // // vesLOG("cavity_prob " << cavity_prob);
    // vesLOG("prob "<< prob);
    // vesLOG("return " << std::boolalpha << (std::min(REAL(1),prob) > enhance::random<REAL>(0.0,1.0)));
    // vesLOG("");
    return std::min(REAL(1),prob) > enhance::random<REAL>(0.0,1.0);

}