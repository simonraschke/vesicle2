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



#include "interaction.hpp"



ves::AngularLennardJonesInteraction::AngularLennardJonesInteraction()
    : kappa (Parameters::getInstance().getOption("system.kappa").as<REAL>())
    , a (1.f + kappa*std::sin(Parameters::getInstance().getOption("system.gamma").as<REAL>()))
    , b (1.f - kappa*std::sin(Parameters::getInstance().getOption("system.gamma").as<REAL>()))
    , c ((cartesian(b,0,0) + cartesian(a-1.f, kappa*std::cos(Parameters::getInstance().getOption("system.gamma").as<REAL>()), 0)).norm())
    , sigma (Parameters::getInstance().getOption("system.ljsigma").as<REAL>())
    , epsilon (Parameters::getInstance().getOption("system.ljepsilon").as<REAL>())
    , cutoff_rez_sq (std::pow(1.f/Parameters::getInstance().getOption("system.cell_min_edge").as<REAL>(), 2))
{
    vesDEBUG(__PRETTY_FUNCTION__);
    vesDEBUG("kappa " << kappa);
    vesDEBUG("a " << a);
    vesDEBUG("b " << b);
    vesDEBUG("c " << c);
    vesDEBUG("sigma " << sigma);
    vesDEBUG("epsilon " << epsilon);
    vesDEBUG("cutoff_rez " << std::sqrt(cutoff_rez_sq));
    vesDEBUG("cutoff_rez_sq " << cutoff_rez_sq);

    box.setLengthX(Parameters::getInstance().getOption("system.box.x").as<REAL>());
    box.setLengthY(Parameters::getInstance().getOption("system.box.y").as<REAL>());
    box.setLengthZ(Parameters::getInstance().getOption("system.box.z").as<REAL>());
}



REAL ves::AngularLennardJonesInteraction::calculate(const particle_t& p1, const particle_t& p2) const
{
    cartesian distance_vec = box.distanceVector(p1, p2);

    const REAL r2 = sigma/distance_vec.squaredNorm();
    if(r2 < cutoff_rez_sq) return 0.f;

    const REAL r6 = r2*r2*r2;

    distance_vec.normalize();
    const cartesian p1_orien_kappa = p1.getOrientation()*kappa/2;
    const cartesian p2_orien_kappa = p2.getOrientation()*kappa/2;

    const REAL chi = 
          std::pow(cartesian( -p1_orien_kappa + distance_vec + p2_orien_kappa ).norm() - a,2)
        + std::pow(cartesian(  p1_orien_kappa + distance_vec - p2_orien_kappa ).norm() - b,2)
        + std::pow(cartesian( -p1_orien_kappa + distance_vec - p2_orien_kappa ).norm() - c,2)
        + std::pow(cartesian(  p1_orien_kappa + distance_vec + p2_orien_kappa ).norm() - c,2);
    
    return 4.f*epsilon*(p1.getLJRejection()*p2.getLJRejection()*r6*r6-(1.f-chi)*p1.getLJAttraction()*p2.getLJAttraction()*r6);
}