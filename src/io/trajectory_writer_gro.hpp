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
#include "parameters.hpp"
#include "particles/base.hpp"
#include "enhance/output_utility.hpp"
#include <filesystem>
#include <memory>
#include <iomanip>



namespace ves 
{ 
    struct TrajectoryWriterGro;
    struct MonteCarloSystem;
};



struct ves::TrajectoryWriterGro
{
    using PATH = Parameters::PATH;
    using FSTREAM = std::ofstream;

    void setup()
    {
        file_path = std::filesystem::path( *(enhance::splitAtDelimiter(file_path.string(), ".").rbegin()+1)+"."+filetype );
        vesDEBUG(__PRETTY_FUNCTION__<< "  " << file_path);
        if(FILE.is_open())
        {
            FILE.close();
        }

        // if(std::filesystem::exists(file_path))
        // {
        //     // splitting the filepath
        //     auto string_parts = enhance::splitAtDelimiter(file_path.string(), ".");

        //     std::string first_part = std::accumulate(string_parts.begin(), std::next(string_parts.rbegin()).base(), std::string(""), [](auto i, auto j){return i+j+".";});
        //     {
        //         std::string name_appendix = "_old";
        //         first_part.insert(std::next(first_part.rbegin()).base(), std::begin(name_appendix), std::end(name_appendix));
        //     }
        //     // std::string filetype = "gro";

        //     // making "trajectory_old.gro" from "trajectory.gro"
        //     PATH destination = std::filesystem::absolute(std::filesystem::path(first_part+filetype));
        //     vesLOG("trajectory file " << file_path.string() << " already exists. will backup to " << destination.string());

        //     // and backup the old trajectory
        //     std::filesystem::copy_file(file_path, destination, std::filesystem::copy_options::overwrite_existing);
        // }

        if(GLOBAL::getInstance().startmode == GLOBAL::STARTMODE::NEW)
        {
            vesLOG("trunc open " << file_path.string());
            FILE.open(file_path, std::ios_base::out);
        }
        else
        {
            vesLOG("app open " << file_path.string());
            FILE.open(file_path, std::ios_base::app);
        }

        FILE.setf(std::ios::fixed | std::ios::showpoint);
    }



    template<typename SYSTEM>
    void write(const SYSTEM& sys)
    {
        const auto lines = sys.getParticles().get().template numType<Particle::TYPE::MOBILE>()*2
                         + sys.getParticles().get().template numType<Particle::TYPE::FRAME>()*2
                         + sys.getParticles().get().template numType<Particle::TYPE::OSMOTIC>();
        const auto kappa = ves::Parameters::getInstance().getOption("system.kappa").as<REAL>();

        FILE << "t=" << std::fixed << sys.getTime() << '\n';
        FILE << lines << '\n';
        
        unsigned long atom = 0;
        for( typename std::remove_const<decltype(lines)>::type residue = 0; residue < lines/2; ++residue)
        {
            const Particle::Base& particle = std::cref(*sys.getParticles().get().data[residue]);
            const auto coords = sys.getBox().get().scaleDown(particle.getCoordinates());
            const auto orientation = /*particle.getType() == Particle::TYPE::OSMOTIC ? decltype(coords)(0,0,0) :*/ particle.getOrientation();
            const auto name = particle.getName();
            
            static std::map<ves::Particle::TYPE, std::string> colorupmap;
            colorupmap[ves::Particle::TYPE::MOBILE] = "S";
            colorupmap[ves::Particle::TYPE::FRAME] = "F";
            colorupmap[ves::Particle::TYPE::OSMOTIC] = "O";

            static std::map<ves::Particle::TYPE, std::string> colordownmap;
            colordownmap[ves::Particle::TYPE::MOBILE] = "C";
            colordownmap[ves::Particle::TYPE::FRAME] = "F";

            if(particle.getType() == Particle::TYPE::OSMOTIC)
            {
                FILE << std::setw(5) <<  residue+1;
                FILE << std::setw(5) <<  particle.getName();
                FILE << std::setw(5) <<  colorupmap.at(particle.getType());
                FILE << std::setw(5) <<  atom+1;
                FILE << std::setprecision(3);
                FILE << std::setw(8) <<  coords(0);
                FILE << std::setw(8) <<  coords(1);
                FILE << std::setw(8) <<  coords(2);
                FILE << std::setprecision(4);
                FILE << std::setw(8) <<  0.f;
                FILE << std::setw(8) <<  0.f;
                FILE << std::setw(8) <<  0.f;
                FILE << '\n';
                ++atom;
            }
            else
            {
                FILE << std::setw(5) <<  residue+1;
                FILE << std::setw(5) <<  particle.getName();
                FILE << std::setw(5) <<  colorupmap.at(particle.getType());
                FILE << std::setw(5) <<  atom+1;
                FILE << std::setprecision(3);
                FILE << std::setw(8) <<  coords(0) + orientation(0)*kappa/2.f;
                FILE << std::setw(8) <<  coords(1) + orientation(1)*kappa/2.f;
                FILE << std::setw(8) <<  coords(2) + orientation(2)*kappa/2.f;
                FILE << std::setprecision(4);
                FILE << std::setw(8) <<  0.f;
                FILE << std::setw(8) <<  0.f;
                FILE << std::setw(8) <<  0.f;
                FILE << '\n';
                ++atom;

                FILE << std::setw(5) <<  residue+1;
                FILE << std::setw(5) <<  particle.getName();
                FILE << std::setw(5) <<  colordownmap.at(particle.getType());
                FILE << std::setw(5) <<  atom+1;
                FILE << std::setprecision(3);
                FILE << std::setw(8) <<  coords(0) - orientation(0)*kappa/2.f;
                FILE << std::setw(8) <<  coords(1) - orientation(1)*kappa/2.f;
                FILE << std::setw(8) <<  coords(2) - orientation(2)*kappa/2.f;
                FILE << std::setprecision(4);
                FILE << std::setw(8) <<  0.f;
                FILE << std::setw(8) <<  0.f;
                FILE << std::setw(8) <<  0.f;
                FILE << '\n';
                ++atom;
            }
        }

        FILE << sys.getBox().get().getLengthX() << ' ' << sys.getBox().get().getLengthY() << ' ' << sys.getBox().get().getLengthZ();
        FILE << '\n';

        makeStartFileVMD(sys);
    }


    template<typename SYSTEM>
    void makeStartFileVMD(const SYSTEM&);
    

protected:
    PATH working_dir {std::filesystem::current_path()};
    PATH file_path {ves::Parameters::getInstance().getOption("output.path").as<PATH>()};
    FSTREAM FILE {};

    const std::string filetype = {"gro"};
};



template<typename SYSTEM>
void ves::TrajectoryWriterGro::makeStartFileVMD(const SYSTEM& sys)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    // OFSTREAM STARTER;
    // STARTER.open(vesicle.vmd);
    // STARTER << "#!/bin/bash" << '\n';
    // STARTER << "vmd" << '\n';
    // STARTER.close();

    const auto lines = sys.getParticles().get().template numType<Particle::TYPE::MOBILE>()*2
                     + sys.getParticles().get().template numType<Particle::TYPE::FRAME>()*2;

    FSTREAM VMD;
    VMD.open(".vmdrc");
    VMD << "mol load gro " << file_path.filename() << '\n';
    VMD << "light 0 on" << '\n';
    VMD << "light 1 on" << '\n';
    VMD << "light 2 on" << '\n';
    VMD << "light 3 on" << '\n';
    VMD << "display nearclip set 0" << '\n';
    
    VMD << "axes location off" << '\n';
    VMD << "stage location off" << '\n';
    
    VMD << "menu main on" << '\n';
    VMD << "menu graphics on" << '\n';
    
    VMD << "display resize 1080 1080" << '\n';
    VMD << "display reposition" << '\n';
    VMD << "display projection perspective" << '\n';
    VMD << "display rendermode GLSL" << '\n';
    VMD << "display cuedensity 0.12" << '\n';
    VMD << "color Display Background white" << '\n';
    // VMD << "draw color black" << '\n';
    VMD << "mol coloring 7 2 ResName" << '\n';
    VMD << "mol modstyle 0 0 Licorice 2 10 30" << '\n';
    VMD << "mol modmaterial 0 0 AOChalky" << '\n';
    
    VMD << "# color definitions" << '\n';
    VMD << "color change rgb  0 0.07 0.20 0.48 ;# blue" << '\n';
    VMD << "color change rgb  1 0.70 0.20 0.10 ;# red" << '\n';
    VMD << "color change rgb  2 0.40 0.40 0.40 ;# gray" << '\n';
    VMD << "color change rgb  3 0.70 0.40 0.00 ;# orange" << '\n';
    VMD << "color change rgb  4 0.80 0.70 0.10 ;# yellow" << '\n';
    VMD << "color change rgb  7 0.13 0.47 0.04 ;# green" << '\n';
    VMD << "color change rgb  8 1.00 1.00 1.00 ;# white" << '\n';
    VMD << "color change rgb 10 0.10 0.70 0.80 ;# cyan" << '\n';
    VMD << "color change rgb 11 0.60 0.10 0.60 ;# purple" << '\n';
    VMD << "color change rgb 16 0.15 0.15 0.15 ;# black" << '\n';

    VMD << "after idle {" << '\n';
    VMD << "  pbc box -color black -width 1" << '\n';
    VMD << "  # set colors" << '\n';
    VMD << "  # create dummy molecule with one atom" << '\n';
    VMD << "  set mol [mol new atoms 1]" << '\n';
    VMD << "  set sel [atomselect $mol all]" << '\n';
    VMD << "  # add items to color categories" << '\n';
    // VMD << "  $sel set name A" << '\n';
    // VMD << "  $sel set type A" << '\n';
    // VMD << "  $sel set name B" << '\n';
    // VMD << "  $sel set type B" << '\n';
    // VMD << "  $sel set name C" << '\n';
    // VMD << "  $sel set type C" << '\n';
    VMD << "  # now we can define colors" << '\n';
    VMD << "  color Name F 1" << '\n';
    VMD << "  color Type F 1" << '\n';
    VMD << "  color Name S 23" << '\n';
    VMD << "  color Type S 23" << '\n';
    VMD << "  color Name C 16" << '\n';
    VMD << "  color Type C 16" << '\n';
    VMD << "  color Name O 2" << '\n';
    VMD << "  color Type O 2" << '\n';
    VMD << "  mol delete $mol" << '\n';

    VMD << "  # clean up" << '\n';
    VMD << "  $sel delete" << '\n';
    // VMD << "  set mol [mol new atoms 1]" << '\n';
    // VMD << "  set sel [atomselect $mol all]" << '\n';
    VMD << "  for {set x 0} {$x < " << lines <<"} {incr x} {" << '\n';
    VMD << "    set y [expr $x+1]" << '\n';
    VMD << "    set sel [atomselect top \"index $x $y\"]" << '\n';
    VMD << "    incr x 1" << '\n';
    VMD << "    set bonds [$sel getbonds]" << '\n';
    VMD << "    set ids [$sel get index]" << '\n';
    VMD << "    lassign $bonds atom1bonds atom2bonds" << '\n';
    VMD << "    lassign $ids atom1id atom2id" << '\n';
    VMD << "    lappend atom1bonds $atom2id" << '\n';
    VMD << "    lappend atom2bonds $atom1id" << '\n';
    VMD << "    $sel setbonds [list $atom1bonds $atom2bonds]"  << '\n';
    VMD << "    }" << '\n';
    
    VMD << "}" << '\n';
    
    VMD.close();
}