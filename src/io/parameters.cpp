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



#include "parameters.hpp"



void ves::Parameters::read(int argc, const char* argv[])
{
    vesDEBUG(__PRETTY_FUNCTION__);
    
    namespace po = boost::program_options;

    optionsMap.clear();
    PATH config_file_name;
    
    po::options_description generalOptions("General Options");
    generalOptions.add_options()
        ("config", po::value<PATH>(&config_file_name), "read from config file")
        ("help,h", "show help")
        ("general.cpu_threads", po::value<std::size_t>()->default_value(tbb::task_scheduler_init::default_num_threads()), "maximum number of threads")
        ("general.acceptance", po::value<std::string>()->default_value("metropolis"), "[metropolis]")
        ("general.ensemble", po::value<std::string>()->default_value("NVT"), "[NVT, uVT]")
        ("general.simulation_mode", po::value<std::string>()->default_value("SA"), "[SA,FGA,OSMOTIC]")
        ("general.fga_mode", po::value<std::string>()->default_value("plane"), "[sphere,plane,tube,pair]")
    ;
    
    po::options_description systemOptions("System Options");
    systemOptions.add_options()
        ("system.mobile", po::value<std::size_t>(), "number of mobile particles")

        ("system.guiding_elements_each", po::value<std::size_t>()->default_value(0), "number of guiding elements per frame guide")
        ("system.guiding_elements_restriction", po::value<std::string>()->default_value("inplace"), "restriction of guiding elemets [inplace,structure]")
        // ("system.guiding_elements_distance", po::value<REAL>(), "preferred distance of guiding elements")
        ("system.frame_guides_grid_edge", po::value<std::size_t>()->default_value(1), "number of frame guides per dimension")
        ("system.plane_edge", po::value<REAL>()->default_value(12), "edge length of planar frame-guide or distance of pair")
        
        ("system.osmotic_density_inside", po::value<REAL>()->default_value(0.01), "density of osmotic particles inside GUV")
        ("system.osmotic_density_outside", po::value<REAL>()->default_value(0.01), "density of osmotic particles in bulk")
        
        ("system.density", po::value<REAL>(), "particle density")
        ("system.mu", po::value<REAL>()->default_value(0), "chemical potential")
        ("system.widom_skip", po::value<std::size_t>()->default_value(1), "chemical potential")
        ("system.box.x", po::value<REAL>(), "box edge x")
        ("system.box.y", po::value<REAL>(), "box edge y")
        ("system.box.z", po::value<REAL>(), "box edge z")
        ("system.temperature", po::value<REAL>(), "temperature")
        ("system.time_max", po::value<std::size_t>()->default_value(std::numeric_limits<std::size_t>::max()), "time step")
        ("system.ljepsilon", po::value<REAL>()->default_value(1.0), "kappa")
        ("system.ljsigma", po::value<REAL>()->default_value(1.0), "kappa")
        ("system.kappa", po::value<REAL>()->default_value(1.0), "kappa")
        ("system.gamma", po::value<REAL>()->default_value(11), "gamma angle")

        ("system.sw_position_min", po::value<REAL>()->default_value(0.05), "position stepwidth min")
        ("system.sw_position_max", po::value<REAL>()->default_value(1.0), "position stepwidth min")
        ("system.acceptance_position_target", po::value<REAL>()->default_value(0.3), "position acceptance percentage target")
        ("system.sw_orientation_min", po::value<REAL>()->default_value(M_PI_4/20), "orientation stepwidth min")
        ("system.sw_orientation_max", po::value<REAL>()->default_value(M_PI_4), "orientation stepwidth max")
        ("system.acceptance_orientation_target", po::value<REAL>()->default_value(0.3), "orientation acceptance percentage target")

        ("system.cell_min_edge", po::value<REAL>()->default_value(3.0), "minimum edge length of cell")
        ("system.max_cells_dim", po::value<std::size_t>()->default_value(20), "maximum cells per dimension")
    ;
    
    po::options_description outputOptions("Output Options");
    outputOptions.add_options()
        ("output.path",  po::value<PATH>()->default_value("trajectory.h5"), "trajectory.h5")
        ("output.skip",  po::value<std::size_t>()->default_value(10000), "print every .. steps")
        ("output.suppress_gro", po::bool_switch()->default_value(true), "no .gro output")
    ;
    
    po::options_description inputOptions("Input Options");
    inputOptions.add_options()
        ("input.path",  po::value<PATH>()->default_value("trajectory.h5"), "trajectory.h5")
    ;

    po::options_description allOptions;
    allOptions.add(generalOptions).add(systemOptions).add(outputOptions).add(inputOptions);

    po::store(po::command_line_parser(argc,argv).options(allOptions).run(),optionsMap);
    po::notify(optionsMap);

    PATH config_file_full_path = fs::current_path() / config_file_name;

    if(optionsMap.count("help"))
    {
        std::clog << '\n' << generalOptions;
        std::clog << '\n' << systemOptions;
        std::clog << '\n' << outputOptions;
        std::clog << '\n' << inputOptions;
        std::exit(EXIT_SUCCESS);
    }
    else if(fs::exists(config_file_full_path) && !config_file_name.empty())
    {
        try
        {
            read_from_file(allOptions,optionsMap);
        }
        catch(boost::bad_any_cast& e)
        {
            vesCRITICAL(e.what())
            std::exit(EXIT_FAILURE);
        }
    }

    std::map<std::string, po::variable_value>& m = optionsMap;
    m["system.gamma"] = po::variable_value(enhance::deg_to_rad(optionsMap["system.gamma"].as<REAL>()), false);

    ves::GLOBAL::getInstance().simulationstatus = ves::GLOBAL::SIMULATIONSTATUS::PREPARATION;
    ves::GLOBAL::getInstance().acceptance = ves::GLOBAL::ACCEPTANCE::METROPOLIS;

    if(fs::exists(getOption("input.path").as<PATH>()))
        ves::GLOBAL::getInstance().startmode = ves::GLOBAL::STARTMODE::RESTART;
    else
        ves::GLOBAL::getInstance().startmode = ves::GLOBAL::STARTMODE::NEW;


    if(getOption("general.ensemble").as<std::string>() == "NVT")
        ves::GLOBAL::getInstance().ensemble = ves::GLOBAL::ENSEMBLE::NVT;
    else if(getOption("general.ensemble").as<std::string>() == "uVT")
        ves::GLOBAL::getInstance().ensemble = ves::GLOBAL::ENSEMBLE::uVT;
    else
        vesWARNING("unable to get general.ensemble, setting to NVT");


    if(getOption("general.fga_mode").as<std::string>() == "sphere")
        ves::GLOBAL::getInstance().fgamode = ves::GLOBAL::FGAMODE::SPHERE;
    else if(getOption("general.fga_mode").as<std::string>() == "plane")
        ves::GLOBAL::getInstance().fgamode = ves::GLOBAL::FGAMODE::PLANE;
    else if(getOption("general.fga_mode").as<std::string>() == "tube")
        ves::GLOBAL::getInstance().fgamode = ves::GLOBAL::FGAMODE::TUBE;
    else if(getOption("general.fga_mode").as<std::string>() == "pair")
    {
        ves::GLOBAL::getInstance().fgamode = ves::GLOBAL::FGAMODE::PAIR;
        m["system.gamma"] = po::variable_value(REAL(0.0), false);
        m["system.guiding_elements_each"] = po::variable_value(std::size_t(2), false);
        m["system.guiding_elements_restriction"] = po::variable_value(std::string("inplace"), false);
        m["system.frame_guides_grid_edge"] = po::variable_value(std::size_t(1), false);
    }
    else
        vesWARNING("unable to get general.fga_mode, setting to plane");


    if(getOption("general.simulation_mode").as<std::string>() == "SA")
        ves::GLOBAL::getInstance().simulationmode = ves::GLOBAL::SIMULATIONMODE::SA;
    else if(getOption("general.simulation_mode").as<std::string>() == "FGA")
        ves::GLOBAL::getInstance().simulationmode = ves::GLOBAL::SIMULATIONMODE::FGA;
    else if(getOption("general.simulation_mode").as<std::string>() == "OSMOTIC")
        ves::GLOBAL::getInstance().simulationmode = ves::GLOBAL::SIMULATIONMODE::OSMOTIC;
    else
        vesWARNING("unable to get general.simulation_mode, setting to self assembly");
}



void ves::Parameters::read_from_file(boost::program_options::options_description& desc, boost::program_options::variables_map& vm )
{
    namespace po = boost::program_options;

    IFSTREAM INPUT( vm["config"].as<PATH>() );
    
    vm.clear();
    
    po::store(po::parse_config_file(INPUT,desc),vm);
    po::notify(vm);
}