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

#include "common/definitions.hpp"
#include "enhance/singleton.hpp"
#include <fstream>
#include <filesystem>
#include <boost/program_options.hpp>
#include <tbb/task_scheduler_init.h>



namespace ves { struct Parameters; }

struct ves::Parameters
    : public enhance::Singleton<Parameters>
{
    using PATH = std::filesystem::path;
    using IFSTREAM = std::ifstream;
    using OFSTREAM = std::ofstream;
    
    void read(int, const char* []);

    const auto& getOption(const std::string& s) const { return optionsMap[s]; }
    const boost::program_options::variables_map & getOptions() const { return optionsMap; }

protected:
    boost::program_options::variables_map optionsMap {};

    static void read_from_file(boost::program_options::options_description&, boost::program_options::variables_map&);
};