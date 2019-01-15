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



#include "trajectory_writer_h5.hpp"




void ves::TrajectoryWriterH5::setup()
{
    file_path = std::filesystem::path( *(enhance::splitAtDelimiter(file_path.string(), ".").rbegin()+1)+"."+filetype );
    vesDEBUG(__PRETTY_FUNCTION__<< "  " << file_path);
    try
    {
        h5file.close();
    }
    catch(...)
    {
        vesWARNING("unable to close h5file");
    }

    if(std::filesystem::exists(file_path))
    {
        // splitting the filepath
        auto string_parts = enhance::splitAtDelimiter(file_path.string(), ".");

        std::string first_part = std::accumulate(string_parts.begin(), std::next(string_parts.rbegin()).base(), std::string(""), [](auto i, auto j){return i+j+".";});
        {
            std::string name_appendix = "_old";
            first_part.insert(std::next(first_part.rbegin()).base(), std::begin(name_appendix), std::end(name_appendix));
        }
        // std::string filetype = "gro";

        // making "trajectory_old.gro" from "trajectory.gro"
        PATH destination = std::filesystem::absolute(std::filesystem::path(first_part+filetype));
        vesLOG("trajectory file " << file_path.string() << " already exists. will backup to " << destination.string());

        // and backup the old trajectory
        std::filesystem::copy_file(file_path, destination, std::filesystem::copy_options::overwrite_existing);
    }

    if(GLOBAL::getInstance().startmode == GLOBAL::STARTMODE::NEW)
    {
        vesLOG("trunc open " << file_path.string());
        h5file.open(file_path.string(), h5xx::file::trunc);
    }
    else
    {
        vesLOG("app open" << file_path.string());
        h5file.open(file_path.string(), h5xx::file::out);
    }
}