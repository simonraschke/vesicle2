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

#include "enhance/singleton.hpp"
#include <atomic>



namespace ves { class GLOBAL; };



class ves::GLOBAL
    : public enhance::Singleton<GLOBAL>
{
public:
    enum class STARTMODE { NEW, RESTART };
    std::atomic<STARTMODE> startmode {STARTMODE::NEW};

    enum class SIMULATIONSTATUS { PREPARATION, RUNNING };
    std::atomic<SIMULATIONSTATUS> simulationstatus {SIMULATIONSTATUS::PREPARATION};

    enum class ENSEMBLE { NVT, uVT };
    std::atomic<ENSEMBLE> ensemble {ENSEMBLE::NVT};

    enum class ACCEPTANCE { METROPOLIS };
    std::atomic<ACCEPTANCE> acceptance {ACCEPTANCE::METROPOLIS};

    enum class SIMULATIONMODE { SA, FGA, OSMOTIC };
    std::atomic<SIMULATIONMODE> simulationmode {SIMULATIONMODE::SA};

    enum class FGAMODE { SPHERE, PLANE };
    std::atomic<FGAMODE> fgamode {FGAMODE::PLANE};
};