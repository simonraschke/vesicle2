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

#include <cstdlib>
#include <exception>
#include <iostream>
#include <tbb/spin_mutex.h>
#include <tbb/atomic.h>



namespace ves { class StepwidhtAlignmentUnit; }



class ves::StepwidhtAlignmentUnit
{
public:
    void setup(std::size_t, float, float, float = 0.3);

    float operator()() const;
    void accepted();
    void rejected();
    float getTarget() const;
    std::size_t getAccepted() const;
    std::size_t getRejected() const;
    float getRatio() const;

protected:
    void alignment_check();
    void do_aligment();

    float stepwidth;
    float stepwidth_min;
    float stepwidth_max;
    float ratio_target;
    float ratio_old {0};
    tbb::atomic<std::size_t> accepted_count {0};
    tbb::atomic<std::size_t> rejected_count {0};
    std::size_t alignment_every {0};
    bool setup_flag {false};
    tbb::spin_mutex mutex;
};