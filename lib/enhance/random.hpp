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

#include <random>



namespace enhance
{
    // The static random Engine to be allocated once
    // if the application is shared memory parallelized
    //      this will not yield the same random numbers on every pull
    // std::random_device is to generate a real random number for the seed
    static struct RandomEngineInit
    {
        // construct and set seed via std::random_device
        explicit RandomEngineInit();

        // construct and set seed manually
        RandomEngineInit(const int);

        // get the seed
        int getSeed() const;

        // the pseudo random engine to be used publicly
        std::mt19937_64 pseudo_engine {};

    private:
        std::random_device true_engine {};
        const int seed;

    // call via enhance::RandomEngine.pseudo_engine
    } thread_local RandomEngine;



    // partially specialized template function to get random numbers
    template<typename T>
    T random(const T, const T);
}