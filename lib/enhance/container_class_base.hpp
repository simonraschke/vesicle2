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



#include <iterator>



namespace enhance
{  
    template<typename T>
    struct ContainerBase
    {
        using container_t = T;
        using element_t = typename container_t::value_type;
        using allocator_t = typename container_t::allocator_type;

        inline auto begin()         { return std::begin(data); };
        inline auto end()           { return std::end(data); };
        inline auto begin() const   { return std::begin(data); };
        inline auto end()   const   { return std::end(data); };
        inline auto cbegin() const  { return std::cbegin(data); };
        inline auto cend()   const  { return std::cend(data); };
        inline auto rbegin()        { return std::rbegin(data); };
        inline auto rend()          { return std::rend(data); };
        inline auto rbegin() const  { return std::rbegin(data); };
        inline auto rend()   const  { return std::rend(data); };
        inline auto crbegin() const { return std::crbegin(data); };
        inline auto crend()   const { return std::crend(data); };

        container_t data;
    };
}