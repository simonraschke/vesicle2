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

#include <iostream>
#include <sstream>
#include <string>
#include <list>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/iter_find.hpp>


namespace enhance
{
    template<typename Bindable_type>
    std::string streamBindableToString( Bindable_type& bindable )
    {
        return static_cast<std::ostringstream&>(std::ostringstream().seekp(0) << bindable).str();
    }



    template<typename Bindable_type>
    const char* streamBindableToChar( Bindable_type& bindable )
    {
        return static_cast<std::ostringstream&>(std::ostringstream().seekp(0) << bindable).str().c_str();
    }



    template <typename... Args>
    #ifdef __clang_major__
    inline std::string toStringViaStream(Args&&... args)
    #elif  __GNUC__
    constexpr inline std::string toStringViaStream(Args&&... args)
    #else
        #error no valid compiler
    #endif
    {
        std::stringstream stream;
        using expander = int[];
        (void) expander {0, (void(stream << std::forward<Args>(args)),0)...};
        return stream.str();
    }



    inline std::deque<std::string> splitAtDelimiter(std::string str, std::string delim)
    {
        std::deque<std::string> stringList;
        boost::iter_split(stringList, str, boost::first_finder(delim));
        return stringList;
    }
    
}
