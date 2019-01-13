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



namespace enhance
{
    // singleton implementaton
    // to derive class from
    // USE ONLY if DERIVED must not be constructed more than once
    // constructing DERIVED twice is impossible!
    template<typename DERIVED>
    struct Singleton
    {
        // get the singleton instance
        // DERIVED will be constructed on first call
        static DERIVED& getInstance()
        {
            if( INSTANCE == nullptr)
                INSTANCE = new DERIVED();
            return *INSTANCE;
        }
        
        // manually destroy DERIVED
        static void destroyInstance()
        {
            delete INSTANCE;
            INSTANCE = nullptr;
        }
        
    protected:
        // only to derive from
        explicit Singleton() { }
        
        // destroy if derived is destroyed
        virtual ~Singleton()
        {
            INSTANCE = nullptr;
        }
        
    private:
        static DERIVED* INSTANCE;
        

        // do not copy
        Singleton(const Singleton&) = delete;
        Singleton& operator= (const Singleton&) = delete;
    };



    // define Singleton<DERIVED> type outside of class
    template<typename DERIVED> DERIVED* Singleton<DERIVED>::INSTANCE = nullptr;
}