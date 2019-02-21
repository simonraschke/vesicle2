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

#include <cmath>
#include <array>
#include <algorithm>
#include <numeric>



namespace enhance 
{ 
    template<typename T = float>
    struct Point
    {
        Point(T value);
        Point(T x, T y, T z);

        template<typename U>
        Point(const Point<U>& other);
        
        template<typename U>
        Point<T>& operator=(const Point<U>& other);

        T operator[](std::size_t i);
        
        template<typename U>
        const Point<T> operator+(const Point<U>& other) const;
        
        template<typename U>
        const Point<T> operator-(const Point<U>& other) const;
        
        const Point<T> operator-() const;
        
        template<typename U>
        Point<T>& operator+=(const Point<U>& other);
        
        template<typename U>
        Point<T>& operator-=(const Point<U>& other);

        float norm() const;
        void normalize();
        Point<T> normalized() const;



    protected:
        std::array<T,3> data;
    };



    template<typename T>
    Point<T>::Point(T value)
        : data({value, value, value})
    {

    }



    template<typename T>
    Point<T>::Point(T x, T y, T z)
        : data({x,y,z})
    {

    }



    template<typename T>
    template<typename U>
    Point<T>::Point(const Point<U>& other)
        : data(std::copy_n( std::begin(other.data), 3, std::begin(data)))
    {
        
    }
    


    template<typename T>
    template<typename U>
    Point<T>& Point<T>::operator=(const Point<U>& other)
    {
        data = other.data;
        return data;
    }
    


    template<typename T>
    T Point<T>::operator[](std::size_t i)
    {
        return data[i];
    }
    


    template<typename T>
    template<typename U>
    const Point<T> Point<T>::operator+(const Point<U>& other) const
    {
        Point<T> copy(0);
        std::transform(std::begin(data), std::end(data), std::begin(other.data), std::begin(copy.data), std::plus<T>());
        return copy;
    }
    


    template<typename T>
    template<typename U>
    const Point<T> Point<T>::operator-(const Point<U>& other) const
    {
        Point<T> copy(0);
        std::transform(std::begin(data), std::end(data), std::begin(other.data), std::begin(copy.data), std::minus<T>());
        return copy;
    }
    


    template<typename T>
    const Point<T> Point<T>::operator-() const
    {
        Point<T> copy = *this;
        std::transform(std::begin(data), std::end(data), std::begin(copy.data), std::negate<T>());
        return copy;
    }
    


    template<typename T>
    template<typename U>
    Point<T>& Point<T>::operator+=(const Point<U>& other)
    {
        *this = *this + other;
        return *this;
    }
    


    template<typename T>
    template<typename U>
    Point<T>& Point<T>::operator-=(const Point<U>& other)
    {
        *this = *this - other;
        return *this;
    }



    template<typename T>
    float Point<T>::norm() const
    {
        return std::sqrt(std::inner_product(std::begin(data), std::end(data), std::begin(data), float(0)));
    }



    template<typename T>
    void Point<T>::normalize()
    {
        const T norm_ = norm();
        std::transform(std::begin(data), std::end(data), std::begin(data), [&norm_](T element){ return element/norm_; } );
    }



    template<typename T>
    Point<T> Point<T>::normalized() const
    {
        Point<T> copy = *this;
        copy.normalize();
        return copy;
    }


} // namespace enhance 