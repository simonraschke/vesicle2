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

#include <cmath>
#include <array>
#include <algorithm>
#include <numeric>



namespace enhance 
{ 
    namespace __internal
    {
        template
        <
            typename T, 
            typename ENABLER = typename std::enable_if
            <
                std::is_arithmetic<T>::value
            , void>::type
        >
        struct ArithmeticEnabler
        {
            using element_t = T;
        };
    } // namespace __internal



    template<typename T = float>
    struct Point
        : public __internal::ArithmeticEnabler<T>
    {
        Point();
        Point(T value);
        Point(T x, T y, T z);

        template<typename U>
        Point(const Point<U>& other);

        Point(const Point<T>& other) = default;
        Point(Point<T>&& other) = default;

        constexpr void swap( Point<T>& other ) noexcept(std::is_nothrow_swappable_v<T>);

        Point<T>& operator=(const Point<T>& other) = default;;
        Point<T>& operator=(Point<T> &&) = default;

        inline auto begin()         { return std::begin(data); };
        inline auto end()           { return std::end(data); };
        constexpr inline auto begin() const   { return std::begin(data); };
        constexpr inline auto end()   const   { return std::end(data); };
        constexpr inline auto cbegin() const  { return std::cbegin(data); };
        constexpr inline auto cend()   const  { return std::cend(data); };
        inline auto rbegin()        { return std::rbegin(data); };
        inline auto rend()          { return std::rend(data); };
        constexpr inline auto rbegin() const  { return std::rbegin(data); };
        constexpr inline auto rend()   const  { return std::rend(data); };
        constexpr inline auto crbegin() const { return std::crbegin(data); };
        constexpr inline auto crend()   const { return std::crend(data); };

        // change data element wise directly
        template<typename UNARY_FUNCTOR>
        void unaryApply(const UNARY_FUNCTOR& unary);

        // make a copy of data, change copy and return copy
        template<typename UNARY_FUNCTOR>
        const Point<T> unaryExpr(const UNARY_FUNCTOR& unary) const;

        template<typename U>
        Point<U> cast() const;

        T& operator[](std::size_t i);
        constexpr T operator[](std::size_t i) const;

        T& operator()(std::size_t i);
        constexpr T operator()(std::size_t i) const;
        
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

        template<typename U>
        T dot(const Point<U>& other) const;

        template<typename U>
        Point<T> cross(const Point<U>& other) const;


    protected:
        std::array<T,3> data;
    };



    template<typename T>
    Point<T>::Point()
        : data(
        {
            static_cast<T>(0), 
            static_cast<T>(0), 
            static_cast<T>(0)
        })
    {

    }



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
    Point<U> Point<T>::cast() const
    {
        return Point<U>
        (
            static_cast<U>(data[0]),
            static_cast<U>(data[1]),
            static_cast<U>(data[2])
        );
    }



    template<typename T>
    constexpr void Point<T>::swap( Point<T>& other ) noexcept(std::is_nothrow_swappable_v<T>)
    {
        data.swap(other.data);
    }



    template<typename T>
    template<typename UNARY_FUNCTOR>
    const Point<T> Point<T>::unaryExpr(const UNARY_FUNCTOR& unary) const
    {
        auto copy = Point<T>(*this);
        std::for_each(std::begin(copy.data), std::end(copy.data), unary);
        return copy;
    }



    template<typename T>
    template<typename UNARY_FUNCTOR>
    void Point<T>::unaryApply(const UNARY_FUNCTOR& unary)
    {
        std::for_each(std::begin(data), std::end(data), unary);
    }



    template<typename T>
    T& Point<T>::operator[](std::size_t i)
    {
        return data[i];
    }
    


    template<typename T>
    constexpr T Point<T>::operator[](std::size_t i) const
    {
        return data[i];
    }
    


    template<typename T>
    T& Point<T>::operator()(std::size_t i) 
    {
        return data[i];
    }
    


    template<typename T>
    constexpr T Point<T>::operator()(std::size_t i) const
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
        auto copy = Point<T>(*this);
        copy.normalize();
        return copy;
    }



    template<typename T>
    template<typename U>
    T Point<T>::dot(const Point<U>& other) const
    {
        return std::inner_product(std::begin(data), std::end(data), std::begin(other.data), static_cast<T>(0.0));
    }



    template<typename T>
    template<typename U>
    Point<T> Point<T>::cross(const Point<U>& other) const
    {
        return Point<T>(
            data[1]*other[2] - data[2]*other[1],
            data[2]*other[0] - data[0]*other[2],
            data[0]*other[1] - data[1]*other[0]
        );
    }
} // namespace enhance 



namespace std
{
    template<class T>
    constexpr void swap( enhance::Point<T>& lhs, enhance::Point<T>& rhs ) noexcept(noexcept(lhs.swap(rhs)))
    {
        lhs.swap(rhs);
    }
}