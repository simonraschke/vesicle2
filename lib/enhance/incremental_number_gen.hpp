#pragma once

/*

CONCEPT

derive from Base class to generate numbers
of integral type incrementally

*/


namespace enhance
{

    template<typename DERIVED
        , typename T = std::size_t
        , typename ENABLER = typename std::enable_if<std::is_integral<T>::value>::type
    >
    class IncrementalNumberGenerator
    {
    public:
        typedef DERIVED type;

        // call and construct derived class at once
        // DERIVED()()  <- will result in an increasing integral number on each call
        T operator()() const
        {
            return i++;
        }

        void reset()
        {
            i = 0;
        }

        virtual ~IncrementalNumberGenerator() = default;

    protected:
        static T i;
        // make it inheritance only
        // IncrementalNumberGenerator() = default;
    };

    template<typename DERIVED, typename T, typename ENABLER>
    T IncrementalNumberGenerator<DERIVED,T,ENABLER>::i = 0;
}
