#include "random.hpp"


namespace enhance
{
    RandomEngineInit::RandomEngineInit()
        : seed(true_engine())
    {
        pseudo_engine.seed(seed);
    }



    RandomEngineInit::RandomEngineInit(const int __seed)
        : seed(__seed)
    {
        pseudo_engine.seed(seed);
    }



    int RandomEngineInit::getSeed() const
    {
        return seed;
    }



    /*
    * specializations of enhance::random<T>(const T,const T)
    */


    // template<>
    // float random(const float a, const float b)
    // {
    //     std::uniform_real_distribution<float> dist(a,b);
    //     return dist(RandomEngine.pseudo_engine);
    // }



    // template<>
    // double random(const double a, const double b)
    // {
    //     std::uniform_real_distribution<double> dist(a,b);
    //     return dist(RandomEngine.pseudo_engine);
    // }



    // template<>
    // unsigned int random(const unsigned int a, const unsigned int b)
    // {
    //     std::uniform_int_distribution<unsigned int> dist(a,b);
    //     return dist(RandomEngine.pseudo_engine);
    // }



    // template<>
    // std::size_t random(const std::size_t a, const std::size_t b)
    // {
    //     std::uniform_int_distribution<std::size_t> dist(a,b);
    //     return dist(RandomEngine.pseudo_engine);
    // }
}
