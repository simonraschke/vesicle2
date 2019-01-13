#pragma once

#define REAL float

constexpr const double PI_4  = 0.785398163397448309616;
constexpr const double PI_2  = 1.57079632679489661923;
constexpr const double PI    = 3.14159265358979323846;
constexpr const double TWOPI = 3.14159265358979323846*2;

#define likely(x)      __builtin_expect(!!(x), 1)
#define unlikely(x)    __builtin_expect(!!(x), 0)


#include <csignal>
#include <iostream>
#ifndef NDEBUG
    #define vesDEBUG(x) {std::cerr << "[DEBUG] "; do { std::cerr << x; } while (0); std::cerr << '\n';}
#else
    #define vesDEBUG(x)
#endif
#define vesLOG(x) {std::clog << "[LOG] "; do { std::clog << x; } while (0); std::clog << '\n';}
#define vesWARNING(x) {std::clog << "[WARNING] "; do { std::clog << x; } while (0); std::clog << '\n';}
#define vesCRITICAL(x) {std::cerr << "[ERROR] "<< __FILE__ <<":" << __LINE__ << "  "; do { std::cerr << x; } while (0); std::cerr <<" raising SIGABRT\n"; std::exit(SIGABRT);}