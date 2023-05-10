/**
 * \file precompiledheader.hpp (description...)
 * 
 * \author A.S. Woodcock
 * 
 * \licence GPL3 (see COPYING.md)
*/

#pragma once
#define NOMINMAX // fix Windows Bullshit

// --- C Standard Library ---

#include <cassert>
//#include <cctype>
//#include <cerrno>
//#include <cfenv>
//#include <cfloat>
//#include <cinttypes>
//#include <ciso646>
//#include <climits>
//#include <clocale>
#include <cmath>
//#include <csetjmp>
//#include <csignal>
//#include <cstdarg>
//#include <cstdbool>
//#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
//#include <ctgmath>
//#include <ctime>
//#include <cuchar>
//#include <cwchar>
//#include <cwctype>


// --- C++ Standard Library ---

// - Containers
#include <array>
// #include <deque>
//#include <forward_list>
//#include <list>
// #include <map>
//#include <queue>
//#include <set>
//#include <stack>
// #include <unordered_map>
// #include <unordered_set>
#include <vector>
#include <initializer_list>

// - Input/Output
#include <fstream>
// #include <iomanip>
//#include <ios>
//#include <iosfwd>
#include <iostream>
// #include <istream>
// #include <ostream>
#include <sstream>

// - Multi-threading
#include <atomic>
#include <condition_variable>
//#include <future>
#include <mutex>
#include <thread>

// - Other
#include <numbers>
#include <algorithm>
// #include <bitset>
#include <chrono>
//#include <codecvt>
// #include <complex>
#include <exception>
#include <functional>
#include <iterator>
#include <limits>
//#include <locale>
#include <memory>
//#include <new>
//#include <numeric>
#include <random>
//#include <ratio>
//#include <regex>
#include <stdexcept>
#include <string>
//#include <system_error>
#include <tuple>
//#include <typeindex>
// #include <typeinfo>
// #include <type_traits>
// #include <utility>
//#include <valarray>
// #include <variant>


// --- My Stuff ---

#include "fasterjet/DetailedException.hpp"
#include "fasterjet/timer.hpp"
// #include "FixedArray.h"
// #include "forwardDeclarations.h"


// -- other important libraries ---

// #define EIGEN_DONT_VECTORIZE
// #include <Eigen/Core>

// --- So that we don't go insane ---

using namespace std::numbers;

template <typename T>
inline auto str(T t)
{
	return std::to_string(t);
}
using std::pair;
using std::string;
using std::shared_ptr;
using std::unique_ptr;
using std::vector;
using std::array;

using int8 = int8_t;
using uint8 = uint8_t;
using int16 = int16_t;
using uint16 = uint16_t;
using int32 = int32_t;
using uint32 = uint32_t;
using int64 = int64_t;
using uint64 = uint64_t;
using uint = unsigned int;

// enum Color
// {
// 	BLACK,
// 	BLUE,
// 	GREEN,
// 	AQUA,
// 	RED,
// 	PURPLE,
// 	YELLOW,
// 	WHITE,
// 	GREY,
// 	LIGHT_BLUE,
// 	LIGHT_GREEN,
// 	LIGHT_AQUA,
// 	LIGHT_RED,
// 	LIGHT_PURPLE,
// 	LIGHT_YELLOW,
// 	BRIGHT_WHITE
// };

// void setStroke(Color color = Color::WHITE);
// void setFill(Color color = Color::BLACK);
