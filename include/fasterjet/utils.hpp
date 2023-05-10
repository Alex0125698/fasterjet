/**
 * \file utils.hpp (description...)
 * 
 * \author A.S. Woodcock
 * 
 * \licence GPL3 (see COPYING.md)
*/

#pragma once
#include "fasterjet/precompiledheader.hpp"

// convert file to string
extern std::string str(const std::ifstream& in);

extern void resetRandom();

// get random number in the interval [a,b)
extern double random(const double a, const double b);

// constants
constexpr float PI = (float)3.141592653589793;
constexpr float MAX_FLOAT = std::numeric_limits<float>::max();

// same as malloc but the address will be an integer multiple of the given alignment
void *aligned_malloc2(size_t alignment, size_t required_bytes);

// aligned free to be used with above
void aligned_free2(void* p);

// square of x
inline float sq(const float x)
{
   return x*x;
}

// find K smallest elements of some array (and corresponding indices). Result is sorted with lowest first.
void findKSmallest(const size_t N, const float* const data, const size_t K, vector<float>& outKSmallest, vector<size_t>& outIndices);

using std::min;
using std::max;
