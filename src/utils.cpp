/**
 * \file utils.cpp (description...)
 * 
 * \author A.S. Woodcock
 * 
 * \licence GPL3 (see COPYING.md)
*/

#include "fasterjet/precompiledheader.hpp"
#include "fasterjet/utils.hpp"

std::string str(const std::ifstream& in)
{
	if (!in.is_open()) throw DETAILEDEXCEPTION("could not open file");
	std::stringstream tmp;
	tmp << in.rdbuf();
	return tmp.str();
}

// source of entropy - mersenne twister engine
static std::mt19937 rng(9);

void resetRandom()
{
   rng.seed(9);
}

double random(const double a, const double b)
{
   // seed for the random number generator
   static int seed = 9;


   // get the random numbers
   std::uniform_real_distribution<double> dist(a,b);
   return dist(rng);
}

void *aligned_malloc2(int64 alignment, int64 required_bytes) 
{
    void *p1;
    void **p2;
    int offset=alignment-1+sizeof(void*);
    p1 = malloc(required_bytes + offset);               // the line you are missing
    p2=(void**)(((int64)(p1)+offset)&~(alignment-1));  //line 5
    p2[-1]=p1; //line 6
    return p2;
}

void aligned_free2( void* p ) 
{
    void* p1 = ((void**)p)[-1];         // get the pointer to the buffer we allocated
    free( p1 );
}

void findKSmallest(const int64 N, const float* const data, const int64 K, std::vector<float>& outKSmallest, std::vector<int64>& outIndices)
{
   // Performance: best O(N+K), worst(NK)

   // initialize arrays with worst possible value
   // will keep this array sorted with lowest value first
   outKSmallest.resize(K);
   std::fill(outKSmallest.begin(), outKSmallest.end(), MAX_FLOAT);
   outIndices.resize(K,0);

   // loop over data (N steps)
   for(int64 i=0; i<N; ++i)
   {
      // check if it is smaller than the biggest value in our subset
      int64 k = K-1;
      if (data[i] < outKSmallest[k])
      {
         // figure out which position in the k-array it belongs to (max K steps)
         for (k = K-2; (int64_t)k >= 0 && data[i] < outKSmallest[k]; --k);
         // we dont want to touch data[k] or below so increment once
         ++k;
         // shift everything up 1 position to make space (discarding last one) (max K steps)
         std::copy_backward(outKSmallest.begin()+k,outKSmallest.end()-1,outKSmallest.end());
         std::copy_backward(outIndices.begin()+k,outIndices.end()-1,outIndices.end());
         // copy new value in place
         outKSmallest[k] = data[i];
         outIndices[k] = i;
      }
   }
}
