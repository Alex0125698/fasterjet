/**
 * \file utils_cuda.hpp (description...)
 * 
 * \author A.S. Woodcock
 * 
 * \licence GPL3 (see COPYING.md)
*/

// set whether cuda is enabled for this project
// TODO: use cmake config for this
#define WITH_CUDA

// if using cuda then just include the appropriate header
#ifdef WITH_CUDA

#include <cuda_runtime.h>

#endif
