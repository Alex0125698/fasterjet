/**
 * \file jetalgorithmsCPU.hpp (description...)
 * 
 * \author A.S. Woodcock
 * 
 * \licence GPL3 (see COPYING.md)
*/

#pragma once
#include "fasterjet/precompiledheader.hpp"
#include "fasterjet/jet.hpp"

extern void makeRandomParticlesGPU(JetDataGPU& d, const int N);
extern void findJetsGPU1(JetDataGPU& ioData, const float inR);
// extern void findJetsGPU2(JetDataGPU& ioData, const float inR);
// extern void findJetsGPU3(JetDataGPU& ioData, const float inR);
// extern void findJetsGPU4(JetDataGPU& ioData, const float inR);
