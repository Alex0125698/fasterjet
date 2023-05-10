/**
 * \file jetalgorithms.hpp (description...)
 * 
 * \author A.S. Woodcock
 * 
 * \licence GPL3 (see COPYING.md)
*/

#pragma once
#include "fasterjet/precompiledheader.hpp"
#include "fasterjet/jet.hpp"

// see documentation in .cpp

extern void makeRandomParticlesCPU(JetData& d, const int N);
extern void findJetsCPU1(JetData& ioData, const float inR);
extern void findJetsCPU2(JetData& ioData, const float inR);
extern void findJetsCPU3(JetData& ioData, const float inR);
extern void findJetsCPU4(JetData& ioData, const float inR);
extern void findJetsCPU5(JetData& ioData, const float inR);
extern void findJetsCPU6(JetData& ioData, const float inR);
extern void jet_recombiner_CPU(JetDataPrivate& d, const float inR, closest_finder_alg alg);
