/**
 * \file jet.hpp (description...)
 * 
 * \author A.S. Woodcock
 * 
 * \licence GPL3 (see COPYING.md)
*/

#pragma once
#include "fasterjet/precompiledheader.hpp"

// define struct to represent particle state
struct Particle
{
   float pt = 0; // transverse momentum
   float eta = 0; // pseudo-rapidity
   float phi = 0; // azimuthal angle
};

// define Jet structure
using Jet = Particle;
struct JetDataPrivate;

// stored data for internal use
struct JetData
{
   // revert particle & jet lists to initial version
   void resetState();
   // get the list of particles
   vector<Particle> getParticles();
   // get the list of jets (you should run a jet algorithm first)
   vector<Jet> getJets();
   // data for internal use only
   JetDataPrivate* p = nullptr;
   // allocates private data
   JetData();
   // deallocates private data
   virtual ~JetData();
};

struct JetDataGPU : public JetData
{
   // allocates private data
   JetDataGPU();
};

// function signature of a closest pair finder
using closest_finder_alg = void (*)(JetDataPrivate&);
