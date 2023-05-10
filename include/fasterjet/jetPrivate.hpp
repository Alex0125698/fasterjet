/**
 * \file jetalgorithmsGPU.hpp (description...)
 * 
 * \author A.S. Woodcock
 * 
 * \licence GPL3 (see COPYING.md)
*/

#pragma once
#include "fasterjet/precompiledheader.hpp"
#include "fasterjet/jet.hpp"

// stored data for internal use
struct JetDataPrivate
{
   // memory alligment for all arrays
   static constexpr size_t ALIGNMENT = 32;

   // todo: make this an input param
   static constexpr float recombinePortion = 0.05;

   // initial number of particles
   size_t nIniParticles = 0;
   // current number of particles
   size_t nParticles = 0;
   // full size of particle arrays
   size_t nReserve = 0;
   // previous full size of particle arrays
   size_t pnReserve = 0;
   // jet "radius"
   float R = 0.0;

   // the inital particle list before recombination
   std::vector<Particle> iniParticles;
   // the current particle list
   Particle* particles = nullptr;
   // mindist list
   float* mindists = nullptr;
   // others list
   int* others = nullptr;
   // list of final jets
   vector<Jet> jets;

   // for some algs its better to have individual phi, eta and pt arrays
   float* v_pt = nullptr;
   float* v_eta = nullptr;
   float* v_phi = nullptr;

   // particle list for the device
   // these are not used for the CPU version
   Particle* gpu_particles = nullptr;
   // mindist list for the device
   float* gpu_mindists = nullptr;
   // others list for the device
   int* gpu_others = nullptr;

   // store the K smallest mins here for each iteration
   std::vector<float> outKSmallest;
   std::vector<size_t> outIndices; 
   std::vector<size_t> doneIndices;

   // free all memory
   virtual void free();
   // init all memory
   virtual void init(const std::vector<Particle>& p);
   // fill in v_pt, v_eta, v_phi
   void copyCPU2CPU();
   // free memory
   virtual ~JetDataPrivate();
};
