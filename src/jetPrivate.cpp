/**
 * \file jetPrivate.cpp (description...)
 * 
 * \author A.S. Woodcock
 * 
 * \licence GPL3 (see COPYING.md)
*/

#include "fasterjet/precompiledheader.hpp"
#include "fasterjet/jetPrivate.hpp"
#include "fasterjet/utils.hpp"

// free all data
void JetDataPrivate::free()
{
   if (particles != nullptr) aligned_free2(particles);
   if (mindists != nullptr) aligned_free2(mindists);
   if (others != nullptr) aligned_free2(others);
   if (v_pt != nullptr) aligned_free2(v_pt);
   if (v_eta != nullptr) aligned_free2(v_eta);
   if (v_phi != nullptr) aligned_free2(v_phi);
   particles = nullptr;
   mindists = nullptr;
   others = nullptr;
   v_pt = nullptr;
   v_eta = nullptr;
   v_phi = nullptr;
}

// init all data
void JetDataPrivate::init(const std::vector<Particle>& p)
{
   iniParticles = p;

   // set the sizes
   nIniParticles = p.size();
   nParticles = nIniParticles;

   // leave some headroom so that writes past
   // the end of array dont have any impact
   constexpr int64 HEADROOM = 10*32;
   
   // get full size of arrays
   int64 fullSize = HEADROOM + p.size();

   // allocate the memory if required
   if (nReserve < fullSize)
   {
      free();

      particles = (Particle*) aligned_malloc2(ALIGNMENT, fullSize*sizeof(particles[0]));
      mindists = (float*) aligned_malloc2(ALIGNMENT, fullSize*sizeof(mindists[0]));
      others = (int*) aligned_malloc2(ALIGNMENT, fullSize*sizeof(others[0]));
      v_pt = (float*) aligned_malloc2(ALIGNMENT, fullSize*sizeof(v_pt[0]));
      v_eta = (float*) aligned_malloc2(ALIGNMENT, fullSize*sizeof(v_eta[0]));
      v_phi = (float*) aligned_malloc2(ALIGNMENT, fullSize*sizeof(v_phi[0]));
      jets.reserve(fullSize);

      pnReserve = nReserve;
      nReserve = fullSize;
   }

   // copy in the new particle data
   for (int64 i=0; i<p.size(); ++i)
   {
      particles[i].eta = p[i].eta;
      particles[i].phi = p[i].phi;
      particles[i].pt = p[i].pt;
   }

   outKSmallest.reserve(1+(int)(nIniParticles*recombinePortion));
   outIndices.reserve(1+(int)(nIniParticles*recombinePortion));
   doneIndices.reserve(2*(1+(int)(nIniParticles*recombinePortion)));

   outKSmallest.clear();
   outIndices.clear();
   doneIndices.clear();

   // clear any old data
   jets.clear();
}

// fill in v_pt, v_eta, v_phi
void JetDataPrivate::copyCPU2CPU()
{
   for (int64 i=0; i<nParticles; ++i)
   {
      v_pt[i] = particles[i].pt;
      v_eta[i] = particles[i].eta;
      v_phi[i] = particles[i].phi;
   }
}

// automatically free data
JetDataPrivate::~JetDataPrivate()
{
   free();
}


// --- JetData Implementation ---


// allocate internal memory
JetData::JetData()
{
   p = new JetDataPrivate;
}

// free internal memory
JetData::~JetData()
{
   if (p != nullptr)
   {
      delete p;
      p = nullptr;
   }
}

// reset to initial state before recombination
void JetData::resetState()
{
   p->init(p->iniParticles);
}

// get particle array
vector<Particle> JetData::getParticles()
{
   return p->iniParticles;
}

// get final jets
vector<Jet> JetData::getJets()
{
   return p->jets;
}
