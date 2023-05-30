/**
 * \file jetalgorithmsCPU.cpp (description...)
 * 
 * \author A.S. Woodcock
 * 
 * \licence GPL3 (see COPYING.md)
*/

#include "fasterjet/precompiledheader.hpp"
#include "fasterjet/jetalgorithmsCPU.hpp"
#include "fasterjet/utils.hpp"
#include "fasterjet/utils_avx.hpp"
#include "fasterjet/jetPrivate.hpp"

// WARNING: this is not actually working yet
// change this if you want to try AVX512
using namespace AVX256;
// using namespace AVX512;

// TODO: can unroll inner loop a bit.
// also, j can be stored as an AVX register
// NB: there are 16 avx registers (32 for avx512)

// N.log(N), closest pair finder using divide & conquer (idiot version) 

// N.log(N), closest pair finder using divide & conquer (no recursion or malloc)
// void closest_pair_finder(JetDataPrivate& d)
// {
//    // sort data in x-direction (phi-dir)
//    std::sort(&(d.particles[0]), &(d.particles[d.nParticles])
//         , [&](const Particle& a, const Particle& b){ return a.phi < b.phi; });

//    // assume data comes already sorted

//    // array for storing y-sorted data (eta-dir)
//    std::vector<Particle> yParts;
//    yParts.reserve(max(2ull,d.nParticles/5));

//    // stack-allocated array for storing sizes stack
//    std::array<int64,32> sizes;
//    int64 send = 0;

//    // current global min dist
//    float minDist = MAX_FLOAT;
//    int minj = 0, minother = 0; // the closest pair



//    int64 currSize = 1;

//    while (currSize < d.nParticles)
//    {
//       currSize *= 2;

//       for (int64 i=0; i<d.nParticles; i+=currSize)
//       {
//          int64 start = i;
//          int64 end = min(d.nParticles, i+currSize)-1;

//          // calc mindist for trivial case
//          if (currSize == 2 && start != end)
//          {
//             auto pi = d.particles[start];
//             auto pj = d.particles[end];

//             float delta_phi = abs(pi.phi - pj.phi); 
//             delta_phi = min(delta_phi, 2*PI-delta_phi);
//             float sqRij = sq(pi.eta - pj.eta) + sq(delta_phi);
//             float dist = (sqRij / sq(d.R)) / sq(pi.pt);

//             if (dist < minDist)
//             {
//                minDist = dist;
//                minj = start;
//                minother = end;
//             }
//          }

//          // recombine with previous
//          {
//             int64 pstart = max(0ul,start-currSize);

//             // keep expanding boundary until we hit current mindist
//             int64 i = start, j = start;
//             while (i <= end && abs(d.particles[i].phi-d.particles[start].phi) < minDist) ++i;
//             while (i >= pstart && abs(d.particles[j].phi-d.particles[start].phi) < minDist) --j;

//             // copy this (hopefully short) subset into temporary buffer
//             yParts.clear();
//             std::copy(d.particles[j+1],d.particles[i], yParts.begin());

//             // now sort in the y-direction
//             std::sort(&(d.particles[0]), &(d.particles[d.nParticles])
//                , [&](const Particle& a, const Particle& b){ return a.eta < b.eta; });

//             // geometric arguments imply that we only need to search each set of up
//             // to 8 consecuative points in y-direction
            
//             int64 start = 0;
//             while (start < yParts.size())
//             {
//                int64 end = min(start+8,yParts.size());

//                for (int64 i=start; i<end; ++i)
//                {
//                   for (int64 j=i+1; j<end; ++j)
//                   {
//                      auto pi = d.particles[i];
//                      auto pj = d.particles[j];

//                      float delta_phi = abs(pi.phi - pj.phi); 
//                      delta_phi = min(delta_phi, 2*PI-delta_phi);
//                      float sqRij = sq(pi.eta - pj.eta) + sq(delta_phi);
//                      float dist = (sqRij / sq(d.R)) / sq(pi.pt);


//                      if (dist < minDist)
//                      {
//                         minDist = dist;
//                         minj = start;
//                         minother = end;
//                      }
//                   }
//                }
//                start += 8;
//             }


//          }

//       }
//    }



// }

// // N.log(N), closest pair finder using divide & conquer (perfect)

// generate N particles with randomized states
void makeRandomParticlesCPU(JetData& d, const int N)
{
   vector<Particle> particles(N);
   for (auto& p : particles)
   {
      p.pt = random(1,10);
      p.eta = random(-6,+6);
      p.phi = random(0,2*PI);
   }
   d.p->init(particles);
}

void test(const float R, const int64 nParticles, float* const __restrict mindists, 
          int* const __restrict others, const float* const __restrict v_pt, 
          const float* const __restrict v_phi, const float* const __restrict v_eta)
{
   // initialize with beam axis distance
   for (int64 i=0; i<nParticles; ++i)
   {
      mindists[i] = 1./sq(v_pt[i]);
      others[i] = -1;
   }

   // check all distances to tother particles and get minimum
   for (int64 i=0; i<nParticles; ++i)
   {
      for (int64 j=0; j<nParticles; ++j)
      {
         if (i == j) continue;

         float delta_phi = abs(v_phi[i] - v_phi[j]); 
         delta_phi = min(delta_phi, 2*PI-delta_phi);
         float sqRij = sq(v_eta[i] - v_eta[j]) + sq(delta_phi);
         float dist = (sqRij / sq(R)) / sq(v_pt[i]);

         if (dist < mindists[i])
         {
            mindists[i] = dist;
            others[i] = j;
         }
      }
   }
}

// same except with sort optimization
void closest_finder_1C(JetDataPrivate& d)
{
   static int counter = 0;
   // initialize with beam axis distance
   for (int64 i=0; i<d.nParticles; ++i)
   {
      d.mindists[i] = 1./sq(d.particles[i].pt);
      d.others[i] = -1;
   }

   d.particles[d.nParticles].eta = 1e5;

   // check all distances to tother particles and get minimum
   for (int64 i=0; i<d.nParticles; ++i)
   {


      int delta = 1;
      for (int x = 0; x<2; ++x)
      {
         for (int64 j=i+delta; j<d.nParticles && j >= 0; j+=delta) // j<d.nParticles
         {
            ++counter;
            auto pi = d.particles[i];
            auto pj = d.particles[j];
            float scale = sq(pi.pt * d.R);
            float delta_eta = sq(pi.eta - pj.eta)/scale;

            // check if eta distance is getting to large
            if (delta_eta >= d.mindists[i]) break;

            float delta_phi = abs(pi.phi - pj.phi); 
            delta_phi = sq(min(delta_phi, 2*PI-delta_phi))/scale;
            float dist = delta_phi + delta_eta;

            if (dist < d.mindists[i])
            {
               d.mindists[i] = dist;
               d.others[i] = j;
            }
         }
         delta = -1;
      }
   }
   // if (d.nParticles < 10)
   // std::cout << "counter " << counter << std::endl;
}

// uses asymmetrical distance and separate Particle struct
void closest_finder_1(JetDataPrivate& d)
{
   static int counter = 0;
   // initialize with beam axis distance
   for (int64 i=0; i<d.nParticles; ++i)
   {
      d.mindists[i] = 1./sq(d.particles[i].pt);
      d.others[i] = -1;
   }

   // check all distances to tother particles and get minimum
   for (int64 i=0; i<d.nParticles; ++i)
   {
      for (int64 j=0; j<d.nParticles; ++j)
      {
         if (i == j) continue;

         ++counter;
         auto pi = d.particles[i];
         auto pj = d.particles[j];

         float delta_phi = abs(pi.phi - pj.phi); 
         delta_phi = min(delta_phi, 2*PI-delta_phi);
         float sqRij = sq(pi.eta - pj.eta) + sq(delta_phi);
         float dist = (sqRij / sq(d.R)) / sq(pi.pt);

         if (dist < d.mindists[i])
         {
            d.mindists[i] = dist;
            d.others[i] = j;
         }
      }
   }

   //    if (d.nParticles < 10)
   // std::cout << "counter " << counter << std::endl;
}

// uses asymmetrical distance and separate eta,phi,pt arrays
void closest_finder_2(JetDataPrivate& d)
{
   // put particle array into separate arrays
   d.copyCPU2CPU();

   // test(d.R, d.nParticles, d.mindists, d.others, d.v_pt, d.v_phi, d.v_eta);
   // return;

   // Somehow unpacking these gives a 25% improvement.
   // compiler must have been doing some stupid indirection before.
   // doesn't seem to be aliasing as restrict makes no difference.
   const float R = d.R;
   const int64 nParticles = d.nParticles;
   float* const __restrict mindists = d.mindists;
   int* const __restrict others = d.others;
   const float* const __restrict v_pt = d.v_pt;
   const float* const __restrict v_phi = d.v_phi;
   const float* const __restrict v_eta = d.v_eta;

   // initialize with beam axis distance
   for (int64 i=0; i<nParticles; ++i)
   {
      mindists[i] = 1./sq(v_pt[i]);
      others[i] = -1;
   }

   // check all distances to tother particles and get minimum
   for (int64 i=0; i<nParticles; ++i)
   {
      for (int64 j=0; j<nParticles; ++j)
      {
         if (i == j) continue;

         float delta_phi = abs(v_phi[i] - v_phi[j]); 
         delta_phi = min(delta_phi, 2*PI-delta_phi);
         float sqRij = sq(v_eta[i] - v_eta[j]) + sq(delta_phi);
         float dist = (sqRij / sq(R)) / sq(v_pt[i]);

         if (dist < mindists[i])
         {
            mindists[i] = dist;
            others[i] = j;
         }
      }
   }
}

// uses symmetrical distance and separate Particle struct
void closest_finder_3(JetDataPrivate& d)
{
   // not needed for now
   // for (int64 i=d.nParticles; i<d.nParticles+5; ++i)
   //    d.particles[i].eta = -10000;

   // initialize with worst possible distance
   for (int64 i=0; i<d.nParticles; ++i)
      d.mindists[i] = MAX_FLOAT;

   // check all distances to tother particles and get minimum
   for (int64 i=0; i<d.nParticles; ++i)
   {
      for (int64 j=0; j<i; ++j)
      {
            auto pi = d.particles[i];
            auto pj = d.particles[j];

            float delta_phi = abs(pi.phi - pj.phi); 
            delta_phi = min(delta_phi, 2*PI-delta_phi);
            float sqRij = sq(pi.eta - pj.eta) + sq(delta_phi);
            float dist = (sqRij / sq(d.R)) / sq(max(pj.pt,pi.pt));

            if (dist < d.mindists[i])
            {
               d.mindists[i] = dist;
               d.others[i] = j;
            }
      }
   }

   // update minima based on the fact that min(other,x) = min(min(i,other)=min(other,i),min(other,j))
   for (int64 i=0; i<d.nParticles; ++i)
   {
      int64 other = d.others[i];
      int64 j = d.others[other];

      float min_other_j = d.mindists[other];
      float min_other_i = d.mindists[i]; // = min_i_other

      if (min_other_j > min_other_i)
      {
         d.mindists[other] = min_other_i;
         d.others[other] = i;
      }
   }

   // check is beam dist is closer
   for (int64 i=0; i<d.nParticles; ++i)
   {
      float beamdist = 1./sq(d.particles[i].pt);
      if (beamdist < d.mindists[i])
      {
         d.others[i] = -1;
         d.mindists[i] = beamdist;
      }
   }
}

// uses symmetrical distance and separate eta,phi,pt arrays
void closest_finder_4(JetDataPrivate& d)
{
   // put particle array into separate arrays
   d.copyCPU2CPU();

   // // deal with pointer aliasing. Gives a 25% improvement
   // const float R = d.R;
   // const int64 nParticles = d.nParticles;
   // float* const __restrict mindists = d.mindists;
   // int* const __restrict others = d.others;
   // const float* const __restrict v_pt = d.v_pt;
   // const float* const __restrict v_phi = d.v_phi;
   // const float* const __restrict v_eta = d.v_eta;

   // initialize with beam axis distance
   for (int64 i=0; i<d.nParticles; ++i)
   {
      d.mindists[i] = 1./sq(d.v_pt[i]);
      d.others[i] = -1;
   }

   // check all distances to tother particles and get minimum
   for (int64 i=0; i<d.nParticles; ++i)
   {
      for (int64 j=0; j<i; ++j)
      {
         float delta_phi = abs(d.v_phi[i] - d.v_phi[j]); 
         delta_phi = min(delta_phi, 2*PI-delta_phi);
         float sqRij = sq(d.v_eta[i] - d.v_eta[j]) + sq(delta_phi);
         float dist = (sqRij / sq(d.R)) / sq(max(d.v_pt[j],d.v_pt[i]));

         if (dist < d.mindists[i])
         {
            d.mindists[i] = dist;
            d.others[i] = j;
         }

         if (dist < d.mindists[j])
         {
            d.mindists[j] = dist;
            d.others[j] = i;
         }
      }
   }
}

// uses symmetrical distance and separate eta,phi,pt arrays
void closest_finder_4B(JetDataPrivate& d)
{
   // put particle array into separate arrays
   d.copyCPU2CPU();

   // initialize with worst possible distance
   for (int64 i=0; i<d.nParticles; ++i)
      d.mindists[i] = MAX_FLOAT;

   // check all distances to tother particles and get minimum
   for (int64 i=0; i<d.nParticles; ++i)
   {
      for (int64 j=0; j<i; ++j)
      {
         float delta_phi = abs(d.v_phi[i] - d.v_phi[j]); 
         delta_phi = min(delta_phi, 2*PI-delta_phi);
         float sqRij = sq(d.v_eta[i] - d.v_eta[j]) + sq(delta_phi);
         float dist = (sqRij / sq(d.R)) / sq(max(d.v_pt[j],d.v_pt[i]));

         if (dist < d.mindists[i])
         {
            d.mindists[i] = dist;
            d.others[i] = j;
         }
      }
   }

   // update minima based on the fact that min(other,x) = min(min(i,other)=min(other,i),min(other,j))
   for (int64 i=0; i<d.nParticles; ++i)
   {
      int64 other = d.others[i];
      int64 j = d.others[other];

      float min_other_j = d.mindists[other];
      float min_other_i = d.mindists[i]; // = min_i_other

      if (min_other_j > min_other_i)
      {
         d.mindists[other] = min_other_i;
         d.others[other] = i;
      }
   }

   // check is beam dist is closer
   for (int64 i=0; i<d.nParticles; ++i)
   {
      float beamdist = 1./sq(d.particles[i].pt);
      if (beamdist < d.mindists[i])
      {
         d.others[i] = -1;
         d.mindists[i] = beamdist;
      }
   }
}

// AVX version; uses asymmetrical distance and separate eta,phi,pt arrays
void closest_finder_5(JetDataPrivate& d)
{
   // put particle array into separate arrays
   d.copyCPU2CPU();

   // constants
   vecXf sign_bit = _mm256_set1_ps(-0.0f);
   vecXf two_pi = _mm256_set1_ps(2*PI);
   vecXf R2 = _mm256_set1_ps(sq(d.R));

   for (int64 i=0; i<d.nParticles; i+=DX)
   {
      // Loads a floating-point vector from an aligned memory address
      vecXf v_pt_i = _mm256_load_ps(&(d.v_pt[i]));
      vecXf v_eta_i = _mm256_load_ps(&(d.v_eta[i]));
      vecXf v_phi_i = _mm256_load_ps(&(d.v_phi[i]));

      // mindists[i] = 1./sq(d.v_pt[i]) (initialize to beam axis dist)
      vecXf mindist = _mm256_set1_ps(1.);
      mindist = _mm256_div_ps(mindist, v_pt_i);
      mindist = _mm256_div_ps(mindist, v_pt_i);

      // others[i] = -1 (initialize to beam axis dist)
      vecXi other = _mm256_set1_epi32(-1);

      // pull this out of loop
      vecXf v_pt_i_m4 = _mm256_div_ps(mindist,R2);

      // get |i|i+1|...|i+7| vector
      vecXi tmp = _mm256_set1_epi32(i);
      vecXi iii = _mm256_setr_epi32(0,1,2,3,4,5,6,7);
      iii = _mm256_add_epi32(iii,tmp);
      
      // find the minimum distance
      for (int j=0; j<d.nParticles; ++j)
      {
         // vecXf v_pt_j = _mm256_set1_ps(v_pt[j]);
         vecXf v_phi_j = _mm256_set1_ps(d.v_phi[j]);
         vecXf v_eta_j = _mm256_set1_ps(d.v_eta[j]);

         // float delta_phi = abs(d.v_phi[i] - d.v_phi[j]); 
         vecXf delta_phi = _mm256_sub_ps(v_phi_i, v_phi_j);
         delta_phi = _mm256_andnot_ps(sign_bit, delta_phi);
         
         // if (j == i) continue; see usage in need2swap mask below
         vecXi jjj = _mm256_set1_epi32(j);
         vecXi mask = _mm256_cmpeq_epi32(jjj,iii);

         // delta_phi = min(delta_phi, 2*PI-delta_phi);
         vecXf tmp1 = _mm256_sub_ps(two_pi,delta_phi);
         delta_phi = _mm256_min_ps(delta_phi,tmp1);

         // float sqRij = sq(d.v_eta - d.v_eta) + sq(delta_phi);
         vecXf tmp2 = _mm256_sub_ps(v_eta_i,v_eta_j);
         vecXf sqRij = _mm256_mul_ps(delta_phi,delta_phi);
         sqRij = _mm256_fmadd_ps(tmp2,tmp2,sqRij);

         // float dist = (sqRij / sq(d.R)) / sq(d.v_pt[i]);
         vecXf dist = _mm256_mul_ps(sqRij,v_pt_i_m4);

         // if (dist < mindist && i != j) other = j;
         vecXf need2swap = _mm256_cmp_ps(dist, mindist, _CMP_LT_OQ);
         need2swap = _mm256_andnot_ps(_mm256_castsi256_ps(mask),need2swap);
         other = _mm256_blendv_epi8(other,jjj,_mm256_castps_si256(need2swap));

         // if (dist < mindist && i != j) mindist = dist;
         mindist = _mm256_blendv_ps(mindist,dist,need2swap);
      }

      // store results
      _mm256_store_ps(&(d.mindists[i]),mindist);
      _mm256_store_si256((vecXi*)(&(d.others[i])),other);
   }
}

// AVX version; uses symmetrical distance and separate eta,phi,pt arrays
void closest_finder_6(JetDataPrivate& d)
{
   // put particle array into separate arrays
   d.copyCPU2CPU();

   // constants
   vecXf sign_bit = _mm256_set1_ps(-0.0f);
   vecXf two_pi = _mm256_set1_ps(2*PI);
   vecXf R2 = _mm256_set1_ps(sq(d.R));

   for (int64 i=0; i<d.nParticles; i+=DX)
   {
      // Loads a floating-point vector from an aligned memory address
      vecXf v_pt_i = _mm256_load_ps(&(d.v_pt[i]));
      vecXf v_eta_i = _mm256_load_ps(&(d.v_eta[i]));
      vecXf v_phi_i = _mm256_load_ps(&(d.v_phi[i]));

      // this time we cannot initialize the mindist to beam axis
      // as we are going to change both i and j in same iteration 
      // which have different beam axes
      vecXf mindist = _mm256_set1_ps(MAX_FLOAT);
      vecXi other = _mm256_set1_epi32(0);

      // get |i|i+1|...|i+7| vector
      vecXi tmp = _mm256_set1_epi32(i);
      vecXi iii = _mm256_setr_epi32(0,1,2,3,4,5,6,7);
      iii = _mm256_add_epi32(iii,tmp);

      // find the minimum distance
      for (int j=0; j<i+DX && j<d.nParticles; ++j)
      {
         vecXf v_pt_j = _mm256_set1_ps(d.v_pt[j]);
         v_pt_j = _mm256_max_ps(v_pt_j,v_pt_i);

         // get (1./sq(max))/R2
         vecXf v_pt_max = _mm256_set1_ps(1.);
         v_pt_max = _mm256_div_ps(v_pt_max, v_pt_j);
         v_pt_max = _mm256_div_ps(v_pt_max, v_pt_j);
         v_pt_max = _mm256_div_ps(v_pt_max,R2);

         vecXf v_phi_j = _mm256_set1_ps(d.v_phi[j]);
         vecXf v_eta_j = _mm256_set1_ps(d.v_eta[j]);

         // float delta_phi = abs(d.v_phi[i] - d.v_phi[j]); 
         vecXf delta_phi = _mm256_sub_ps(v_phi_i, v_phi_j);
         delta_phi = _mm256_andnot_ps(sign_bit, delta_phi);
         
         // if (j == i) continue; see usage in need2swap mask below
         vecXi jjj = _mm256_set1_epi32(j);
         vecXi mask = _mm256_cmpeq_epi32(jjj,iii);

         // delta_phi = min(delta_phi, 2*PI-delta_phi);
         vecXf tmp1 = _mm256_sub_ps(two_pi,delta_phi);
         delta_phi = _mm256_min_ps(delta_phi,tmp1);

         // float sqRij = sq(d.v_eta - d.v_eta) + sq(delta_phi);
         vecXf tmp2 = _mm256_sub_ps(v_eta_i,v_eta_j);
         vecXf sqRij = _mm256_mul_ps(delta_phi,delta_phi);
         sqRij = _mm256_fmadd_ps(tmp2,tmp2,sqRij);

         // float dist = (sqRij / sq(d.R)) / sq(d.v_pt[i]);
         vecXf dist = _mm256_mul_ps(sqRij,v_pt_max);

         // if (dist < mindist && i != j) other = j;
         vecXf need2swap = _mm256_cmp_ps(dist, mindist, _CMP_LT_OQ);
         need2swap = _mm256_andnot_ps(_mm256_castsi256_ps(mask),need2swap);
         other = _mm256_blendv_epi8(other,jjj,_mm256_castps_si256(need2swap));

         // if (dist < mindist && i != j) mindist = dist;
         mindist = _mm256_blendv_ps(mindist,dist,need2swap);
      }

      // store results
      _mm256_store_ps(&(d.mindists[i]),mindist);
      _mm256_store_si256((vecXi*)(&(d.others[i])),other);
   }

   // update minima based on the fact that min(other,x) = min(min(i,other)=min(other,i),min(other,j))
   for (int64 i=0; i<d.nParticles; ++i)
   {
      int64 other = d.others[i];
      int64 j = d.others[other];

      float min_other_j = d.mindists[other];
      float min_other_i = d.mindists[i]; // = min_i_other

      if (min_other_j > min_other_i)
      {
         d.mindists[other] = min_other_i;
         d.others[other] = i;
      }
   }

   // check is beam dist is closer
   for (int64 i=0; i<d.nParticles; ++i)
   {
      float beamdist = 1./sq(d.particles[i].pt);
      if (beamdist < d.mindists[i])
      {
         d.others[i] = -1;
         d.mindists[i] = beamdist;
      }
   }

}

// run the recombination steps for a generic closest pair algorithm
void jet_recombiner_CPU(JetDataPrivate& d, const float inR, closest_finder_alg alg)
{
   d.R = inR;

   // keep running until all particles are recombined
   while (d.nParticles > 0)
   {
      // setup
      int64 nRecombineSteps = max(1, (int)(d.nParticles*d.recombinePortion));
      d.doneIndices.clear();

      // run closest pair finder
      alg(d);

      // get the top K closest pairs, with K = nRecombineSteps
      findKSmallest(d.nParticles, d.mindists, nRecombineSteps, d.outKSmallest, d.outIndices);

      // recombine top K closest pairs
      for (int64 k=0; k<nRecombineSteps; ++k)
      {
         float mindist = d.outKSmallest[k];
         int64 imin = d.outIndices[k];
         int64 other = d.others[imin];

         // we know that the i's are unique but the other may be equal
         // to another i or other. So we keep track of which ones we have done
         for (int64 i=0; i<d.doneIndices.size(); ++i) // max 2K steps; can change to O(1)
            if (d.doneIndices[i] == imin || d.doneIndices[i] == other) mindist = -1000.0;
         if (mindist == -1000.0) continue;
         
         d.doneIndices.push_back(imin);
         d.doneIndices.push_back(other);

         // if its the beam axis then start a new jet
         if (other == -1)
         {
            std::cout << ' ' << d.nParticles;
            d.jets.push_back(d.particles[imin]);
            d.particles[imin].phi = -2000.;
            // std::cout << "size " << d.nParticles << std::endl;
         }
         // otherwise recombine and delete old particles
         else
         {
            d.particles[imin].pt = (d.particles[imin].pt + d.particles[other].pt)/2.;
            d.particles[imin].eta = (d.particles[imin].eta + d.particles[other].eta)/2.;

            float deltaPhi = abs(d.particles[imin].phi - d.particles[other].phi);
            if (2*PI - deltaPhi < deltaPhi)
               d.particles[imin].phi = (d.particles[imin].phi + d.particles[other].phi)/2. + PI;
            else
               d.particles[imin].phi = (d.particles[imin].phi + d.particles[other].phi)/2.;
            if (d.particles[imin].phi > 2*PI) d.particles[imin].phi -= 2*PI;
            d.particles[other].phi = -2000.;
         }
      }

      // finally, erase the invalid particles
      int64 i=0, j=d.nParticles-1;
      while(i<j)
      {
         // find next valid particle
         while(d.particles[j].phi == -2000.) --j;
         // find next invalid particle
         while(i<j && d.particles[i].phi != -2000.) ++i;
         // swap the particles
         if (i<j) std::swap(d.particles[i],d.particles[j]);
      }
      
      // chop off any remaining invalid particles (TODO is it needed?)
      while(d.particles[j].phi == -2000.) 
         --j;
      d.nParticles = j+1;
   }

   std::cout << "\n";
}

// run the recombination steps for a generic closest pair algorithm
void jet_recombiner_CPU_sorted(JetDataPrivate& d, const float inR, closest_finder_alg alg)
{
   d.R = inR;
   

   // sort data in the eta-direction, N.log(N) steps
   std::sort(&(d.particles[0]), &(d.particles[d.nParticles])
      , [&](const Particle& a, const Particle& b){ return a.eta < b.eta; });

   // keep running until all particles are recombined
   while (d.nParticles > 0)
   {
      // setup
      int64 nRecombineSteps = max(1, (int)(d.nParticles*d.recombinePortion));
      d.doneIndices.clear();

      // run closest pair finder
      alg(d);

      // get the top K closest pairs, with K = nRecombineSteps
      findKSmallest(d.nParticles, d.mindists, nRecombineSteps, d.outKSmallest, d.outIndices);

      int64 toDeleteCount = 0;

      // recombine top K closest pairs
      for (int64 k=0; k<nRecombineSteps; ++k)
      {
         float mindist = d.outKSmallest[k];
         int64 imin = d.outIndices[k];
         int64 other = d.others[imin];

         // we know that the i's are unique but the other may be equal
         // to another i or other. So we keep track of which ones we have done
         for (int64 i=0; i<d.doneIndices.size(); ++i) // max 2K steps; can change to O(1)
            if (d.doneIndices[i] == imin || d.doneIndices[i] == other) mindist = -1000.0;
         if (mindist == -1000.0) continue;
         
         d.doneIndices.push_back(imin);
         d.doneIndices.push_back(other);
         ++toDeleteCount;

         // if its the beam axis then start a new jet
         if (other == -1)
         {
            std::cout << ' ' << d.nParticles;
            d.jets.push_back(d.particles[imin]);
            d.particles[imin].phi = -2000.;
         }
         // otherwise recombine and delete old particles
         else
         {
            float newpt, neweta, newphi;
            newpt = (d.particles[imin].pt + d.particles[other].pt)/2.;
            neweta = (d.particles[imin].eta + d.particles[other].eta)/2.;

            float deltaPhi = abs(d.particles[imin].phi - d.particles[other].phi);
            if (2*PI - deltaPhi < deltaPhi)
               newphi = (d.particles[imin].phi + d.particles[other].phi)/2. + PI;
            else
               newphi = (d.particles[imin].phi + d.particles[other].phi)/2.;
            if (newphi > 2*PI) newphi -= 2*PI;

            // get the min and max of the two indices
            int64 ilow = min(other,imin);
            int64 ihigh = max(other,imin);

            // // figure out where new particle will sit in array (to remain sorted)
            // int64 pos = ilow+1;
            // while(pos < ihigh && neweta > d.particles[pos].eta) ++pos;

            // // we kept going until pos was gt 
            // --pos;

            // // shift all memory down to make room
            // if (pos != ilow)
            // {
            //    std::move(&(d.particles[ilow+1]), &(d.particles[pos+1]), &(d.particles[ilow]));
            // }

            // insert the new particle
            d.particles[ihigh].pt = newpt;
            d.particles[ihigh].eta = neweta;
            d.particles[ihigh].phi = newphi;
            d.particles[ilow].phi = -2000.;
            // d.particles[ilow].eta = 1e30; // !!!!

            // // shift everything above the 'high' particle down 1 position to eliminate it
            // std::move(&(d.particles[ihigh+1]), &(d.particles[d.nParticles]), &(d.particles[ihigh]));

            // set the new array size
            // d.nParticles -= 1;

            // // check sorted
            // for (int i=0; i<(int)(d.nParticles-1); ++i)
            // {
            //    if (d.particles[i].eta > d.particles[i+1].eta)
            //    {
            //       assert(false);
            //    }
            // }


         }
      }

      // -- eliminate holes
      
      // deal with end condition so that it always stops
      d.particles[d.nParticles].phi = -2000.;
      d.particles[d.nParticles+1].phi = 2000.;

      // get first particle to remove
      int64 swapHere = 0;
      while (d.particles[swapHere].phi != -2000.) ++swapHere;

      // shift each valid range of particles down to replace invalid ones
      int64 i=swapHere;
      for ( ; ; )
      {
         while (d.particles[i].phi == -2000.) ++i;
         int64 valid_start = i;
         if (valid_start >= d.nParticles) break;
         while (d.particles[i].phi != -2000.) ++i;
         int64 valid_end = i;
         std::copy(&(d.particles[valid_start]), &(d.particles[valid_end]), &(d.particles[swapHere]));
         swapHere += valid_end-valid_start;
         if (valid_end >= d.nParticles) break;
      }


      d.nParticles -= toDeleteCount;

      // ----------------------------

      d.particles[d.nParticles].eta = -MAX_FLOAT;

      // the last index of the sorted segment
      int isorted = 0;

      while (true)
      {
         // expand sorted segment until we find next unsorted element
         while (d.particles[isorted].eta <= d.particles[isorted+1].eta) ++isorted;

         // if we reached the end of array then done
         if (isorted+1 >= d.nParticles) break;

         // now go backwards until we find the unsorted element's correct position
         int j = isorted;
         while (j >= 0 && d.particles[j].eta > d.particles[isorted+1].eta) --j;

         // make temporary since we will overwrite it
         Particle tmp = d.particles[isorted+1];

         // now shift everything right to make space
         std::copy_backward(&(d.particles[j+1]), &(d.particles[isorted+1]), &(d.particles[isorted+1+1]));

         // now insert particle into proper position
         d.particles[j+1] = tmp;
      }

      // // sort data in the eta-direction, N.log(N) steps
      // if (d.nParticles > 0)
      //    std::sort(&(d.particles[0]), &(d.particles[d.nParticles])
      //       , [&](const Particle& a, const Particle& b){ return a.eta < b.eta; });

   }

   std::cout << "\n";
}

// wrap 'em up

void findJetsCPU1(JetData& d, const float R) 
{
   jet_recombiner_CPU(*d.p,R,closest_finder_1);
}
void findJetsCPU2(JetData& d, const float R) 
{
   jet_recombiner_CPU(*d.p,R,closest_finder_2);
}
void findJetsCPU3(JetData& d, const float R) 
{
   jet_recombiner_CPU(*d.p,R,closest_finder_3);
}
void findJetsCPU4(JetData& d, const float R) 
{
   // jet_recombiner_CPU(*d.p,R,closest_finder_4);
   jet_recombiner_CPU(*d.p,R,closest_finder_4B);
}
void findJetsCPU5(JetData& d, const float R) 
{
   jet_recombiner_CPU(*d.p,R,closest_finder_5);
}
void findJetsCPU6(JetData& d, const float R) 
{
   jet_recombiner_CPU(*d.p,R,closest_finder_6);
}
void findJetsCPU7(JetData& d, const float R) 
{
   jet_recombiner_CPU_sorted(*d.p,R,closest_finder_1C);
}
