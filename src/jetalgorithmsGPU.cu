/**
 * \file jetalgorithmsGPU.cu (description...)
 * 
 * \author A.S. Woodcock
 * 
 * \licence GPL3 (see COPYING.md)
*/

#include "fasterjet/precompiledheader.hpp"
#include "fasterjet/jetalgorithmsGPU.hpp"
#include "fasterjet/jetalgorithmsCPU.hpp"
#include "fasterjet/utils.hpp"
#include "fasterjet/jetPrivate.hpp"
#include "fasterjet/utils_cuda.hpp"

// cuda block size
constexpr int BLOCK_SIZE = 256;

// generate N particles with randomized states
void makeRandomParticlesGPU(JetDataGPU& d, const int N)
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

// square of x
__device__ inline float sqr(const float x)
{
   return x*x;
}

// extend JetDataPrivate with GPU-specific stuff

struct JetDataPrivateGPU : public JetDataPrivate
{ 
   void free()
   {
      JetDataPrivate::free();
      if (gpu_particles != nullptr) cudaFree(gpu_particles);
      if (gpu_mindists != nullptr) cudaFree(gpu_mindists);
      if (gpu_others != nullptr) cudaFree(gpu_others);
      gpu_particles = nullptr;
      gpu_mindists = nullptr;
      gpu_others = nullptr;
   }

   void init(const std::vector<Particle>& p)
   {
      JetDataPrivate::init(p);

      if (pnReserve < nReserve)
      {
         // free();
         assert (gpu_particles == nullptr);
         assert (gpu_mindists == nullptr);
         assert (gpu_others == nullptr);
         cudaMalloc(&gpu_particles, nReserve*sizeof(Particle));
         cudaMalloc(&gpu_mindists, nReserve*sizeof(float));
         cudaMalloc(&gpu_others, nReserve*sizeof(int));
         pnReserve = nReserve;
      }
   }

   // fill in gpu_particles from particles
   void copyCPU2GPU()
   {
      cudaMemcpy(gpu_particles, &(particles[0]), nParticles * sizeof(Particle), cudaMemcpyHostToDevice);
   }

   // fill in mindists, others from gpu_mindists, gpu_others
   void copyGPU2CPU()
   {
      cudaMemcpy(&(mindists[0]), gpu_mindists, nParticles * sizeof(float), cudaMemcpyDeviceToHost);
      cudaMemcpy(&(others[0]), gpu_others, nParticles * sizeof(int), cudaMemcpyDeviceToHost);
   }

};

JetDataGPU::JetDataGPU()
{
   p = new JetDataPrivateGPU;
}

// uses asymmetrical distance and separate Particle struct
__global__ void closest_finder_1(const int nParticles, const float R, const Particle* __restrict const particles, 
                                 float* __restrict const mindists, int* __restrict const others)
{
   int index = blockIdx.x * blockDim.x + threadIdx.x;
   if (index >= nParticles) return; 
   
   size_t i = index;

   // initialize with beam axis distance
   mindists[i] = 1./sqr(particles[index].pt);
   others[i] = -1;

   // check all distances to tother particles and get minimum
   // no improvement from this...
   float mindist = mindists[i];
   int other = others[i];
   auto pi = particles[i];

   for (size_t j=0; j<nParticles; ++j)
   {
      auto pj = particles[j];

      float delta_phi = abs(pi.phi - pj.phi); 
      delta_phi = min(delta_phi, 2*PI-delta_phi);
      float sqRij = sqr(pi.eta - pj.eta) + sqr(delta_phi);
      float dist = (sqRij / sqr(R)) / sqr(pi.pt);

      // some improvement...
      bool shouldSwap = dist < mindist && i != j;
      mindist = shouldSwap ? dist : mindist;
      other   = shouldSwap ?    j : other;
   }

   mindists[i] = mindist;
   others[i] = other;
}

// ...
void closest_finder_wrapper(JetDataPrivateGPU& d)
{
      // copy particle array to device
      d.copyCPU2GPU();

      // launch the kernel
      const int N_BLOCKS = 1+(d.nIniParticles + BLOCK_SIZE - 1) / BLOCK_SIZE;
      closest_finder_1<<<N_BLOCKS,BLOCK_SIZE>>>(d.nParticles, d.R, d.gpu_particles, d.gpu_mindists, d.gpu_others);

      // Wait for GPU to finish before accessing on host
      // apparently we dont need it here
      // cudaDeviceSynchronize();

      // copy results to host
      d.copyGPU2CPU();
}


void findJetsGPU1(JetDataGPU& d, const float R)
{
   jet_recombiner_CPU(*d.p,R,(closest_finder_alg)closest_finder_wrapper);
}
