/**
 * \file main.cpp (description...)
 * 
 * \author A.S. Woodcock
 * 
 * \licence GPL3 (see COPYING.md)
*/

#include "fasterjet/precompiledheader.hpp"
#include "fasterjet/jetalgorithmsCPU.hpp"
#include "fasterjet/jetalgorithmsGPU.hpp"
#include "fasterjet/utils.hpp"
#include <iomanip>

// number of particles (actually its 1 less since we include the beam axis)
constexpr size_t N_PARTICLES = 1024;
// "radius" for jet finder
constexpr float R = 1.5;

void runJetFinderCPUSingleCore()
{
   // all data for our algorithms is stored here
   JetDataGPU data;
   // JetData data;

   // initialize particle array randomly
   makeRandomParticlesCPU(data, N_PARTICLES);
   // makeRandomParticlesGPU(data, N_PARTICLES);

   Timer timer;

   auto testIt = [&](std::function<void()> fcn, string name)
   {
      std::cout << "testing " << name <<std::endl;
      const size_t nTrials = 10;
      double minTime = MAX_FLOAT;
      double aveTime = 0;
      for (size_t i = 0; i<nTrials; ++i)
      {
         // reset state
         data.resetState();
         // start timer
         timer.restart();
         // run alg
         fcn();
         // check the duration
         double dur = timer.getDuration();
         minTime = min(minTime,dur);
         aveTime += dur;
      }

      aveTime /= nTrials;
      std::cout << "   run-time: " << std::setprecision(3) << minTime << "(min) " << aveTime << "(ave) seconds" << std::endl;
      std::cout << "   nJets: " << data.getJets().size() << std::endl;
   };

   testIt([&](){findJetsCPU1(data,R);}, "ALG1: cpu, asymmetrical, struct");
   testIt([&](){findJetsCPU2(data,R);}, "ALG2: cpu, asymmetrical, arrays");
   testIt([&](){findJetsCPU3(data,R);}, "ALG3: cpu, symmetrical, struct");
   testIt([&](){findJetsCPU4(data,R);}, "ALG4: cpu, symmetrical, arrays");
   testIt([&](){findJetsCPU5(data,R);}, "ALG5: cpu, avx2, asymmetrical, arrays");
   testIt([&](){findJetsCPU6(data,R);}, "ALG6: cpu, avx2, symmetrical, arrays");
   // testIt([&](){findJetsGPU1(data,R);}, "ALG7: gpu, asymmetrical, arrays");

}

int main()
{
	try
	{
      runJetFinderCPUSingleCore();
	}
	catch(DetailedException& e)
	{
		std::cout << e.getDetails() << std::endl;
	}
	catch(std::exception& e)
	{
		std::cout << "std::exception  " << e.what() << std::endl;
	}
	return 0;
}
