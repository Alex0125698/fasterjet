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
#include "fastjet/ClusterSequence.hh"

bool compareFloat(const float a, const float b)
{
    return abs(a-b) < 1e-4;
}

fastjet::PseudoJet partA_to_partB(const Particle& in)
{
    float eta = in.eta;
    float pt = in.pt;
    float phi = in.phi;

    float px = pt*cos(phi);
    float py = pt*sin(phi);

    float theta = 2.*atan(std::exp(-eta));

    float pz = pt / tan(theta);


    auto result = fastjet::PseudoJet(px,py,pz,0.0);

    float eta2 = result.pseudorapidity();
    float pt2 = result.perp();
    float phi2 = result.phi();

    assert(compareFloat(eta,eta2));
    assert(compareFloat(pt,pt2));
    assert(compareFloat(phi,phi2));

    return result;
}

// int testFastjet () {
//     cout << " pt y phi" << endl;
//     for (unsigned i = 0; i < jets.size(); i++) {
//         cout << "jet " << i << ": "<< jets[i].pt() << " "
//         << jets[i].rap() << " " << jets[i].phi() << endl;
//         vector<PseudoJet> constituents = jets[i].constituents();
//         for (unsigned j = 0; j < constituents.size(); j++) {
//             cout << " constituent " << j << "'s pt: "<< constituents[j].pt() << endl;
//         }
//     }
// }

void findJetsFJ1(JetDataGPU& data, const float R)
{
    using namespace fastjet;
    std::vector<PseudoJet> particles;
    auto pin = data.getParticles();
    for (auto& p : pin) particles.push_back(partA_to_partB(p));

    JetDefinition jet_def(antikt_algorithm, R);

    ClusterSequence cs(particles, jet_def);

    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

    std::cout << "nJets " << jets.size() << std::endl;
}

// number of particles
constexpr int64 N_PARTICLES = 512;
// "radius" for jet finder
constexpr float R = 1.5;

void runJetFinderCPUSingleCore()
{
   // all data for our algorithms is stored here
   JetDataGPU data;

   // initialize particle array randomly
   // makeRandomParticlesCPU(data, N_PARTICLES);
   makeRandomParticlesGPU(data, N_PARTICLES);

   Timer timer;

   auto testIt = [&](std::function<void()> fcn, string name)
   {
      std::cout << "\ntesting " << name <<std::endl;
      const int64 nTrials = 10;
      double minTime = MAX_FLOAT;
      double aveTime = 0;
      for (int64 i = 0; i<nTrials; ++i)
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

   testIt([&](){findJetsFJ1(data,R);}, "fastjet"); 
   testIt([&](){findJetsCPU1(data,R);}, "ALG1: cpu, asymmetrical, struct");
   // testIt([&](){findJetsCPU2(data,R);}, "ALG2: cpu, asymmetrical, arrays");
   // testIt([&](){findJetsCPU3(data,R);}, "ALG3: cpu, symmetrical, struct");
   // testIt([&](){findJetsCPU4(data,R);}, "ALG4: cpu, symmetrical, arrays");
   testIt([&](){findJetsCPU5(data,R);}, "ALG5: cpu, avx2, asymmetrical, arrays");
   // testIt([&](){findJetsCPU6(data,R);}, "ALG6: cpu, avx2, symmetrical, arrays");
   testIt([&](){findJetsCPU7(data,R);}, "ALG7: cpu, asymmetrical, struct, sorted"); 
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
