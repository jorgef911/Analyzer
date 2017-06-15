#ifndef Systematics_hh
#define Systematics_hh

#include <string>
#include <iostream>
#include <functional>
#include <unordered_map>
#include "Particle.h"
//#include <boost/unordered_map.hpp>

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#pragma GCC diagnostic pop


class TRandom3;

class Systematics {

public:
  Systematics();
  Systematics(std::unordered_map<string, PartStats> const &distats);
  ~Systematics();

  void init();

  void shiftParticle(Particle& jet, TLorentzVector recJet, double const& ratio, double& dPx, double& dPy, string syst);

private:
};
#endif /*Systematics_hh*/
