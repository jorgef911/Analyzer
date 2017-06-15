#ifndef JETRESOLUTION
#define JETRESOLUTION

#ifndef STANDALONE
#define STANDALONE
#endif

#include "TRandom3.h"
#include <unordered_map>
#include "Particle.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
//#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#pragma GCC diagnostic pop


namespace pxl {
    class Particle;
}

class JetResolution{
    public:
        JetResolution();
        JetResolution( std::unordered_map<string, PartStats> &distats );
        ~JetResolution() {}

        double getJetResolutionCorrFactor(  TLorentzVector& recJet,
                                            TLorentzVector& genJet,
                                            double const npv,
                                            double const rho,
                                            int const updown
                                            );

    private:
        double getResolution( double const pt, double const eta, double const rho, double const npv ) const;
        double getResolutionSF( double const pt, double const eta, double const rho, double const npv, int const updown ) const;

        TRandom3 m_rand;

        JME::JetResolution m_resolutionPt;
        JME::JetResolutionScaleFactor m_resolutionPt_sf;
};

#endif /*JETRESOLUTION*/
