#include "JetResolution.h"
#include <sstream>
#include <cmath>


JetResolution::JetResolution( std::unordered_map<string, PartStats> &distats ) :
    m_rand( 0 ),
    //m_resolutionPt( distats["Jet_systematics"].smap["jetResolutionFileUnmatched"] ),
    //m_resolutionPt_sf( distats["Jet_systematics"].smap["jetResolutionSFFile" ]  )
    m_resolutionPt( "Pileup/Spring16_25nsV6_MC_PtResolution_AK4PFchs.txt" ),
    m_resolutionPt_sf( "Pileup/Spring16_25nsV6_MC_SF_AK4PFchs.txt"  )
{
}


double JetResolution::getResolution( double const pt, double const eta, double const rho, double const npv ) const {
   JME::JetParameters parameters;
   parameters.setJetPt( pt );
   parameters.setJetEta( eta );
   parameters.setRho( rho );
   parameters.setNPV( npv );
   return m_resolutionPt.getResolution( parameters );
}

double JetResolution::getResolutionSF( double const pt, double const eta, double const rho, double const npv, int const updown ) const {
   JME::JetParameters parameters;
   parameters.setJetPt( pt );
   parameters.setJetEta( eta );
   parameters.setRho( rho );
   parameters.setNPV( npv );
   if ( updown == 0 )
      return m_resolutionPt_sf.getScaleFactor( parameters );
   else if ( updown == -1 )
      return m_resolutionPt_sf.getScaleFactor( parameters, Variation::UP );
   else if ( updown == 1 )
      return m_resolutionPt_sf.getScaleFactor( parameters, Variation::DOWN );
   else {
        std::stringstream err;
        err << "In file " << __FILE__ << ", function " << __func__ << ", line " << __LINE__ << std::endl;
        err << "Variaton of jet resolution scale factor for " << updown << " not possible. Choose:" << std::endl;
        err << "\t 0 - no variation\n\t 1 - scale up\n\t-1 - scale down " << std::endl;
        return 0;
   }
}

double JetResolution::getJetResolutionCorrFactor(
                                            TLorentzVector& recJet,
                                            TLorentzVector& genJet,
                                            double const npv=0,
                                            double const rho=0,
                                            int const updown=0
                                                        ) {
    double const recJetPt  = recJet.Pt();
    double recJetEta = recJet.Eta();
    //no eta>5 possible set it to 5 times the sign of eta
    if ( fabs( recJetEta ) >= 5 ) {
        recJetEta = 4.99 * recJetEta / fabs( recJetEta );
    }

    double scaling_factor = 1.0;
    scaling_factor = getResolutionSF( recJetPt, recJetEta, rho, npv, updown );

    double jetCorrFactor = 1.0;

    // Found a match?
    if( genJet == TLorentzVector(0,0,0,0) ) {
        double const genJetPt  = genJet.Pt();
        double const corrJetPt = std::max( 0.0, genJetPt +
                                                             scaling_factor *
                                                             ( recJetPt - genJetPt )
                                                             );

        jetCorrFactor = corrJetPt / recJetPt;

    // If not, just smear with a Gaussian.
    } else {

        double const sigma_MC = getResolution( recJetPt, recJetEta, rho, npv );

        double const sigma = sigma_MC *
                                    std::sqrt( scaling_factor * scaling_factor - 1.0 );
        double const corrJetPt = m_rand.Gaus( recJetPt, sigma );

        jetCorrFactor = corrJetPt / recJetPt;
    }

    // WARNING: 0 can be returned! Catch this case at the place this function is
    // used!
    return jetCorrFactor;
}
