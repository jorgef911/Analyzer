#ifndef Particle_h
#define Particle_h

// system include files
#include <memory>

// user include files
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <regex>

#include <TTree.h>
#include <TBranch.h>
#include <TLorentzVector.h>

#include "tokenizer.hpp"
#include "Cut_enum.h"

using namespace std;
typedef unsigned int uint;

struct PartStats {
  unordered_map<string,double> dmap;
  unordered_map<string,string> smap;
  unordered_map<string,pair<double,double> > pmap;
  //  unordered_map<string,bool> bmap;
  vector<string> bset;

  bool bfind(string cut) const {
    return find(bset.begin(), bset.end(), cut) != bset.end();
  }
  
};


enum class PType { Electron, Muon, Tau, Jet, FatJet, None};


class Particle {

public:
  Particle();
  Particle(TTree*, string, string, vector<string>);
  virtual ~Particle() {}

  virtual vector<CUTS> findExtraCuts() {return vector<CUTS>();}
  void init();
  void unBranch();
  double pt(uint) const;
  double eta(uint) const;
  double phi(uint) const;
  double energy(uint) const;
  virtual double charge(uint) const;
  TLorentzVector p4(uint) const;
  TLorentzVector& p4(uint);
  TLorentzVector RecoP4(uint) const;
  TLorentzVector& RecoP4(uint);

  uint size() const;
  vector<TLorentzVector>::iterator begin();
  vector<TLorentzVector>::iterator end();
  vector<TLorentzVector>::const_iterator begin() const;
  vector<TLorentzVector>::const_iterator end() const;

  bool needSyst(int) const;

  void addPtEtaPhiESyst(double, double, double, double, int);
  void addP4Syst(TLorentzVector, int);
  void setOrigReco();
  void setCurrentP(int);
  string getName() {return GenName;};

  bool findCut(const vector<string>&, string);
  
  PType type;
  unordered_map<string, PartStats> pstats;
  const map<PType,CUTS> cutMap = {{PType::Electron, CUTS::eGElec}, {PType::Muon, CUTS::eGMuon},
				  {PType::Tau, CUTS::eGTau}};


protected:
  void getPartStats(string);
  TTree* BOOM;
  string GenName;
  unordered_map<CUTS, string, EnumHash> jetNameMap = {
    {CUTS::eRJet1, "Jet1"},               {CUTS::eRJet2, "Jet2"},
    {CUTS::eRCenJet, "CentralJet"},      {CUTS::eRBJet, "BJet"},
    {CUTS::eR1stJet, "FirstLeadingJet"},  {CUTS::eR2ndJet, "SecondLeadingJet"},
    {CUTS::eRWjet, "WJet"}
  };

 private:
  //vector<double>* mpt = 0;
  //vector<double>* meta = 0;
  //vector<double>* mphi = 0;
  //vector<double>* menergy = 0;
  uint  m_n;
  float m_pt[];
  float m_phi[];
  float m_eta[];
  float m_mass[];
  

  vector<TLorentzVector> Reco;
  vector<TLorentzVector> *cur_P;
  vector<string> syst_names;
  vector<vector<TLorentzVector>* > systVec;

  string activeSystematic;
};

class Photon : public Particle {
public:
  Photon();
  Photon(TTree*, string, vector<string>);

  float hoverE[];
  float phoR[];
  float sigmaIEtaIEta[];
  float pfIso_all[];
  float pfIso_chg[];
  bool eleVeto[];
  bool hasPixelSeed[];
};


/////////////////////////////////////////////////////////////////
class Generated : public Particle {

public:
  Generated();
  Generated(TTree*, string, vector<string>);
  
  
  int  genPartIdxMother[];
  int  pdg_id[];
  int  status[];
  int  statusFlags[];
  
};

/////////////////////////////////////////////////////////////////////////
class Jet : public Particle {

public:
  Jet(TTree*, string, vector<string>);

  vector<CUTS> findExtraCuts();
  vector<CUTS> overlapCuts(CUTS);
  bool passedLooseJetID(int);
   
  float area[];
  float bDiscriminator[];
  float chargedEmEnergyFraction[];
  float chargedHadronEnergyFraction[];
  float neutralEmEmEnergyFraction[];
  float neutralHadEnergyFraction[];
  int jetId[];
  int nMuons[];
  int numberOfConstituents[];
  int puID[];
  int partonFlavour[];

 protected:

};

class FatJet : public Particle {

public:
  FatJet(TTree*, string, vector<string>);

  vector<CUTS> findExtraCuts();
  vector<CUTS> overlapCuts(CUTS);
  
  float tau1[];
  float tau2[];
  float tau3[];
  float tau4[];
  float PrunedMass[];
  float SoftDropMass[];

};

class Lepton : public Particle {

public:
  Lepton(TTree*, string, string, vector<string>);

  vector<CUTS> findExtraCuts();

  double charge(uint)const;
  int _charge[];


  virtual bool get_Iso(int, double, double) const {return false;}
};

class Electron : public Lepton {

public:
  Electron(TTree*, string, vector<string>);

  bool get_Iso(int, double, double) const;
  
  bitset<8> cbIDele1;
  bitset<8> cbIDele2;
  bitset<8> cbHLTIDele1;
  bitset<8> cbHLTIDele2;
  
  float miniPFRelIso_all[];
  float miniPFRelIso_chg[];
  float mvaFall17Iso[];
  float mvaFall17noIso[];
  float pfRelIso03_all[];
  float pfRelIso03_chg[];
  bool cutBased[];
  bool cutBased_HLTPreSel[];
  bool mvaIso_90[];
  bool mvanoIso_WP90[];
  bool mvaIso_80[];
  bool mvanoIso_WP80[];
  bool mvaIso_WPL[];
  bool mvanoIso_WPL[];
  bool isPassHEEPId[];

};



class Muon : public Lepton {

public:
  Muon(TTree*, string, vector<string>);

  bool get_Iso(int, double, double) const;

  bool tight[];
  bool soft[];
  float miniPFRelIso_all[];
  float miniPFRelIso_chg[];
  float pfRelIso03_all[];
  float pfRelIso03_chg[];
  float pfRelIso04_all[];
};

class Taus : public Lepton {

public:
  Taus(TTree*, string, vector<string>);

  //  void findExtraCuts();
  vector<CUTS> findExtraCuts();
  bool get_Iso(int, double, double) const;
  bool pass_against_Elec(CUTS, int);
  bool pass_against_Muon(CUTS, int);
  
  bitset<8> tau1minIso;
  bitset<8> tau1maxIso;
  
  bitset<8> tau2minIso;
  bitset<8> tau2maxIso;
  
  bitset<8> tau1ele;
  bitset<8> tau1mu;
  
  bitset<8> tau2ele;
  bitset<8> tau2mu;
  
  
  uint8_t againstElectron[];
  uint8_t againstMuon[];
  bool DecayMode[];
  bool DecayModeNewDMs[];
  //uint8_t Tau_idMVAnewDM2017v2[];
  uint8_t MVAoldDM[];
  uint8_t MVAnewDM[];
  //uint8_t Tau_idMVAoldDM2017v1[];
  //uint8_t Tau_idMVAoldDM2017v2[];
  //uint8_t Tau_idMVAoldDMdR032017v2[];
   
  int decayMode[];
  float leadTkDeltaEta[];
  float leadTkDeltaPhi[];
  float leadTkPtOverTauPt[];
  float dz[];
  float dxy[];
  float chargedIsoPtSum[];
  float neutralIso[];
  float puCorr[];
  
};



#endif
