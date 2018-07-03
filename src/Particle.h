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

const Int_t MAXINDEX = 500;

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
  float m_pt[MAXINDEX];
  float m_phi[MAXINDEX];
  float m_eta[MAXINDEX];
  float m_mass[MAXINDEX];
  

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

  float hoverE[MAXINDEX];
  float phoR[MAXINDEX];
  float sigmaIEtaIEta[MAXINDEX];
  float pfIso_all[MAXINDEX];
  float pfIso_chg[MAXINDEX];
  bool eleVeto[MAXINDEX];
  bool hasPixelSeed[MAXINDEX];
};


/////////////////////////////////////////////////////////////////
class Generated : public Particle {

public:
  Generated();
  Generated(TTree*, string, vector<string>);
  
  
  int  genPartIdxMother[MAXINDEX];
  int  pdg_id[MAXINDEX];
  int  status[MAXINDEX];
  int  statusFlags[MAXINDEX];
  
};

/////////////////////////////////////////////////////////////////////////
class Jet : public Particle {

public:
  Jet(TTree*, string, vector<string>);

  vector<CUTS> findExtraCuts();
  vector<CUTS> overlapCuts(CUTS);
  bool passedLooseJetID(int);
   
  float area[MAXINDEX];
  float bDiscriminator[MAXINDEX];
  float chargedEmEnergyFraction[MAXINDEX];
  float chargedHadronEnergyFraction[MAXINDEX];
  float neutralEmEmEnergyFraction[MAXINDEX];
  float neutralHadEnergyFraction[MAXINDEX];
  int jetId[MAXINDEX];
  int nMuons[MAXINDEX];
  int numberOfConstituents[MAXINDEX];
  int puID[MAXINDEX];
  int partonFlavour[MAXINDEX];

 protected:

};

class FatJet : public Particle {

public:
  FatJet(TTree*, string, vector<string>);

  vector<CUTS> findExtraCuts();
  vector<CUTS> overlapCuts(CUTS);
  
  float tau1[MAXINDEX];
  float tau2[MAXINDEX];
  float tau3[MAXINDEX];
  float tau4[MAXINDEX];
  float PrunedMass[MAXINDEX];
  float SoftDropMass[MAXINDEX];

};

class Lepton : public Particle {

public:
  Lepton(TTree*, string, string, vector<string>);

  vector<CUTS> findExtraCuts();

  double charge(uint)const;
  int _charge[MAXINDEX];


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
  
  float miniPFRelIso_all[MAXINDEX];
  float miniPFRelIso_chg[MAXINDEX];
  float mvaFall17Iso[MAXINDEX];
  float mvaFall17noIso[MAXINDEX];
  float pfRelIso03_all[MAXINDEX];
  float pfRelIso03_chg[MAXINDEX];
  bool cutBased[MAXINDEX];
  bool cutBased_HLTPreSel[MAXINDEX];
  bool mvaIso_90[MAXINDEX];
  bool mvanoIso_WP90[MAXINDEX];
  bool mvaIso_80[MAXINDEX];
  bool mvanoIso_WP80[MAXINDEX];
  bool mvaIso_WPL[MAXINDEX];
  bool mvanoIso_WPL[MAXINDEX];
  bool isPassHEEPId[MAXINDEX];

};



class Muon : public Lepton {

public:
  Muon(TTree*, string, vector<string>);

  bool get_Iso(int, double, double) const;

  bool tight[MAXINDEX];
  bool soft[MAXINDEX];
  float miniPFRelIso_all[MAXINDEX];
  float miniPFRelIso_chg[MAXINDEX];
  float pfRelIso03_all[MAXINDEX];
  float pfRelIso03_chg[MAXINDEX];
  float pfRelIso04_all[MAXINDEX];
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
  
  
  uint8_t againstElectron[MAXINDEX];
  uint8_t againstMuon[MAXINDEX];
  bool DecayMode[MAXINDEX];
  bool DecayModeNewDMs[MAXINDEX];
  //uint8_t Tau_idMVAnewDM2017v2[MAXINDEX];
  uint8_t MVAoldDM[MAXINDEX];
  uint8_t MVAnewDM[MAXINDEX];
  //uint8_t Tau_idMVAoldDM2017v1[MAXINDEX];
  //uint8_t Tau_idMVAoldDM2017v2[MAXINDEX];
  //uint8_t Tau_idMVAoldDMdR032017v2[MAXINDEX];
   
  int decayMode[MAXINDEX];
  float leadTkDeltaEta[MAXINDEX];
  float leadTkDeltaPhi[MAXINDEX];
  float leadTkPtOverTauPt[MAXINDEX];
  float dz[MAXINDEX];
  float dxy[MAXINDEX];
  float chargedIsoPtSum[MAXINDEX];
  float neutralIso[MAXINDEX];
  float puCorr[MAXINDEX];
  
};



#endif
