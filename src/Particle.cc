#include "Particle.h"
#include <signal.h>

#define SetBranch(name, variable) BOOM->SetBranchStatus(name, 1);  BOOM->SetBranchAddress(name, &variable);

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////    PARTICLE   ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////



//particle is a objet that stores multiple versions of the particle candidates
Particle::Particle(TTree* _BOOM, string _GenName, string filename, vector<string> _syst_names) : BOOM(_BOOM), GenName(_GenName), syst_names(_syst_names) {
  type = PType::None;
  getPartStats(filename);

  for( auto item : syst_names) {
    systVec.push_back(new vector<TLorentzVector>());
  }

  //set the pt to an empty vector if the branch does not exist
  //backward compatible yeah!!
  if( _BOOM->GetListOfBranches()->FindObject((GenName+"_pt").c_str()) ==0){
    mpt=new vector<double>();
    cout<<"no "<<GenName<<" will deactivate branch!"<<endl;
  }else{
    SetBranch((GenName+"_pt").c_str(), mpt);
    SetBranch((GenName+"_eta").c_str(), meta);
    SetBranch((GenName+"_phi").c_str(), mphi);
    SetBranch((GenName+"_energy").c_str(), menergy);
  }

  //  activeSystematic="orig";
}

double Particle::pt(uint index)const         {return cur_P->at(index).Pt();}
double Particle::eta(uint index)const        {return cur_P->at(index).Eta();}
double Particle::phi(uint index)const        {return cur_P->at(index).Phi();}
double Particle::energy(uint index)const     {return cur_P->at(index).E();}
double Particle::charge(uint index)const     {return 0;}

uint Particle::size()const                   {return Reco.size();}
vector<TLorentzVector>::iterator Particle::begin(){ return cur_P->begin();}
vector<TLorentzVector>::iterator Particle::end(){ return cur_P->end();}
vector<TLorentzVector>::const_iterator Particle::begin()const { return cur_P->begin();}
vector<TLorentzVector>::const_iterator Particle::end()const { return cur_P->end();}

TLorentzVector Particle::p4(uint index)const {return (cur_P->at(index));}
TLorentzVector& Particle::p4(uint index) {return cur_P->at(index);}
TLorentzVector Particle::RecoP4(uint index)const {return Reco.at(index);}
TLorentzVector& Particle::RecoP4(uint index) {return Reco.at(index);}




void Particle::addPtEtaPhiESyst(double ipt,double ieta, double iphi, double ienergy, int syst){
  TLorentzVector mp4;
  mp4.SetPtEtaPhiE(ipt,ieta,iphi,ienergy);
  systVec.at(syst)->push_back(mp4);
}


void Particle::addP4Syst(TLorentzVector mp4, int syst){
  systVec.at(syst)->push_back(mp4);
}


void Particle::init(){
    //cleanup of the particles
  Reco.clear();
  for(auto it: systVec){
    it->clear();
  }
  TLorentzVector tmp;
  for(uint i=0; i < mpt->size(); i++) {
    tmp.SetPtEtaPhiE(mpt->at(i),meta->at(i),mphi->at(i),menergy->at(i));
    Reco.push_back(tmp);

  }
  setCurrentP(-1);

}


void Particle::setOrigReco() {
  /////memory loss here if no smear and new vector. Only once, so ignore for now...
  systVec.at(0) = &Reco;
}


void Particle::setCurrentP(int syst){
  if(syst == -1) {
    cur_P = &Reco;
  } else if( systVec.at(syst)->size() == 0) {
    cur_P = systVec.at(0);  //orig
  } else {
    cur_P = systVec.at(syst);
  }
  //  activeSystematic=syst;
}


void Particle::unBranch() {
  BOOM->SetBranchStatus((GenName+"*").c_str(), 0);
}


void Particle::getPartStats(string filename) {
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  ifstream info_file(filename);
  boost::char_separator<char> sep(", \t");

  if(!info_file) {
    std::cout << "could not open file " << filename <<std::endl;
    return;
  }

  vector<string> stemp;
  string group,line;
  while(getline(info_file, line)) {
    tokenizer tokens(line, sep);
    stemp.clear();
    for(tokenizer::iterator iter = tokens.begin();iter != tokens.end(); iter++) {
      if( ((*iter)[0] == '/' && (*iter)[0] == '/') || ((*iter)[0] == '#') ) break;
      stemp.push_back(*iter);
    }
    if(stemp.size() == 0) continue;
    else if(stemp.size() == 1) {
      group = stemp[0];
      continue;
    } else if(group == "") {
      cout << "error in " << filename << "; no groups specified for data" << endl;
      exit(1);
    } else if(stemp.size() == 2) {

      if(stemp[1] == "1" || stemp[1] == "true" ) pstats[group].bset.insert(stemp[0]);
      //      else if(stemp[1] == "0"  || stemp[1] == "false" ) pstats[group]bmap[stemp[0]]=false;

      else if(stemp[1].find_first_not_of("0123456789+-.") == string::npos) pstats[group].dmap[stemp[0]]=stod(stemp[1]);
      else pstats[group].smap[stemp[0]] = stemp[1];
    } else  pstats[group].pmap[stemp[0]] = make_pair(stod(stemp[1]), stod(stemp[2]));
  }
  info_file.close();
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////    PHOTON  //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


Photon::Photon(TTree* _BOOM, string filename, vector<string> syst_names) : Particle(_BOOM, "Photon", filename, syst_names) {
  SetBranch("Photon_et", et);
  SetBranch("Photon_HoverE", hoverE);
  SetBranch("Photon_phoR9", phoR);
  SetBranch("Photon_SigmaIEtaIEta", sigmaIEtaIEta);
  SetBranch("Photon_SigmaIPhiIPhi", sigmaIPhiIPhi);
  SetBranch("Photon_PFChIso", pfChIso);
  SetBranch("Photon_PFPhoIso", pfPhoIso);
  SetBranch("Photon_PFNeuIso", pfNeuIso);
  SetBranch("Photon_EleVeto", eleVeto);
  SetBranch("Photon_hasPixelSeed", hasPixelSeed);
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////    GENERATED   ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


Generated::Generated(TTree* _BOOM, string filename, vector<string> syst_names) : Particle(_BOOM, "Gen", filename, syst_names) {

  SetBranch("Gen_pdg_id", pdg_id);
  SetBranch("Gen_motherpdg_id", motherpdg_id);
  SetBranch("Gen_status", status);
  SetBranch("Gen_BmotherIndex", BmotherIndex);
}


///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////    JET  ////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


Jet::Jet(TTree* _BOOM, string filename, vector<string> syst_names) : Particle(_BOOM, "Jet", filename, syst_names) {
  type = PType::Jet;
  SetBranch("Jet_neutralHadEnergyFraction", neutralHadEnergyFraction);
  SetBranch("Jet_neutralEmEmEnergyFraction", neutralEmEmEnergyFraction);
  SetBranch("Jet_numberOfConstituents", numberOfConstituents);
  SetBranch("Jet_muonEnergyFraction", muonEnergyFraction);
  SetBranch("Jet_chargedHadronEnergyFraction", chargedHadronEnergyFraction);
  SetBranch("Jet_chargedMultiplicity", chargedMultiplicity);
  SetBranch("Jet_chargedEmEnergyFraction", chargedEmEnergyFraction);
  SetBranch("Jet_partonFlavour", partonFlavour);
  SetBranch("Jet_bDiscriminator_pfCISVV2", bDiscriminator);
}

vector<CUTS> Jet::findExtraCuts() {
  vector<CUTS> return_vec;
  if(pstats["Smear"].bset.find("SmearTheJet") != pstats["Smear"].bset.end()) {
    return_vec.push_back(CUTS::eGen);
  }
  return return_vec;
}

vector<CUTS> Jet::overlapCuts(CUTS ePos) {
  vector<CUTS> returnCuts;
  unordered_set<string>& tmpset = pstats[jetNameMap.at(ePos)].bset;
  if(tmpset.find("RemoveOverlapWithMuon1s") != tmpset.end()) returnCuts.push_back(CUTS::eRMuon1);
  if(tmpset.find("RemoveOverlapWithMuon2s") != tmpset.end()) returnCuts.push_back(CUTS::eRMuon2);
  if(tmpset.find("RemoveOverlapWithElectron1s") != tmpset.end()) returnCuts.push_back(CUTS::eRElec1);
  if(tmpset.find("RemoveOverlapWithElectron2s") != tmpset.end()) returnCuts.push_back(CUTS::eRElec2);
  if(tmpset.find("RemoveOverlapWithTau1s") != tmpset.end()) returnCuts.push_back(CUTS::eRTau1);
  if(tmpset.find("RemoveOverlapWithTau2s") != tmpset.end()) returnCuts.push_back(CUTS::eRTau2);

  return returnCuts;
}


///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////    FATJET   ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


FatJet::FatJet(TTree* _BOOM, string filename, vector<string> syst_names) : Particle(_BOOM, "Jet_toptag", filename, syst_names) {
  type = PType::FatJet;
  SetBranch("Jet_toptag_tau1", tau1);
  SetBranch("Jet_toptag_tau2", tau2);
  SetBranch("Jet_toptag_tau3", tau3);
  SetBranch("Jet_toptag_PrunedMass", PrunedMass);
  SetBranch("Jet_toptag_SoftDropMass", SoftDropMass);
}

vector<CUTS> FatJet::findExtraCuts() {
  vector<CUTS> return_vec;
  if(pstats["Smear"].bset.find("SmearTheJet") != pstats["Smear"].bset.end()) {
    return_vec.push_back(CUTS::eGen);
  }
  return return_vec;
}

vector<CUTS> FatJet::overlapCuts(CUTS ePos) {
  vector<CUTS> returnCuts;
  unordered_set<string>& tmpset = pstats[jetNameMap.at(ePos)].bset;
  if(tmpset.find("RemoveOverlapWithMuon1s") != tmpset.end()) returnCuts.push_back(CUTS::eRMuon1);
  if(tmpset.find("RemoveOverlapWithMuon2s") != tmpset.end()) returnCuts.push_back(CUTS::eRMuon2);
  if(tmpset.find("RemoveOverlapWithElectron1s") != tmpset.end()) returnCuts.push_back(CUTS::eRElec1);
  if(tmpset.find("RemoveOverlapWithElectron2s") != tmpset.end()) returnCuts.push_back(CUTS::eRElec2);
  if(tmpset.find("RemoveOverlapWithTau1s") != tmpset.end()) returnCuts.push_back(CUTS::eRTau1);
  if(tmpset.find("RemoveOverlapWithTau2s") != tmpset.end()) returnCuts.push_back(CUTS::eRTau2);

  return returnCuts;
}


///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////    LEPTON   ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


Lepton::Lepton(TTree* _BOOM, string GenName, string EndName, vector<string> syst_names) : Particle(_BOOM, GenName, EndName, syst_names) {
  SetBranch((GenName+"_charge").c_str(), _charge);
}

vector<CUTS> Lepton::findExtraCuts() {
  vector<CUTS> return_vec;
  unordered_set<string>& tmpset = pstats["Smear"].bset;
  if(tmpset.find("SmearTheParticle") != tmpset.end() || tmpset.find("MatchToGen") != tmpset.end()) {
    return_vec.push_back(cutMap.at(type));
  }
  return return_vec;
}

double Lepton::charge(uint index)const     {return _charge->at(index);}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////    ELECTRON   ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


Electron::Electron(TTree* _BOOM, string filename, vector<string> syst_names) : Lepton(_BOOM, "patElectron", filename, syst_names) {
  type = PType::Electron;
  unordered_set<string>& elec1 = pstats["Elec1"].bset;
  unordered_set<string>& elec2 = pstats["Elec2"].bset;
  if(elec1.find("DoDiscrByIsolation") != elec1.end() || elec2.find("DoDiscrByIsolation") != elec2.end()) {
    SetBranch("patElectron_isoChargedHadrons", isoChargedHadrons);
    SetBranch("patElectron_isoNeutralHadrons", isoNeutralHadrons);
    SetBranch("patElectron_isoPhotons", isoPhotons);
    SetBranch("patElectron_isoPU", isoPU);
  }
  if(elec1.find("DoDiscrByVetoID") != elec1.end() || elec2.find("DoDiscrByVetoID") != elec2.end()) {
    SetBranch("patElectron_isPassVeto", isPassVeto);
  }
  if(elec1.find("DoDiscrByLooseID") != elec1.end() || elec2.find("DoDiscrByLooseID") != elec2.end()) {
    SetBranch("patElectron_isPassLoose", isPassLoose);
  }
  if(elec1.find("DoDiscrByMediumID") != elec1.end() || elec2.find("DoDiscrByMediumID") != elec2.end()) {
    SetBranch("patElectron_isPassMedium", isPassMedium);
  }
  if(elec1.find("DoDiscrByTightID") != elec1.end() || elec2.find("DoDiscrByTightID") != elec2.end()) {
    SetBranch("patElectron_isPassTight", isPassTight);
  }
  if(elec1.find("DoDiscrByHEEPID") != elec1.end() || elec2.find("DoDiscrByHEEPID") != elec2.end()) {
    SetBranch("patElectron_isPassHEEPId", isPassHEEPId);
  }
}


bool Electron::get_Iso(int index, double min, double max) const {
  double maxIsoval = std::max(0.0, isoNeutralHadrons->at(index) + isoPhotons->at(index) - 0.5 * isoPU->at(index));
  double isoSum = (isoChargedHadrons->at(index) + maxIsoval) / pt(index);
  return (isoSum >= min && isoSum < max);
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////   MUON  // ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


Muon::Muon(TTree* _BOOM, string filename, vector<string> syst_names) : Lepton(_BOOM, "Muon", filename, syst_names) {
  type = PType::Muon;
  unordered_set<string>& mu1 = pstats["Muon1"].bset;
  unordered_set<string>& mu2 = pstats["Muon2"].bset;

  if(mu1.find("DoDiscrByTightID") != mu1.end() || mu2.find("DoDiscrByTightID") != mu2.end()) {
    SetBranch("Muon_tight", tight);
     }
  if(mu1.find("DoDiscrBySoftID") != mu1.end() || mu2.find("DoDiscrBySoftID") != mu2.end()) {
    SetBranch("Muon_soft", soft);
  }
  if(mu1.find("DoDiscrByIsolation") != mu1.end() || mu2.find("DoDiscrByIsolation") != mu2.end()) {

    SetBranch("Muon_isoCharged", isoCharged);
    SetBranch("Muon_isoNeutralHadron", isoNeutralHadron);
    SetBranch("Muon_isoPhoton", isoPhoton);
    SetBranch("Muon_isoPU", isoPU);
  }
}

bool Muon::get_Iso(int index, double min, double max) const {
  double maxIsoval = std::max(0.0, isoNeutralHadron->at(index) + isoPhoton->at(index) - 0.5 * isoPU->at(index));
  double isoSum = (isoCharged->at(index) + maxIsoval) / pt(index);
  return (isoSum >= min && isoSum < max);
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////    TAUS    ///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


Taus::Taus(TTree* _BOOM, string filename, vector<string> syst_names) : Lepton(_BOOM, "Tau", filename, syst_names) {
  type = PType::Tau;

  ////Electron discrimination  
  string string_tau1 = pstats["Tau1"].smap["DiscrAgainstElectron"];
  string string_tau2 = pstats["Tau2"].smap["DiscrAgainstElectron"];

  SetBranch(("Tau_" + string_tau1).c_str(), againstElectron.first);
  if(string_tau1 != string_tau2) {
    SetBranch(("Tau_" + string_tau2).c_str(), againstElectron.second);
  } else {
    againstElectron.second = againstElectron.first;
  }

  ////Muon discrimination
  string_tau1 = pstats["Tau1"].smap["DiscrAgainstMuon"];
  string_tau2 = pstats["Tau2"].smap["DiscrAgainstMuon"];

  SetBranch(("Tau_" + string_tau1).c_str(), againstMuon.first);
  if(string_tau1 != string_tau2) {
    SetBranch(("Tau_" + string_tau2).c_str(), againstMuon.second);
  } else {
    againstMuon.second = againstMuon.first;
  }

  string_tau1 = pstats["Tau1"].smap["DiscrByMaxIsolation"];
  string_tau2 = pstats["Tau2"].smap["DiscrByMaxIsolation"];
  SetBranch(("Tau_" + string_tau1).c_str(), maxIso.first);
  if(string_tau1 != string_tau2) {
    SetBranch(("Tau_" + string_tau2).c_str(), maxIso.second);
  } else {
    maxIso.second = maxIso.first;
  }

  string_tau1 = pstats["Tau1"].smap["DiscrByMinIsolation"];
  string_tau2 = pstats["Tau2"].smap["DiscrByMinIsolation"];
  if(string_tau1 != "ZERO") {
    SetBranch(("Tau_" + string_tau1).c_str(), maxIso.first);
  }
  if(string_tau2 != "ZERO" && string_tau1 != string_tau2) {
    SetBranch(("Tau_"+ string_tau2).c_str(), maxIso.second);
  } else {
    maxIso.second = maxIso.first;
  }

  SetBranch("Tau_decayModeFindingNewDMs", decayModeFindingNewDMs);
  SetBranch("Tau_nProngs", nProngs);
  SetBranch("Tau_decayMode", decayMode);
  SetBranch("Tau_leadChargedCandPt", leadChargedCandPt);
  SetBranch("Tau_leadChargedCandTrack_ptError", leadChargedCandPtError);
  SetBranch("Tau_leadChargedCandValidHits", leadChargedCandValidHits);

}

vector<CUTS> Taus::findExtraCuts() {
  vector<CUTS> return_vec = Lepton::findExtraCuts();

  unordered_set<string>& tau1 = pstats["Tau1"].bset;
  unordered_set<string>& tau2 = pstats["Tau2"].bset;
  if(tau1.find("RemoveOverlapWithMuon1s") != tau1.end() || tau2.find("RemoveOverlapWithMuon1s") != tau2.end())
    return_vec.push_back(CUTS::eRMuon1);
  if(tau1.find("RemoveOverlapWithMuon2s") != tau1.end() || tau2.find("RemoveOverlapWithMuon2s") != tau2.end())
    return_vec.push_back(CUTS::eRMuon2);
  if(tau1.find("RemoveOverlapWithElectron1s") != tau1.end() || tau2.find("RemoveOverlapWithElectron1s") != tau2.end())
    return_vec.push_back(CUTS::eRElec1);
  if(tau1.find("RemoveOverlapWithElectron2s") != tau1.end() || tau2.find("RemoveOverlapWithElectron2s") != tau2.end())
    return_vec.push_back(CUTS::eRElec2);

  return return_vec;
}

bool Taus::get_Iso(int index, double onetwo, double max) const {
  double maxIsoval = (onetwo == 1) ? maxIso.first->at(index) : maxIso.second->at(index);
  vector<int>* minIsotmp = (onetwo == 1) ? minIso.first : minIso.second;
  double minIsoval = (minIsotmp != 0) ? minIsotmp->at(index) : true;
  return (maxIsoval > 0.5 && minIsoval > 0.5);
}

bool Taus::pass_against_Elec(CUTS ePos, int index) {
  return (ePos == CUTS::eRTau1) ? againstElectron.first->at(index) : againstElectron.second->at(index);
}

bool Taus::pass_against_Muon(CUTS ePos, int index) {
  return (ePos == CUTS::eRTau1) ? againstMuon.first->at(index) : againstMuon.second->at(index);
}
