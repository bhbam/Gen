#include "Gen/GenAnalyzer/interface/GenAnalyzer.h"
std::set<int> uniqueJetIdxs;
std::set<int> uniqueTauIdxs;

float V_A1_mass_jet_reco;
float V_A2_mass_jet_reco;
float V_H_mass_jet_reco;
float V_A1_mass_tau_reco;
float V_A2_mass_tau_reco;
float V_H_mass_tau_reco;
float V_dR_jets_reco;
float V_dR_taus_reco;
float V_N_matched_jets_reco;
float V_N_matched_taus_reco;

float V_met_e_reco;
float V_met_pt_reco;
float V_met_px_reco;
float V_met_py_reco;
float V_met_phi_reco;
float V_met_significance_reco;

float V_pupimet_e_reco;
float V_pupimet_pt_reco;
float V_pupimet_px_reco;
float V_pupimet_py_reco;
float V_pupimet_phi_reco;
float V_pupimet_significance_reco;

float V_genmet_e_reco;
float V_genmet_pt_reco;
float V_genmet_px_reco;
float V_genmet_py_reco;
float V_genmet_phi_reco;
float V_genmet_significance_reco;
//-----------------------now do what ever initialization is needed

void GenAnalyzer::branchesReco(TTree* tree, edm::Service<TFileService> &fs)
{
  // H_reco_mass     = fs->make<TH1D>("H_mass_jet_reco"   , "H_mass_jet_reco;Events"                 ,  30,  160, 100);
  tree->Branch("A1_mass_jet_reco",  &V_A1_mass_jet_reco);
  tree->Branch("A2_mass_jet_reco",  &V_A2_mass_jet_reco);
  tree->Branch("H_mass_jet_reco",   &V_H_mass_jet_reco);
  tree->Branch("A1_mass_tau_reco",  &V_A1_mass_tau_reco);
  tree->Branch("A2_mass_tau_reco",  &V_A2_mass_tau_reco);
  tree->Branch("H_mass_tau_reco",   &V_H_mass_tau_reco);
  tree->Branch("dR_jets_reco",      &V_dR_jets_reco);
  tree->Branch("dR_taus_reco",      &V_dR_taus_reco);
  tree->Branch("N_matched_jets_reco",      &V_N_matched_jets_reco);
  tree->Branch("N_matched_taus_reco",      &V_N_matched_taus_reco);

  tree->Branch("met_e_reco",      &V_met_e_reco);
  tree->Branch("met_pt_reco",      &V_met_pt_reco);
  tree->Branch("met_px_reco",      &V_met_px_reco);
  tree->Branch("met_py_reco",      &V_met_py_reco);
  tree->Branch("met_phi_reco",      &V_met_phi_reco);
  tree->Branch("met_significance_reco",      &V_met_significance_reco);

  tree->Branch("pupimet_e_reco",      &V_pupimet_e_reco);
  tree->Branch("pupimet_pt_reco",      &V_pupimet_pt_reco);
  tree->Branch("pupimet_px_reco",      &V_pupimet_px_reco);
  tree->Branch("pupimet_py_reco",      &V_pupimet_py_reco);
  tree->Branch("pupimet_phi_reco",      &V_pupimet_phi_reco);
  tree->Branch("pupimet_significance_reco",      &V_pupimet_significance_reco);

  tree->Branch("genmet_e_reco",      &V_genmet_e_reco);
  tree->Branch("genmet_pt_reco",      &V_genmet_pt_reco);
  tree->Branch("genmet_px_reco",      &V_genmet_px_reco);
  tree->Branch("genmet_py_reco",      &V_genmet_py_reco);
  tree->Branch("genmet_phi_reco",      &V_genmet_phi_reco);
  tree->Branch("genmet_significance_reco",      &V_genmet_significance_reco);



}

// ---------------------- Fill tree with reco info  ------------
void GenAnalyzer::fillReco(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  V_A1_mass_jet_reco = -1111.1111;
  V_A2_mass_jet_reco = -1111.1111;
  V_H_mass_jet_reco  = -1111.1111;
  V_A1_mass_tau_reco = -1111.1111;
  V_A2_mass_tau_reco = -1111.1111;
  V_H_mass_tau_reco  = -1111.1111;
  V_dR_jets_reco  = -1111.1111;
  V_dR_taus_reco  = -1111.1111;
  V_N_matched_jets_reco  = -1111.1111;
  V_N_matched_taus_reco  = -1111.1111;

  V_met_e_reco  = -1111.1111;
  V_met_pt_reco  = -1111.1111;
  V_met_px_reco  = -1111.1111;
  V_met_py_reco  = -1111.1111;
  V_met_phi_reco  = -1111.1111;
  V_met_significance_reco  = -1111.1111;

  V_pupimet_e_reco  = -1111.1111;
  V_pupimet_pt_reco  = -1111.1111;
  V_pupimet_px_reco  = -1111.1111;
  V_pupimet_py_reco  = -1111.1111;
  V_pupimet_phi_reco  = -1111.1111;
  V_pupimet_significance_reco  = -1111.1111;

  V_genmet_e_reco  = -1111.1111;
  V_genmet_pt_reco  = -1111.1111;
  V_genmet_px_reco  = -1111.1111;
  V_genmet_py_reco  = -1111.1111;
  V_genmet_phi_reco  = -1111.1111;
  V_genmet_significance_reco  = -1111.1111;

  vJetIdxs.clear();
  vTauIdxs.clear();
  uniqueJetIdxs.clear();

  edm::Handle<std::vector<reco::GenParticle> > genParticles;
  iEvent.getByToken(genParticlesToken_,   genParticles);

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  edm::Handle<reco::PFTauCollection> taus;
  iEvent.getByToken(tauCollectionT_, taus);

  edm::Handle<reco::PFMETCollection> metHandle;
  iEvent.getByToken(metToken_, metHandle);

  edm::Handle<reco::PFMETCollection> pupimetHandle;
  iEvent.getByToken(pupimetToken_, pupimetHandle);

  edm::Handle<std::vector<reco::GenMET>> genmetHandle;
  iEvent.getByToken(genmetToken_, genmetHandle);

  if (metHandle.isValid() && !metHandle->empty())
    {

      float met_e = metHandle->begin()->sumEt();
      float met_pt = metHandle->begin()->pt();
      float met_px = metHandle->begin()->px();
      float met_py = metHandle->begin()->py();
      float met_phi = metHandle->begin()->phi();
      float met_significance = metHandle->begin()->significance();
      V_met_e_reco  = met_e;
      V_met_pt_reco  = met_pt;
      V_met_px_reco  = met_px;
      V_met_py_reco  = met_py;
      V_met_phi_reco  = met_phi;
      V_met_significance_reco  = met_significance;
      if (isDebug) std::cout<< "MET: pt = " << met_pt<< " , phi = " << met_phi<< ", E = " << met_e<< ", px  "<< met_px<<std::endl;
    }

  else
    {
     std::cout<<"LogWarning(METAnalyzer)" << "MET collection not found!"<<std::endl;
    }


    if (pupimetHandle.isValid() && !pupimetHandle->empty())
      {

        float pupimet_e = pupimetHandle->begin()->sumEt();
        float pupimet_pt = pupimetHandle->begin()->pt();
        float pupimet_px = pupimetHandle->begin()->px();
        float pupimet_py = pupimetHandle->begin()->py();
        float pupimet_phi = pupimetHandle->begin()->phi();
        float pupimet_significance = pupimetHandle->begin()->significance();
        V_pupimet_e_reco  = pupimet_e;
        V_pupimet_pt_reco  = pupimet_pt;
        V_pupimet_px_reco  = pupimet_px;
        V_pupimet_py_reco  = pupimet_py;
        V_pupimet_phi_reco  = pupimet_phi;
        V_pupimet_significance_reco  = pupimet_significance;
        if (isDebug) std::cout<< "PUPI MET: pt = " << pupimet_pt<< " , phi = " << pupimet_phi<< ", E = " << pupimet_e<< ", px  "<< pupimet_px<<std::endl;
      }

      else
        {
         std::cout<<"LogWarning(METAnalyzer)" << "PUPI MET collection not found!"<<std::endl;
        }


  float genmet_e = (genmetHandle->front()).sumEt();
  float genmet_pt = (genmetHandle->front()).pt();
  float genmet_px = (genmetHandle->front()).px();
  float genmet_py = (genmetHandle->front()).py();
  float genmet_phi = (genmetHandle->front()).phi();
  float genmet_significance = (genmetHandle->front()).significance();
  V_genmet_e_reco  = genmet_e;
  V_genmet_pt_reco  = genmet_pt;
  V_genmet_px_reco  = genmet_px;
  V_genmet_py_reco  = genmet_py;
  V_genmet_phi_reco  = genmet_phi;
  V_genmet_significance_reco  = genmet_significance;
  if (isDebug) std::cout<< "GEN MET: pt = " << genmet_pt<< " , phi = " << genmet_phi<< ", E = " << genmet_e<< ", px  "<< genmet_px<<std::endl;


  bool pass_reco= false;
  for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen)
  {

    if ( abs(iGen->pdgId()) != 35 || iGen->numberOfDaughters() != 2 || iGen->daughter(0)->pdgId() != 25 || iGen->daughter(1)->pdgId() != 25 ) continue;
    if ( abs(iGen->daughter(0)->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(0)->daughter(1)->pdgId()) != 15 || abs(iGen->daughter(1)->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(1)->daughter(1)->pdgId()) != 15 ) continue;
    if ( abs(iGen->daughter(0)->daughter(0)->status()) != 2 || abs(iGen->daughter(0)->daughter(1)->status()) != 2 || abs(iGen->daughter(1)->daughter(0)->status()) != 2 || abs(iGen->daughter(1)->daughter(1)->status()) != 2 ) continue;

      for ( unsigned iJ(0); iJ != jets->size(); ++iJ )
      {
        reco::PFJetRef iJet( jets, iJ );
        if ( std::abs(iJet->pt())  < minJetPt_ ) continue;
        if ( std::abs(iJet->eta()) > maxJetEta_ ) continue;
        float dR_jet_A1 = reco::deltaR( iJet->eta(),iJet->phi(), iGen->daughter(0)->eta(),iGen->daughter(0)->phi() );
        float dR_jet_A2 = reco::deltaR( iJet->eta(),iJet->phi(), iGen->daughter(1)->eta(),iGen->daughter(1)->phi() );

        if (!((dR_jet_A1 < 0.4) ^ (dR_jet_A2 < 0.4))) continue;
        if (isDebug) std::cout << "dR_jet_A1: " << dR_jet_A1 << ", dR_jet_A2: " << dR_jet_A2 << std::endl;
        if (isDebug) std::cout << "dR_jet_passed: " << std::endl;
        uniqueJetIdxs.insert( iJ );
      }
      vJetIdxs.assign(uniqueJetIdxs.begin(), uniqueJetIdxs.end());
      if (isDebug) std::cout << "   Number of Selected Jet  "<< vJetIdxs.size() << std::endl;

      if (vJetIdxs.size()>1)
      {
          for ( unsigned iJ(0); iJ != vJetIdxs.size(); ++iJ )
          {
            reco::PFJetRef iJet( jets, vJetIdxs[iJ] );
            for ( unsigned iT(0); iT != taus->size(); ++iT )
            {
              reco::PFTauRef iTau( taus, iT );
              float dR_jet_tau = reco::deltaR( iJet->eta(),iJet->phi(), iTau->eta(),iTau->phi());
              if (dR_jet_tau > 0.4) continue;
              if (isDebug) std::cout << "dR_jet_tau: " << dR_jet_tau<<std::endl;
              uniqueTauIdxs.insert( iT );
            }
          }
          if (isDebug) std::cout << "   Number of Selected Tau  "<<uniqueTauIdxs.size() << std::endl;
          vTauIdxs.assign(uniqueTauIdxs.begin(), uniqueTauIdxs.end());
          if (vTauIdxs.size() > 1)
          {
            pass_reco = true;
          }

        }

  }
  // Taking first two Jets and Tau
  reco::PFJetRef iJet_0( jets, vJetIdxs[0] );
  reco::PFJetRef iJet_1( jets, vJetIdxs[1] );
  reco::PFTauRef iTau_0( taus, vTauIdxs[0] );
  reco::PFTauRef iTau_1( taus, vTauIdxs[1] );
  float dR_reco_jets = reco::deltaR( iJet_0->eta(),iJet_0->phi(), iJet_1->eta(),iJet_1->phi() );
  float dR_reco_taus = reco::deltaR( iTau_0->eta(),iTau_0->phi(), iTau_1->eta(),iTau_1->phi() );
  TLorentzVector A1_jet_inv  = SetTaus(iJet_0->pt(), iJet_0->eta(), iJet_0->phi(), iJet_0->mass());
  TLorentzVector A2_jet_inv  = SetTaus(iJet_1->pt(), iJet_1->eta(), iJet_1->phi(), iJet_1->mass());
  TLorentzVector H_jet_inv = A1_jet_inv + A2_jet_inv;
  TLorentzVector A1_tau_inv  = SetTaus(iTau_0->pt(), iTau_0->eta(), iTau_0->phi(), iTau_0->mass());
  TLorentzVector A2_tau_inv  = SetTaus(iTau_1->pt(), iTau_1->eta(), iTau_1->phi(), iTau_1->mass());
  TLorentzVector H_tau_inv = A1_tau_inv + A2_tau_inv;
  if (pass_reco && (dR_reco_jets > 0.5) && (dR_reco_taus >0.5) )
  {

    V_A1_mass_jet_reco = A1_jet_inv.M();
    V_A2_mass_jet_reco = A2_jet_inv.M();
    V_H_mass_jet_reco  = H_jet_inv.M();
    V_A1_mass_tau_reco = A1_tau_inv.M();
    V_A2_mass_tau_reco = A2_tau_inv.M();
    V_H_mass_tau_reco  = H_tau_inv.M();
    V_dR_jets_reco     = dR_reco_jets;
    V_dR_taus_reco     = dR_reco_taus;
    V_N_matched_jets_reco     = vJetIdxs.size();
    V_N_matched_taus_reco     = vTauIdxs.size();

    if (isDebug) std::cout << "A1_jet_inv: " << V_A1_mass_jet_reco << "  , A2_jet_inv: " << V_A2_mass_jet_reco<< " ,  H_jet_inv: "  << V_H_mass_jet_reco << std::endl;
    if (isDebug) std::cout << "A1_tau_inv: " << V_A1_mass_tau_reco << "  , A2_tau_inv: " << V_A2_mass_tau_reco<< " ,  H_tau_inv: "  << V_H_mass_tau_reco << std::endl;

  }

}
