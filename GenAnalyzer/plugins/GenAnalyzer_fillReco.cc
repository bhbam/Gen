#include "Gen/GenAnalyzer/interface/GenAnalyzer.h"
std::set<int> uniqueJetIdxs;
std::set<int> uniqueTauIdxs;

float V_A1_mass_jet_reco;
float V_A2_mass_jet_reco;
float V_H_mass_jet_reco;
float V_A1_mass_tau_reco;
float V_A2_mass_tau_reco;
float V_H_mass_tau_reco;
// TH1D *H_reco_mass;
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
  vJetIdxs.clear();
  vTauIdxs.clear();
  uniqueJetIdxs.clear();
  edm::Handle<std::vector<reco::GenParticle> > genParticles;
  iEvent.getByToken(genParticlesToken_,   genParticles);
  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);
  edm::Handle<reco::PFTauCollection> taus;
  iEvent.getByToken(tauCollectionT_, taus);
  bool pass_reco= false;
  for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen)
  {

    if ( abs(iGen->pdgId()) != 35 || iGen->numberOfDaughters() != 2 || iGen->daughter(0)->pdgId() != 25 || iGen->daughter(1)->pdgId() != 25 ) continue;
    if ( abs(iGen->daughter(0)->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(0)->daughter(1)->pdgId()) != 15 || abs(iGen->daughter(1)->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(1)->daughter(1)->pdgId()) != 15 ) continue;
    if ( abs(iGen->daughter(0)->daughter(0)->status()) != 2 || abs(iGen->daughter(0)->daughter(1)->status()) != 2 || abs(iGen->daughter(1)->daughter(0)->status()) != 2 || abs(iGen->daughter(1)->daughter(1)->status()) != 2 ) continue;

      for ( unsigned iJ(0); iJ != jets->size(); ++iJ )
      {
        reco::PFJetRef iJet( jets, iJ );
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

  if (pass_reco)
  {
    reco::PFJetRef iJet_0( jets, vJetIdxs[0] );
    reco::PFJetRef iJet_1( jets, vJetIdxs[1] );
    reco::PFTauRef iTau_0( taus, vTauIdxs[0] );
    reco::PFTauRef iTau_1( taus, vTauIdxs[1] );
    TLorentzVector A1_jet_inv  = SetTaus(iJet_0->pt(), iJet_0->eta(), iJet_0->phi(), iJet_0->mass());
    TLorentzVector A2_jet_inv  = SetTaus(iJet_1->pt(), iJet_1->eta(), iJet_1->phi(), iJet_1->mass());
    TLorentzVector H_jet_inv = A1_jet_inv + A2_jet_inv;
    TLorentzVector A1_tau_inv  = SetTaus(iTau_0->pt(), iTau_0->eta(), iTau_0->phi(), iTau_0->mass());
    TLorentzVector A2_tau_inv  = SetTaus(iTau_1->pt(), iTau_1->eta(), iTau_1->phi(), iTau_1->mass());
    TLorentzVector H_tau_inv = A1_tau_inv + A2_tau_inv;
    V_A1_mass_jet_reco = A1_jet_inv.M();
    V_A2_mass_jet_reco = A2_jet_inv.M();
    V_H_mass_jet_reco  = H_jet_inv.M();
    V_A1_mass_tau_reco = A1_tau_inv.M();
    V_A2_mass_tau_reco = A2_tau_inv.M();
    V_H_mass_tau_reco  = H_tau_inv.M();
    if (isDebug) std::cout << "A1_jet_inv: " << V_A1_mass_jet_reco << "  , A2_jet_inv: " << V_A2_mass_jet_reco<< " ,  H_jet_inv: "  << V_H_mass_jet_reco << std::endl;
    if (isDebug) std::cout << "A1_tau_inv: " << V_A1_mass_tau_reco << "  , A2_tau_inv: " << V_A2_mass_tau_reco<< " ,  H_tau_inv: "  << V_H_mass_tau_reco << std::endl;
    // H_reco_mass->Fill(V_H_mass_jet_reco);

  }

}
