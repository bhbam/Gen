#ifndef GenAnalyzer_h
#define GenAnalyzer_h
// -*- C++ -*-
//
// Package:    Gen/GenAnalyzer
// Class:      GenAnalyzer
//
/**\class GenAnalyzer GenAnalyzer.cc Gen/GenAnalyzer/plugins/GenAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Bhim Bam
//         Created:  Tue, 03 Sep 2024 18:24:42 GMT
//
//


// system include files
#include <memory>
#include <iostream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//for the GenParticleCollection and GenParticles
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Utilities/interface/RegexMatch.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include "DataFormats/METReco/interface/GenMET.h"


#include "TLorentzVector.h"
#include "TH2D.h"
#include "TTree.h"
//
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TVector3.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TLorentzVector.h"
#include <string>
#include <cstring>
#include <set>
using std::vector;



// class declaration

class GenAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
   public:
      explicit GenAnalyzer(const edm::ParameterSet&);
      ~GenAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //  flags
      bool isDebug;
      bool print_trigger;
      double minJetPt_;
      double maxJetEta_;
      // ----------member data ---------------------------

      //Tokens
      edm::Service<TFileService> fs;
      edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticlesToken_;
      edm::InputTag genParticles_;
      edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_ ;
      edm::EDGetTokenT<reco::PFJetCollection> jetCollectionT_;
      edm::EDGetTokenT<reco::PFTauCollection> tauCollectionT_;
      edm::EDGetTokenT<reco::PFMETCollection> metToken_;
      edm::EDGetTokenT<reco::PFMETCollection> pupimetToken_;
      edm::EDGetTokenT<std::vector<reco::GenMET>> genmetToken_;

       // Main TTree
      TTree *RHTree;


      // Selection and filling functions
      void branchesTrigger         ( TTree*, edm::Service<TFileService>& );
      void fillTrigger             ( const edm::Event&, const edm::EventSetup& );
      void branchesReco         ( TTree*, edm::Service<TFileService>& );
      void fillReco            ( const edm::Event&, const edm::EventSetup& );
      std::vector<int> vJetIdxs;
      std::vector<int> vTauIdxs;
      // function used in plugins
      TLorentzVector SetTaus(Float_t tau_pt, Float_t tau_eta, Float_t tau_phi, Float_t tau_mass){
        TLorentzVector Tau_Candidate;
        Tau_Candidate.SetPtEtaPhiM(tau_pt, tau_eta, tau_phi, tau_mass);
        return Tau_Candidate;
      }


};
#endif
