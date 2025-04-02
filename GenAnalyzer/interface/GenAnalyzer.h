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
      
      // ----------member data ---------------------------

      //Tokens
      edm::Service<TFileService> fs;
      edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticlesToken_;
      edm::InputTag genParticles_;


       // Main TTree
      TTree *RHTree;




};
#endif
