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
// #include <iostream>
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
//
#include "Gen/GenAnalyzer/inference/unbaising_function.h"

// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.



int pdgid_ = 25; // A/H to tautau
// int pdgid_ = 553; // Upsilon 1S to tauatau


using std::vector;
using reco::GenParticle;

int ntotal_event ;
int npassed_event ;

unsigned int runId_;
unsigned int lumiId_;
unsigned long long eventId_;

TH1D *H_tau_att_genA1_M_inv;
TH1D *H_tau_att_genA1_M;
TH1D *H_tau_att_dR_A1_Tau1;
TH1D *H_tau_att_dR_A1_Tau2;
TH1D *H_tau_att_dR_Tau1_Tau2;

TH1D *H_tau_att_A1_pt;
TH1D *H_tau_att_Tau1_pt;
TH1D *H_tau_att_Tau2_pt;
TH1D *H_tau_att_A1_eta;

TH1D *H_tau_att_Tau1_eta;
TH1D *H_tau_att_Tau2_eta;
TH1D *H_tau_att_A1_phi;
TH1D *H_tau_att_Tau1_phi;
TH1D *H_tau_att_Tau2_phi;

TH1D *H_tau_att_Tau1_Tau2_deta;
TH1D *H_tau_att_Tau1_Tau2_dphi;

TH2D *H_tau_att_Tau1_Tau2_dphi_deta;

float V_att_genA1_M_inv;
float V_att_genA1_M;
float V_att_dR_A1_Tau1;
float V_att_dR_A1_Tau2;
float V_att_dR_Tau1_Tau2;

float V_att_A1_pt;
float V_att_Tau1_pt;
float V_att_Tau2_pt;
float V_att_A1_eta;
float V_att_Tau1_eta;
float V_att_Tau2_eta;
float V_att_A1_phi;
float V_att_Tau1_phi;
float V_att_Tau2_phi;

float V_att_Tau1_Tau2_deta;
float V_att_Tau1_Tau2_dphi;

// vector<int> vAIdxs;
vector<float> V_att_genA1_M_inv_;
vector<float> V_att_genA1_M_;
vector<float> V_att_dR_A1_Tau1_;
vector<float> V_att_dR_A1_Tau2_;
vector<float> V_att_dR_Tau1_Tau2_;

vector<float> V_att_A1_pt_;
vector<float> V_att_Tau1_pt_;
vector<float> V_att_Tau2_pt_;
vector<float> V_att_A1_eta_;
vector<float> V_att_Tau1_eta_;
vector<float> V_att_Tau2_eta_;
vector<float> V_att_A1_phi_;
vector<float> V_att_Tau1_phi_;
vector<float> V_att_Tau2_phi_;

vector<float> V_att_Tau1_Tau2_deta_;
vector<float> V_att_Tau1_Tau2_dphi_;


TLorentzVector SetTaus(Float_t tau_pt, Float_t tau_eta, Float_t tau_phi, Float_t tau_mass){
  TLorentzVector Tau_Candidate;
  Tau_Candidate.SetPtEtaPhiM(tau_pt, tau_eta, tau_phi, tau_mass);
  return Tau_Candidate;
}

// TLorentzVector SetMETs(Float_t met, Float_t metphi){
//   TLorentzVector Met;
//   Met.SetPtEtaPhiM(met, 0, metphi, 0.);
//   return Met;
// }

class GenAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit GenAnalyzer(const edm::ParameterSet&);
      ~GenAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      // edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      bool unbiasing;
      edm::Service<TFileService> fs;
      edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticlesToken_;
      edm::InputTag genParticles_;

      TTree *RHTree;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GenAnalyzer::GenAnalyzer(const edm::ParameterSet& iConfig)
 // :
  // tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))

{
   //now do what ever initialization is needed
   RHTree = fs->make<TTree>("RHTree","Gen info Tree");

   H_tau_att_genA1_M_inv     = fs->make<TH1D>("h_genA1_M_inv"   , "m^{gen_inv A1};m^{gen_inv A1};Events"                     ,  42,  3.6, 20);
   H_tau_att_genA1_M     = fs->make<TH1D>("h_genA1_M"   , "m^{gen A1};m^{gen A1};Events"                     ,  42,  3.6, 20);
   H_tau_att_dR_A1_Tau1     = fs->make<TH1D>("h_dR_A1_Tau1"   , "dR^{gen A1_Tau1};dR^{gen A1_Tau1};Events"                     ,  10,  0, 5);
   H_tau_att_dR_A1_Tau2     = fs->make<TH1D>("h_dR_A1_Tau2"   , "dR^{gen A1_Tau2};dR^{gen A1_Tau2};Events"                     ,  10,  0, 5);
   H_tau_att_dR_Tau1_Tau2     = fs->make<TH1D>("h_dR_Tau1_Tau2"   , "dR^{gen Tau1_Tau2};dR^{gen Tau1_Tau2};Events"                     ,  10,  0, 2);
   H_tau_att_A1_pt     = fs->make<TH1D>("h_A1_pt"   , "pt^{gen A1};pt^{gen A1};Events"                     ,  35,  30, 200);
   H_tau_att_Tau1_pt     = fs->make<TH1D>("h_Tau1_pt"   , "pt^{gen Tau1};pt^{gen Tau1};Events"                     ,  70,  10, 150);
   H_tau_att_Tau2_pt     = fs->make<TH1D>("h_Tau2_pt"   , "pt^{gen Tau2};pt^{gen Tau2};Events"                     ,  70,  10, 150);
   H_tau_att_A1_eta     = fs->make<TH1D>("h_A1_eta"   , "eta^{gen A1};eta^{gen A1};Events"                     ,  20,  -6, 6);
   H_tau_att_Tau1_eta     = fs->make<TH1D>("h_Tau1_eta"   , "eta^{gen Tau1};eta^{gen Tau1};Events"                     ,  20,  -6, 6);
   H_tau_att_Tau2_eta     = fs->make<TH1D>("h_Tau2_eta"   , "eta^{gen Tau2};eta^{gen Tau2};Events"                     ,  20,  -6, 6);
   H_tau_att_A1_phi     = fs->make<TH1D>("h_A1_phi"   , "phi^{gen A1};phi^{gen A1};Events"                     ,  20,  -3.2, 3.2);
   H_tau_att_Tau1_phi     = fs->make<TH1D>("h_Tau1_phi"   , "phi^{gen Tau1};phi^{gen Tau1};Events"                     ,  20,  -3.2, 3.2);
   H_tau_att_Tau2_phi     = fs->make<TH1D>("h_Tau2_phi"   , "phi^{gen Tau2};phi^{gen Tau2};Events"                     ,  20,  -3.2, 3.2);

   H_tau_att_Tau1_Tau2_deta     = fs->make<TH1D>("h_Tau1_Tau2_deta"   , "deta^{gen Tau1_Tau2};deta^{gen Tau1_Tau2};Events"                     ,  20,  0, .8);
   H_tau_att_Tau1_Tau2_dphi     = fs->make<TH1D>("h_Tau1_Tau2_dphi"   , "dphi^{gen Tau1_Tau2};dphi^{gen Tau1_Tau2};Events"                     ,  20,  0, .8);

   H_tau_att_Tau1_Tau2_dphi_deta     = fs->make<TH2D>("h_Tau1_Tau2_dphi_deta"   , "dphi vs deta^{gen Tau1_Tau2};dphi;deta" , 20,  0, .5 , 20,  0, .5);

   RHTree->Branch("Event",  &eventId_);
   RHTree->Branch("Run",  &runId_);
   RHTree->Branch("LumiSection",  &lumiId_);

   RHTree->Branch("GenA1_inv",  &V_att_genA1_M_inv);
   RHTree->Branch("GenA1",  &V_att_genA1_M);
   RHTree->Branch("dR_A1_Tau1",  &V_att_dR_A1_Tau1);
   RHTree->Branch("dR_A1_Tau2",  &V_att_dR_A1_Tau2);
   RHTree->Branch("dR_Tau1_Tau2",  &V_att_dR_Tau1_Tau2);

   RHTree->Branch("A1_pt",  &V_att_A1_pt);
   RHTree->Branch("Tau1_pt",  &V_att_Tau1_pt);
   RHTree->Branch("Tau2_pt",  &V_att_Tau2_pt);
   RHTree->Branch("A1_eta",  &V_att_A1_eta);
   RHTree->Branch("Tau1_eta",  &V_att_Tau1_eta);
   RHTree->Branch("Tau2_eta",  &V_att_Tau2_eta);
   RHTree->Branch("A1_phi",  &V_att_A1_phi);
   RHTree->Branch("Tau1_phi",  &V_att_Tau1_phi);
   RHTree->Branch("Tau2_phi",  &V_att_Tau2_phi);

   RHTree->Branch("Tau1_Tau2_deta",  &V_att_Tau1_Tau2_deta);
   RHTree->Branch("Tau1_Tau2_dphi",  &V_att_Tau1_Tau2_dphi);

   genParticlesToken_   = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"));
   unbiasing  = iConfig.getParameter<bool>("unbiasing");
   std::cout<<" unbiasing----> "<< unbiasing <<std::endl;
}


GenAnalyzer::~GenAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  eventId_ = iEvent.id().event();
  runId_ = iEvent.id().run();
  lumiId_ = iEvent.id().luminosityBlock();
  // std::cout<<"Event info"<<eventId_<<":  "<<runId_<<":  "<<lumiId_<<std::endl;

   V_att_genA1_M_inv    = -1111.1111;
   V_att_genA1_M    = -1111.1111;
   V_att_dR_A1_Tau1    = -1111.1111;
   V_att_dR_A1_Tau2    = -1111.1111;
   V_att_dR_Tau1_Tau2    = -1111.1111;
   V_att_A1_pt    = -1111.1111;
   V_att_Tau1_pt    = -1111.1111;
   V_att_Tau2_pt    = -1111.1111;
   V_att_A1_eta    = -1111.1111;
   V_att_Tau1_eta    = -1111.1111;
   V_att_Tau2_eta    = -1111.1111;
   V_att_A1_phi    = -1111.1111;
   V_att_Tau1_phi    = -1111.1111;
   V_att_Tau2_phi    = -1111.1111;

   V_att_Tau1_Tau2_deta    = -1111.1111;
   V_att_Tau1_Tau2_dphi    = -1111.1111;

   edm::Handle<std::vector<reco::GenParticle> > genParticles;
   iEvent.getByToken(genParticlesToken_,   genParticles);



  // unsigned int NAs = 0;
  // // unsigned int NTau_fromA = 0;
  // // for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
  // //   if ( abs(iGen->pdgId()) != 15 || abs(iGen->mother()->pdgId()) != 25) continue;
  // //   NTau_fromA++;
  // // }
  // // std::cout << "  >>>>>> Number Tau from  A <<<<<"<<"    "  <<  NTau_fromA << std::endl;
  // // vAIdxs.clear();
  // for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
  //   if ( iGen->pdgId() != -25 || abs(iGen->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(1)->pdgId()) != 15) continue;
  //   NAs++;
  //
  //   // vAIdxs.push_back(NAs-1);
  //   // std::cout<<"Size------------"<<vAIdxs.size()<<std::endl;
  // }
  // std::cout<<" Anti particle "<<NAs<<std::endl;
  // std::cout << "  >>>>>> Number of A giving Tau <<<<<"<<"    "  <<  NAs << std::endl;

float genA1_mass_inv = -1111.1111;
float genA1_mass = -1111.1111;
float A1_Tau1_dR = -1111.1111;
float A1_Tau2_dR = -1111.1111;
float Tau1_Tau2_dR = -1111.1111;
float A1_pt = -1111.1111;
float Tau1_pt = -1111.1111;
float Tau2_pt = -1111.1111;
float A1_eta = -1111.1111;
float Tau1_eta = -1111.1111;
float Tau2_eta = -1111.1111;
float A1_phi = -1111.1111;
float Tau1_phi = -1111.1111;
float Tau2_phi = -1111.1111;

float Tau1_Tau2_deta = -1111.1111;
float Tau1_Tau2_dphi = -1111.1111;



for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) { //Gen loop
  bool pass_gen = false;
  bool pass_unbaising = true;
  if (abs(iGen->pdgId()) != pdgid_ || iGen->numberOfDaughters() != 2 || abs(iGen->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(1)->pdgId()) != 15 ) continue;
  ntotal_event++;
  if ( abs(iGen->daughter(0)->status()) != 2 || abs(iGen->daughter(1)->status()) != 2 || iGen->daughter(0)->numberOfMothers() < 1 || iGen->daughter(1)->numberOfMothers() < 1) continue;


  float dR_A1_Tau1 = reco::deltaR( iGen->daughter(0)->eta(), iGen->daughter(0)->phi(), iGen->eta(), iGen->phi());
  float dR_A1_Tau2 = reco::deltaR( iGen->daughter(1)->eta(), iGen->daughter(1)->phi(), iGen->eta(), iGen->phi());
  float dR_Tau1_Tau2 = reco::deltaR( iGen->daughter(0)->eta(), iGen->daughter(0)->phi(), iGen->daughter(1)->eta(), iGen->daughter(1)->phi());

  // if(dR_Tau1_Tau2 > 0.4) continue;
  pass_gen = true;

  if (unbiasing)
  {

    std::vector <int> pT_bins   = {35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300};
    std::vector <double> m_bins = {4.0, 4.4, 4.8, 5.2, 5.6, 6.0, 6.4, 6.8, 7.2, 7.6, 8.0, 8.4, 8.8, 9.2, 9.6, 10.0, 10.4, 10.8, 11.2, 11.6, 12.0, 12.4, 12.8, 13.2, 13.6, 14.0, 14.4, 14.8, 15.2, 15.6, 16.0, 16.4, 16.8, 17.2, 17.6, 18};
    std::vector <int> occ =
    {/*35    40    45    50   55    60     65    70    75    80    85    90    95    100   105   110  115    120   125   130   135   140   145   150   155   160   165  170    175   180   185   190   195  200    205    210  215   220   225   230   235   240   245   250   255   260  265    270   275   280  285    290   295   300 */
      8756, 8931, 8730, 8679, 8789, 8674, 8647, 8708, 8829, 8688, 9185, 8729, 8870, 8608, 8551, 8674, 8848, 8805, 8782, 8764, 8713, 8745, 8783, 8699, 8819, 8773, 8742, 8752, 8778, 8706, 8801, 9036, 8937, 8827, 8764, 8807, 8702, 8994, 8736, 8913, 9018, 8722, 8965, 8838, 9024, 8684, 8897, 8713, 8700, 8857, 8829, 8924, 8825, 8699,
      10100, 9858, 10010, 10379, 10171, 10122, 10323, 10158, 10019, 10337, 10143, 10015, 9940, 10131, 10054, 10121, 9862, 9863, 10175, 10237, 10169, 10077, 10249, 10038, 10057, 10096, 10116, 9971, 10141, 10226, 10078, 10092, 10209, 10190, 10111, 9868, 9880, 9984, 10131, 10203, 9951, 9905, 10136, 10230, 10098, 9975, 10076, 10160, 10115, 10092, 10350, 9978, 9903, 10115,
      9891, 9937, 9801, 10044, 10093, 9796, 9820, 9892, 9832, 9843, 9819, 9936, 9631, 10031, 10024, 9832, 10216, 10292, 9782, 10231, 9938, 9738, 10111, 10050, 9776, 9780, 9872, 9814, 9928, 10024, 10231, 9967, 9879, 9884, 10135, 10040, 9788, 9988, 10024, 10083, 9975, 9788, 9901, 9774, 9593, 9917, 10073, 9917, 10008, 9850, 9769, 9733, 10064, 9916,
      9803, 9802, 9807, 9738, 10072, 9752, 9727, 9517, 9890, 9944, 9858, 9704, 9452, 9827, 10039, 9937, 9671, 9870, 9759, 9882, 9881, 9926, 9719, 10015, 9489, 9750, 9662, 9839, 9835, 9901, 9842, 9763, 9968, 9561, 9624, 9790, 10084, 9816, 10049, 9947, 9958, 9706, 9957, 9745, 9822, 9803, 9824, 9998, 9586, 9841, 9560, 9777, 9507, 10143,
      9763, 9610, 9613, 9592, 9966, 9729, 9550, 9750, 9688, 9716, 9423, 9956, 10104, 9867, 9726, 9585, 9703, 9701, 9793, 9799, 9713, 9729, 9631, 9640, 9398, 9456, 9943, 10014, 10005, 9571, 9755, 9697, 9651, 9567, 9674, 9749, 9525, 9644, 9564, 9715, 9498, 9759, 9571, 9797, 9607, 9948, 9600, 9810, 9611, 9695, 9821, 9763, 9825, 9773,
      9447, 9498, 9724, 9376, 9505, 9804, 9956, 9509, 9623, 9660, 9838, 9534, 9544, 9753, 9597, 9441, 9684, 9649, 9905, 9703, 9647, 9716, 9442, 9628, 9817, 9577, 9528, 9626, 9621, 9356, 9705, 9627, 9575, 9739, 9653, 9813, 9747, 9717, 9615, 9759, 9731, 9655, 9331, 9625, 9585, 9573, 9707, 9790, 9591, 9429, 9612, 9557, 9584, 9469,
      9578, 9659, 9631, 9378, 9736, 9556, 9598, 9509, 9671, 9747, 9606, 9571, 9418, 9174, 9607, 9360, 9532, 9417, 9417, 9535, 9537, 9482, 9503, 9297, 9464, 9654, 9724, 9488, 9639, 9484, 9479, 9643, 9816, 9605, 9668, 9698, 9494, 9421, 9495, 9474, 9423, 9501, 9535, 9566, 9738, 9458, 9604, 9539, 9712, 9481, 9396, 9527, 9514, 9507,
      9408, 9186, 9498, 9249, 9475, 9492, 9445, 9289, 9385, 9602, 9485, 9336, 9368, 9479, 9360, 9353, 9381, 9323, 9416, 9302, 9519, 9619, 9343, 9477, 9545, 9534, 9358, 9494, 9688, 9535, 9348, 9251, 9419, 9327, 9459, 9389, 9325, 9238, 9440, 9326, 9469, 9268, 9450, 9484, 9280, 9539, 9743, 9464, 9387, 9440, 9259, 9393, 9360, 9241,
      9294, 9535, 9506, 9297, 9440, 9313, 9328, 9380, 9093, 9271, 9389, 9323, 9518, 9251, 9357, 9439, 9317, 9348, 9114, 9330, 9259, 9111, 9115, 9391, 9427, 9509, 9429, 9415, 9301, 9340, 9324, 9278, 9280, 9160, 9350, 9176, 9347, 9227, 9292, 9539, 9484, 9366, 9367, 9503, 9349, 9253, 9328, 9384, 9349, 9012, 9419, 9280, 9503, 9457,
      9222, 9175, 9196, 9359, 9078, 9167, 9186, 9098, 9447, 9370, 9527, 9245, 9278, 9280, 9247, 9094, 9332, 9364, 9334, 9298, 9437, 9371, 9132, 9463, 9098, 9181, 9477, 9131, 8977, 9251, 9246, 9262, 9394, 9352, 9286, 9363, 9180, 9370, 9484, 9271, 9333, 9301, 9168, 9308, 9224, 9547, 9345, 9171, 9454, 9228, 8998, 9206, 9365, 9191,
      9091, 9093, 9346, 9222, 9152, 9302, 9046, 9263, 9344, 9408, 9238, 8800, 9297, 9437, 9112, 9083, 9092, 9094, 9273, 9164, 9320, 9319, 9181, 9214, 9537, 9313, 9016, 9020, 9199, 9091, 9548, 9156, 9223, 9205, 9199, 8930, 9206, 9259, 9426, 9341, 9157, 9139, 9178, 9140, 9374, 9032, 9395, 9221, 9255, 9311, 9231, 9222, 9168, 9148,
      9173, 8925, 9217, 9191, 9267, 9365, 9075, 9237, 9187, 8901, 9085, 8917, 9309, 9135, 9080, 9084, 9280, 8954, 9162, 9088, 9060, 9028, 9109, 9193, 9126, 9338, 8864, 9164, 9223, 9144, 9115, 9128, 9260, 9082, 9330, 9082, 9507, 9167, 9127, 8917, 9350, 8838, 9001, 9309, 8992, 9295, 9067, 9231, 9122, 9081, 9081, 9087, 8979, 9107,
      9217, 8731, 9207, 9108, 9215, 9096, 9010, 9107, 9045, 9280, 9004, 9205, 9135, 9254, 8988, 9258, 9083, 8930, 9039, 9086, 9325, 8962, 9146, 9261, 8990, 8892, 8881, 9085, 9308, 8934, 9162, 9112, 9241, 9317, 9184, 8974, 9060, 9164, 9063, 8831, 8991, 9112, 9123, 9028, 9158, 9066, 8997, 9136, 8999, 8950, 9148, 8971, 9246, 8949,
      9151, 9126, 9101, 8964, 8946, 8981, 9096, 8973, 9054, 8980, 8933, 9192, 8962, 9097, 8838, 9086, 9103, 8703, 9008, 9009, 9161, 8790, 9086, 9076, 9037, 8914, 8880, 9138, 8847, 8856, 8912, 8859, 9023, 9095, 8954, 9175, 8813, 9149, 9042, 8841, 9020, 9087, 9027, 9024, 9126, 8848, 9225, 9028, 8969, 9349, 9004, 9059, 8960, 9207,
      8923, 9169, 8865, 8998, 8945, 9090, 9149, 8969, 9153, 9287, 8969, 8896, 8978, 8731, 8909, 9074, 9079, 8864, 9175, 8984, 8994, 9164, 8913, 9253, 8936, 9050, 8997, 9036, 8968, 9099, 8997, 9205, 8971, 9128, 8959, 8870, 9038, 8908, 9007, 9089, 8662, 8942, 8730, 8947, 8937, 8967, 8827, 9062, 8957, 8972, 9064, 9167, 8939, 8956,
      8840, 8841, 9058, 8898, 8993, 8902, 8900, 8919, 8873, 8807, 8883, 8915, 9079, 8955, 8951, 8850, 9062, 8806, 8801, 9239, 8878, 9017, 9096, 9103, 8652, 8795, 9001, 8691, 8829, 9001, 8785, 8777, 8859, 8956, 8956, 9028, 8680, 9074, 8913, 8964, 9004, 8641, 9121, 8866, 8777, 8867, 8819, 8932, 8889, 8833, 9044, 8938, 8885, 8981,
      8934, 8684, 8670, 8971, 8871, 8941, 8924, 8903, 9053, 8876, 8973, 8847, 9099, 8587, 8919, 8818, 9035, 8772, 9071, 8795, 8933, 8895, 8896, 8812, 8737, 8977, 8798, 8939, 8708, 8731, 8838, 8809, 9022, 8983, 8945, 8817, 8838, 9032, 9015, 8943, 8695, 8779, 8839, 8873, 8628, 8770, 8989, 9022, 9001, 9141, 8943, 9031, 8882, 8813,
      8800, 8790, 8835, 8936, 8960, 8683, 8720, 8800, 9049, 8686, 8651, 9002, 8927, 8874, 8896, 8805, 8870, 8792, 8626, 8783, 9082, 8973, 8859, 8920, 8611, 8807, 8789, 8926, 8976, 8823, 8653, 8673, 8689, 8729, 8679, 8906, 8701, 8781, 8810, 8928, 8629, 8860, 8705, 8767, 8660, 8702, 8736, 8700, 8781, 8986, 9037, 9018, 8970, 9006,
      8714, 8737, 8855, 8815, 8855, 8632, 8634, 8864, 8965, 8636, 8732, 8851, 8911, 8848, 8910, 8696, 8762, 8800, 8797, 8980, 8827, 8612, 8812, 8530, 8756, 8722, 8926, 8686, 8600, 8778, 8992, 8646, 8844, 8852, 8826, 8640, 8877, 8790, 8665, 8656, 8892, 8688, 8989, 8599, 8664, 8795, 8553, 8832, 8791, 8658, 8790, 8950, 8748, 8897,
      8728, 8899, 8993, 8802, 8691, 8606, 8552, 8836, 8462, 8753, 8771, 8850, 8659, 8872, 8936, 8558, 8514, 8813, 8683, 8654, 8638, 8447, 8722, 8674, 8774, 8850, 8946, 8656, 9014, 8591, 8582, 8714, 8760, 8600, 8559, 8807, 8638, 8667, 8534, 9086, 8541, 8595, 8656, 8686, 8650, 8814, 8681, 8625, 8902, 8708, 8943, 8767, 8772, 8713,
      8623, 8777, 8564, 8772, 8668, 8865, 8753, 8531, 8648, 8854, 8716, 8476, 8364, 8697, 8711, 8830, 8729, 8748, 8678, 8688, 8816, 8644, 8713, 8454, 8659, 8670, 8712, 8375, 8670, 8811, 8662, 8547, 8666, 8952, 8793, 8659, 8743, 8894, 8663, 8564, 8726, 8983, 8674, 8698, 8721, 8733, 8552, 8580, 8950, 8600, 8538, 8933, 8664, 8535,
      8749, 8839, 8668, 8570, 8585, 8529, 8662, 8676, 8603, 8632, 8749, 8779, 8684, 8657, 8627, 8686, 8551, 8668, 8727, 8734, 8868, 8696, 8562, 8797, 8613, 8744, 8644, 8558, 8632, 8694, 8659, 8852, 8565, 8623, 8730, 8615, 8762, 8677, 8485, 8802, 8616, 8619, 8660, 8759, 8636, 8583, 8717, 8671, 8612, 8457, 8610, 8450, 8550, 8632,
      8462, 8599, 8664, 8384, 8593, 8591, 8399, 8497, 8517, 8614, 8628, 8644, 8601, 8668, 8794, 8705, 8743, 8831, 8542, 8572, 8664, 8632, 8446, 8514, 8762, 8608, 8569, 8607, 8437, 8751, 8580, 8606, 8557, 8474, 8536, 8739, 8408, 8651, 8778, 8580, 8575, 8606, 8631, 8517, 8951, 8784, 8606, 8715, 8589, 8541, 8625, 8791, 8451, 8607,
      8736, 8500, 8644, 8547, 8758, 8457, 8386, 8655, 8669, 8502, 8633, 8716, 8682, 8622, 8648, 8917, 8710, 8774, 8568, 8641, 8467, 8589, 8667, 8407, 8576, 8858, 8696, 8678, 8790, 8641, 8676, 8736, 8672, 8550, 8516, 8387, 8744, 8470, 8733, 8637, 8644, 8524, 8755, 8604, 8512, 8560, 8729, 8572, 8556, 8489, 8437, 8552, 8514, 8565,
      8831, 8483, 8611, 8498, 8532, 8591, 8629, 8536, 8452, 8459, 8523, 8681, 8613, 8628, 8557, 8648, 8522, 8383, 8361, 8618, 8572, 8512, 8436, 8577, 8595, 8397, 8711, 8607, 8626, 8619, 8549, 8583, 8812, 8518, 8317, 8671, 8537, 8746, 8506, 8648, 8447, 8366, 8673, 8588, 8647, 8793, 8613, 8528, 8341, 8540, 8462, 8401, 8492, 8472,
      8441, 8569, 8487, 8436, 8434, 8348, 8670, 8408, 8379, 8502, 8433, 8572, 8273, 8440, 8417, 8571, 8666, 8487, 8588, 8482, 8460, 8547, 8519, 8585, 8626, 8347, 8565, 8464, 8337, 8677, 8459, 8430, 8593, 8134, 8565, 8541, 8416, 8422, 8484, 8709, 8550, 8363, 8585, 8435, 8593, 8499, 8454, 8455, 8345, 8718, 8607, 8527, 8445, 8515,
      8528, 8459, 8549, 8495, 8351, 8465, 8522, 8561, 8321, 8443, 8284, 8414, 8482, 8461, 8341, 8718, 8384, 8390, 8462, 8417, 8697, 8353, 8617, 8576, 8519, 8483, 8700, 8670, 8241, 8635, 8554, 8657, 8448, 8570, 8650, 8510, 8547, 8458, 8341, 8628, 8586, 8539, 8479, 8360, 8490, 8730, 8452, 8598, 8542, 8255, 8429, 8485, 8596, 8465,
      8372, 8591, 8510, 8377, 8359, 8509, 8380, 8437, 8337, 8701, 8354, 8341, 8446, 8763, 8450, 8498, 8539, 8725, 8588, 8168, 8580, 8620, 8456, 8535, 8564, 8576, 8464, 8513, 8445, 8241, 8726, 8416, 8358, 8693, 8541, 8560, 8482, 8532, 8366, 8536, 8479, 8333, 8493, 8358, 8423, 8466, 8503, 8430, 8440, 8503, 8759, 8317, 8531, 8722,
      8492, 8367, 8488, 8601, 8340, 8521, 8717, 8579, 8410, 8383, 8409, 8295, 8322, 8382, 8691, 8355, 8336, 8583, 8452, 8449, 8518, 8454, 8384, 8290, 8495, 8537, 8149, 8381, 8418, 8391, 8578, 8379, 8315, 8440, 8651, 8410, 8361, 8411, 8190, 8484, 8539, 8485, 8541, 8430, 8536, 8578, 8359, 8295, 8408, 8423, 8553, 8269, 8810, 8347,
      8278, 8701, 8246, 8522, 8537, 8552, 8378, 8281, 8460, 8387, 8517, 8352, 8346, 8407, 8308, 8334, 8508, 8250, 8333, 8374, 8592, 8443, 8297, 8412, 8190, 8431, 8330, 8492, 8431, 8311, 8398, 8291, 8219, 8379, 8733, 8544, 8404, 8725, 8277, 8585, 8483, 8297, 8354, 8311, 8468, 8469, 8454, 8317, 8362, 8363, 8580, 8505, 8469, 8227,
      8176, 8632, 8505, 8281, 8751, 8413, 8434, 8234, 8280, 8540, 8442, 8492, 8415, 8439, 8448, 8599, 8400, 8389, 8395, 8520, 8598, 8530, 8366, 8353, 8363, 8272, 8612, 8343, 8463, 8339, 8515, 8456, 8179, 8399, 8399, 8501, 8414, 8342, 8233, 8195, 8539, 8428, 8366, 8563, 8440, 8183, 8285, 8373, 8329, 8486, 8497, 8653, 8342, 8226,
      8405, 8396, 8233, 8424, 8425, 8405, 8419, 8337, 8344, 8231, 8313, 8539, 8387, 8339, 8307, 8403, 8412, 8040, 8473, 8433, 8294, 8308, 8561, 8219, 8341, 8446, 8290, 8492, 8224, 8305, 8432, 8317, 8448, 8285, 8575, 8305, 8313, 8585, 8275, 8375, 8292, 8445, 8463, 8457, 8282, 8400, 8516, 8314, 8378, 8584, 8338, 8362, 8438, 8193,
      8225, 8254, 8256, 8530, 8312, 8278, 7995, 8283, 8511, 8361, 8471, 8420, 8271, 8385, 8200, 8427, 8299, 8359, 8409, 8064, 8601, 8299, 8304, 8365, 8264, 8372, 8379, 8341, 8210, 8465, 8470, 8229, 8497, 8265, 8163, 8304, 8432, 8285, 8452, 8407, 8134, 8409, 8337, 8286, 8439, 8330, 8389, 8245, 8403, 8184, 7970, 8297, 8193, 8334,
      8490, 8215, 8247, 8295, 8272, 8360, 8116, 8250, 8220, 8365, 8439, 8423, 8282, 8393, 8171, 8195, 8038, 8085, 8278, 8440, 8523, 8241, 8187, 8343, 8354, 8302, 8176, 8092, 8199, 8427, 8414, 8220, 8371, 8239, 8443, 8228, 8049, 8372, 8475, 8398, 8079, 8380, 8361, 8277, 8128, 8234, 8329, 8419, 8285, 8339, 8354, 8418, 8440, 8483,
      8425, 8474, 8393, 8255, 8342, 8303, 8211, 8200, 8228, 8232, 8126, 8153, 8352, 8552, 8330, 8287, 8343, 8188, 8627, 8480, 8343, 8206, 8525, 8173, 8314, 8273, 8260, 8384, 8396, 8207, 8328, 8353, 8283, 8208, 8353, 8298, 8252, 8341, 8377, 8474, 8395, 8470, 8371, 8376, 8342, 8103, 8235, 8269, 8234, 8287, 8310, 8353, 8437, 8257,
      8230, 8275, 8414, 8450, 8368, 8313, 8192, 8208, 8071, 8174, 8307, 8081, 8113, 8326, 8138, 8329, 8341, 8365, 8139, 8353, 8188, 8157, 8395, 8297, 8092, 8163, 8369, 8150, 8288, 8399, 8037, 8077, 8127, 8119, 8066, 8207, 8149, 8271, 8258, 8343, 8267, 8204, 8305, 8313, 8359, 8124, 8274, 8268, 8215, 8197, 8086, 8449, 8267, 8229
    };

    vector <double> invpdf    = get_inverse_pdf(occ);

    float rand_sampler  = rand() / float(RAND_MAX);
    float pT_gen        = iGen->pt();
    float m_gen         = iGen->mass();
    // std::cout << "   Testing with pT = " << pT_gen << " and m = " << m_gen << std::endl;
    float wgt           = lookup_invpdf(m_gen, m_bins, pT_gen, pT_bins, invpdf);
    // if (debug) std::cout << "      wgt " << wgt  << " | rand_sampler " << rand_sampler << std::endl;
    if (rand_sampler > wgt)
      {
      pass_unbaising = false;
      }
   } // unbaising loop



  TLorentzVector GenTau1  = SetTaus(iGen->daughter(0)->pt(), iGen->daughter(0)->eta(), iGen->daughter(0)->phi(), iGen->daughter(0)->mass());
  TLorentzVector GenTau2  = SetTaus(iGen->daughter(1)->pt(), iGen->daughter(1)->eta(), iGen->daughter(1)->phi(), iGen->daughter(1)->mass());
  TLorentzVector GenA1 = GenTau1 + GenTau2;


  genA1_mass_inv = GenA1.M();

  genA1_mass = iGen->mass();




  A1_Tau1_dR = dR_A1_Tau1;
  A1_Tau2_dR = dR_A1_Tau2;
  Tau1_Tau2_dR = dR_Tau1_Tau2;

  A1_pt = iGen->pt();
  Tau1_pt = iGen->daughter(0)->pt();
  Tau2_pt = iGen->daughter(1)->pt();
  A1_eta = iGen->eta();
  Tau1_eta = iGen->daughter(0)->eta();
  Tau2_eta = iGen->daughter(1)->eta();
  A1_phi = iGen->phi();
  Tau1_phi = iGen->daughter(0)->phi();
  Tau2_phi = iGen->daughter(1)->phi();

  Tau1_Tau2_deta = abs(Tau1_eta-Tau2_eta);
  Tau1_Tau2_dphi = abs(Tau1_phi-Tau2_phi);


  std::cout << "  >>>>>> A (25)/ U (553) gen  <<<<<"<<"<<< status: "<<iGen->status()<<"<<<pt:  "<<iGen->pt()<<"<<<eta:  "<<iGen->eta()<<"<<<phi:  "<<iGen->phi()<<"<<<mass:  "<<iGen->mass() <<"<<< pdgid:  "<<iGen->pdgId()<< std::endl;
  std::cout << "  >>>>>> Gen A (25)/ U (553) LorentzVector  <<<<<"<<"<<<mass:  "<< GenA1.M() << std::endl;
  std::cout << "  >>>>>> Tau1 gen (15) <<<<<"<<"<<< status: "<<iGen->daughter(0)->status()<<"<<<pt:  "<<iGen->daughter(0)->pt()<<"<<<eta:  "<<iGen->daughter(0)->eta()<<"<<<phi:  "<<iGen->daughter(0)->phi()<<"<<<mass:  "<<iGen->daughter(0)->mass() << std::endl;
  std::cout << "  >>>>>> Tau2 gen (15) <<<<<"<<"<<< status: "<<iGen->daughter(1)->status()<<"<<<pt:  "<<iGen->daughter(1)->pt()<<"<<<eta:  "<<iGen->daughter(1)->eta()<<"<<<phi:  "<<iGen->daughter(1)->phi()<<"<<<mass:  "<<iGen->daughter(1)->mass() << std::endl;
  std::cout << "  >>>>>> dR_A1_Tau1:  "<<dR_A1_Tau1<< std::endl;
  std::cout << "  >>>>>> dR_A1_Tau2:  "<<dR_A1_Tau2<< std::endl;
  std::cout << "  >>>>>> dR_Tau1_Tau2:  "<<dR_Tau1_Tau2<< std::endl;

  // }

if (pass_gen && pass_unbaising){
npassed_event++;


V_att_genA1_M_inv    = genA1_mass_inv;
V_att_genA1_M    = genA1_mass;
V_att_dR_A1_Tau1    = A1_Tau1_dR;
V_att_dR_A1_Tau2    = A1_Tau2_dR;
V_att_dR_Tau1_Tau2    = Tau1_Tau2_dR;

V_att_A1_pt    = A1_pt;
V_att_Tau1_pt    = Tau1_pt;
V_att_Tau2_pt    = Tau2_pt;
V_att_A1_eta    = A1_eta;
V_att_Tau1_eta    = Tau1_eta;
V_att_Tau2_eta    = Tau2_eta;
V_att_A1_phi    = A1_phi;
V_att_Tau1_phi    = Tau1_phi;
V_att_Tau2_phi    = Tau2_phi;

V_att_Tau1_Tau2_deta    = Tau1_Tau2_deta;
V_att_Tau1_Tau2_dphi    = Tau1_Tau2_dphi;

V_att_genA1_M_inv_.push_back( V_att_genA1_M_inv );
V_att_genA1_M_.push_back( V_att_genA1_M );
V_att_dR_A1_Tau1_.push_back( V_att_dR_A1_Tau1 );
V_att_dR_A1_Tau2_.push_back( V_att_dR_A1_Tau2 );
V_att_dR_Tau1_Tau2_.push_back( V_att_dR_Tau1_Tau2 );

V_att_A1_pt_.push_back( V_att_A1_pt );
V_att_Tau1_pt_.push_back( V_att_Tau1_pt );
V_att_Tau2_pt_.push_back( V_att_Tau2_pt );
V_att_A1_eta_.push_back( V_att_A1_eta );
V_att_Tau1_eta_.push_back( V_att_Tau1_eta );
V_att_Tau2_eta_.push_back( V_att_Tau2_eta );
V_att_A1_phi_.push_back( V_att_A1_phi );
V_att_Tau1_phi_.push_back( V_att_Tau1_phi );
V_att_Tau2_phi_.push_back( V_att_Tau2_phi );

V_att_Tau1_Tau2_deta_.push_back( V_att_Tau1_Tau2_deta );
V_att_Tau1_Tau2_dphi_.push_back( V_att_Tau1_Tau2_dphi );




// // for ( unsigned iG(0); iG != vAIdxs.size(); ++iG ) {
//
H_tau_att_genA1_M_inv->Fill( V_att_genA1_M_inv);
H_tau_att_genA1_M->Fill( V_att_genA1_M);
H_tau_att_dR_A1_Tau1->Fill( V_att_dR_A1_Tau1);
H_tau_att_dR_A1_Tau2->Fill( V_att_dR_A1_Tau2);
H_tau_att_dR_Tau1_Tau2->Fill( V_att_dR_Tau1_Tau2);

H_tau_att_A1_pt->Fill( V_att_A1_pt);
H_tau_att_Tau1_pt->Fill( V_att_Tau1_pt);
H_tau_att_Tau2_pt->Fill( V_att_Tau2_pt);
H_tau_att_A1_eta->Fill( V_att_A1_eta);
H_tau_att_Tau1_eta->Fill( V_att_Tau1_eta);
H_tau_att_Tau2_eta->Fill( V_att_Tau2_eta);
H_tau_att_A1_phi->Fill( V_att_A1_phi);
H_tau_att_Tau1_phi->Fill( V_att_Tau1_phi);
H_tau_att_Tau2_phi->Fill( V_att_Tau2_phi);

H_tau_att_Tau1_Tau2_deta->Fill( V_att_Tau1_Tau2_deta);
H_tau_att_Tau1_Tau2_dphi->Fill( V_att_Tau1_Tau2_dphi);

H_tau_att_Tau1_Tau2_dphi_deta->Fill(V_att_Tau1_Tau2_dphi,V_att_Tau1_Tau2_deta);

RHTree->Fill();

}//pass loop
}//Gen loop


V_att_genA1_M_inv_.clear();
V_att_genA1_M_.clear();
V_att_dR_A1_Tau1_.clear();
V_att_dR_A1_Tau2_.clear();
V_att_dR_Tau1_Tau2_.clear();

V_att_A1_pt_.clear();
V_att_Tau1_pt_.clear();
V_att_Tau2_pt_.clear();
V_att_A1_eta_.clear();
V_att_Tau1_eta_.clear();
V_att_Tau2_eta_.clear();
V_att_A1_phi_.clear();
V_att_Tau1_phi_.clear();
V_att_Tau2_phi_.clear();

V_att_Tau1_Tau2_deta_.clear();
V_att_Tau1_Tau2_dphi_.clear();


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
GenAnalyzer::beginJob()
{
  ntotal_event = 0;
  npassed_event = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void
GenAnalyzer::endJob()

{


  std::cout << "  >>>>>> Total events selected events <<<<<  "<<npassed_event<<"/"<<ntotal_event<<std::endl;


}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenAnalyzer);
