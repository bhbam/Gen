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
bool unbiasing = false;

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

    std::vector <int> pT_bins   = {35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200};
    std::vector <double> m_bins = {4.0, 4.4, 4.8, 5.2, 5.6, 6.0, 6.4, 6.8, 7.2, 7.6, 8.0, 8.4, 8.8, 9.2, 9.6, 10.0, 10.4, 10.8, 11.2, 11.6, 12.0, 12.4, 12.8, 13.2, 13.6, 14.0, 14.4, 14.8, 15.2, 15.6, 16.0, 16.4, 16.8, 17.2, 17.6, 18, 18.4, 18.8, 19.2, 19.6, 20 };
    std::vector <int> occ =
    {
      0, 0, 63, 572, 1012, 1150, 1226, 1193, 1227, 1272, 1157, 1256, 1199, 1257, 1189, 1232, 1255, 1185, 1198, 1235, 1248, 1242, 1238, 1199, 1316, 1216, 1207, 1242, 1276, 1274, 1212, 1220, 1123, 1273,
      0, 45, 536, 1007, 1342, 1405, 1427, 1427, 1354, 1366, 1397, 1366, 1334, 1341, 1381, 1425, 1466, 1485, 1363, 1422, 1374, 1394, 1419, 1371, 1382, 1458, 1409, 1389, 1463, 1443, 1440, 1422, 1436, 1420,
      0, 208, 761, 1156, 1305, 1378, 1336, 1433, 1381, 1441, 1384, 1349, 1381, 1472, 1418, 1472, 1341, 1263, 1409, 1364, 1393, 1410, 1358, 1449, 1427, 1353, 1414, 1366, 1399, 1428, 1361, 1357, 1348, 1395,
      0, 369, 913, 1239, 1300, 1404, 1287, 1401, 1420, 1342, 1344, 1288, 1366, 1337, 1421, 1282, 1363, 1290, 1339, 1402, 1314, 1401, 1324, 1327, 1304, 1335, 1428, 1284, 1336, 1408, 1321, 1380, 1401, 1306,
      8, 472, 852, 1221, 1300, 1338, 1247, 1404, 1229, 1347, 1298, 1279, 1269, 1420, 1295, 1353, 1325, 1307, 1343, 1362, 1364, 1313, 1392, 1341, 1434, 1287, 1314, 1363, 1368, 1275, 1381, 1291, 1329, 1436,
      35, 556, 986, 1190, 1200, 1303, 1320, 1316, 1455, 1437, 1270, 1302, 1347, 1317, 1289, 1333, 1311, 1394, 1297, 1385, 1297, 1339, 1325, 1246, 1339, 1321, 1347, 1350, 1379, 1321, 1317, 1392, 1307, 1297,
      84, 590, 908, 1209, 1254, 1330, 1418, 1313, 1322, 1333, 1391, 1294, 1385, 1306, 1317, 1393, 1353, 1388, 1380, 1269, 1405, 1237, 1249, 1381, 1343, 1286, 1340, 1366, 1270, 1478, 1375, 1258, 1264, 1380,
      153, 637, 966, 1147, 1354, 1293, 1347, 1356, 1302, 1333, 1298, 1454, 1308, 1353, 1245, 1275, 1320, 1208, 1330, 1349, 1382, 1428, 1281, 1348, 1250, 1287, 1332, 1278, 1307, 1298, 1323, 1380, 1239, 1279,
      128, 690, 1037, 1153, 1172, 1281, 1300, 1278, 1356, 1240, 1235, 1355, 1348, 1241, 1322, 1317, 1347, 1303, 1310, 1253, 1239, 1325, 1288, 1378, 1246, 1346, 1290, 1309, 1369, 1278, 1319, 1420, 1246, 1302,
      151, 646, 1041, 1201, 1234, 1306, 1283, 1394, 1307, 1267, 1275, 1313, 1268, 1318, 1320, 1309, 1312, 1275, 1295, 1263, 1305, 1266, 1249, 1290, 1351, 1304, 1255, 1220, 1212, 1324, 1467, 1321, 1253, 1298,
      198, 678, 911, 1189, 1260, 1277, 1220, 1241, 1321, 1321, 1197, 1294, 1267, 1256, 1268, 1231, 1289, 1270, 1281, 1259, 1272, 1216, 1316, 1294, 1298, 1186, 1357, 1352, 1311, 1325, 1292, 1251, 1312, 1358,
      227, 746, 942, 1120, 1215, 1346, 1247, 1292, 1247, 1247, 1306, 1267, 1224, 1331, 1275, 1240, 1311, 1271, 1216, 1303, 1342, 1218, 1373, 1347, 1294, 1293, 1267, 1344, 1320, 1276, 1347, 1254, 1242, 1271,
      202, 746, 955, 1193, 1214, 1230, 1173, 1263, 1318, 1316, 1205, 1223, 1294, 1233, 1237, 1287, 1209, 1291, 1299, 1232, 1243, 1305, 1253, 1275, 1195, 1155, 1315, 1223, 1274, 1324, 1231, 1160, 1300, 1346,
      230, 713, 1049, 1106, 1180, 1243, 1277, 1221, 1240, 1248, 1216, 1225, 1218, 1231, 1234, 1237, 1282, 1220, 1199, 1261, 1252, 1307, 1283, 1234, 1262, 1188, 1192, 1237, 1269, 1218, 1235, 1333, 1271, 1183,
      289, 710, 969, 1166, 1183, 1283, 1268, 1276, 1348, 1270, 1226, 1297, 1238, 1219, 1272, 1324, 1207, 1283, 1284, 1328, 1183, 1127, 1312, 1316, 1225, 1139, 1291, 1233, 1359, 1242, 1119, 1327, 1264, 1254,
      295, 780, 934, 1089, 1250, 1245, 1245, 1277, 1164, 1214, 1304, 1154, 1158, 1246, 1197, 1291, 1242, 1174, 1256, 1313, 1259, 1264, 1279, 1281, 1271, 1268, 1268, 1237, 1293, 1203, 1274, 1269, 1145, 1200,
      290, 744, 959, 1114, 1226, 1183, 1196, 1152, 1235, 1250, 1237, 1266, 1249, 1284, 1164, 1310, 1210, 1186, 1247, 1230, 1268, 1180, 1219, 1346, 1206, 1211, 1254, 1309, 1208, 1261, 1246, 1328, 1216, 1169,
      289, 744, 1003, 1104, 1145, 1199, 1297, 1156, 1296, 1161, 1256, 1093, 1249, 1175, 1183, 1194, 1274, 1209, 1249, 1189, 1333, 1311, 1250, 1235, 1230, 1241, 1167, 1276, 1248, 1158, 1203, 1243, 1198, 1293,
      312, 728, 1034, 1142, 1239, 1144, 1286, 1299, 1237, 1193, 1233, 1240, 1204, 1248, 1203, 1124, 1184, 1186, 1177, 1172, 1222, 1219, 1163, 1188, 1148, 1177, 1229, 1247, 1255, 1268, 1198, 1294, 1211, 1217,
      348, 727, 981, 1104, 1204, 1217, 1172, 1230, 1290, 1192, 1136, 1187, 1232, 1266, 1276, 1168, 1234, 1231, 1328, 1256, 1186, 1152, 1207, 1205, 1141, 1225, 1199, 1156, 1238, 1182, 1230, 1167, 1152, 1211,
      357, 774, 965, 1089, 1125, 1144, 1161, 1205, 1204, 1217, 1228, 1139, 1265, 1216, 1238, 1160, 1201, 1196, 1217, 1237, 1128, 1193, 1204, 1245, 1233, 1248, 1261, 1187, 1131, 1209, 1205, 1190, 1210, 1189,
      337, 745, 957, 1092, 1218, 1130, 1204, 1101, 1157, 1200, 1196, 1271, 1258, 1183, 1153, 1152, 1171, 1230, 1226, 1204, 1241, 1273, 1209, 1188, 1235, 1151, 1079, 1239, 1207, 1227, 1137, 1101, 1146, 1163,
      383, 811, 950, 1034, 1170, 1141, 1182, 1275, 1217, 1204, 1264, 1152, 1114, 1185, 1129, 1145, 1324, 1231, 1166, 1199, 1194, 1183, 1131, 1205, 1193, 1129, 1141, 1138, 1154, 1262, 1168, 1101, 1242, 1205,
      382, 708, 989, 1092, 1256, 1195, 1130, 1098, 1207, 1165, 1213, 1164, 1188, 1238, 1113, 1213, 1173, 1260, 1128, 1235, 1232, 1238, 1128, 1255, 1170, 1169, 1229, 1166, 1188, 1243, 1141, 1155, 1187, 1307,
      379, 715, 985, 1053, 1127, 1094, 1167, 1078, 1173, 1154, 1101, 1169, 1228, 1180, 1199, 1253, 1170, 1207, 1204, 1160, 1123, 1103, 1176, 1183, 1116, 1147, 1160, 1209, 1166, 1195, 1153, 1189, 1244, 1191,
      397, 671, 928, 1099, 1154, 1248, 1198, 1153, 1174, 1226, 1163, 1145, 1174, 1164, 1260, 1195, 1188, 1072, 1144, 1200, 1174, 1233, 1182, 1214, 1152, 1104, 1155, 1171, 1111, 1133, 1167, 1225, 1168, 1125,
      380, 751, 944, 1069, 1079, 1205, 1200, 1241, 1175, 1137, 1144, 1222, 1182, 1173, 1166, 1193, 1138, 1144, 1156, 1228, 1187, 1156, 1118, 1149, 1236, 1188, 1200, 1131, 1163, 1208, 1193, 1186, 1246, 1218,
      340, 748, 961, 1073, 1073, 1144, 1134, 1185, 1142, 1233, 1137, 1256, 1262, 1273, 1206, 1107, 1246, 1186, 1137, 1152, 1114, 1147, 1194, 1226, 1131, 1199, 1203, 1182, 1157, 1171, 1169, 1171, 1174, 1202,
      398, 769, 924, 1128, 1072, 1161, 1163, 1259, 1198, 1095, 1161, 1222, 1241, 1184, 1185, 1206, 1200, 1134, 1167, 1202, 1225, 1180, 1168, 1189, 1156, 1103, 1188, 1206, 1170, 1150, 1137, 1198, 1190, 1126,
      448, 756, 974, 1081, 1072, 1120, 1160, 1151, 1190, 1136, 1185, 1140, 1205, 1159, 1135, 1159, 1120, 1139, 1141, 1128, 1106, 1222, 1138, 1148, 1140, 1155, 1167, 1230, 1135, 1144, 1151, 1198, 1196, 1118,
      449, 727, 949, 1004, 1095, 1057, 1188, 1166, 1129, 1225, 1162, 1128, 1213, 1136, 1169, 1198, 1169, 1183, 1106, 1155, 1089, 1153, 1225, 1174, 1137, 1154, 1137, 1186, 1131, 1077, 1129, 1181, 1142, 1257,
      405, 767, 967, 1143, 1181, 1162, 1237, 1133, 1126, 1192, 1166, 1163, 1094, 1253, 1103, 1167, 1100, 1159, 1184, 1173, 1075, 1129, 1193, 1185, 1135, 1160, 1160, 1149, 1131, 1127, 1192, 1085, 1277, 1124,
      412, 767, 956, 1100, 1131, 1127, 1157, 1046, 1175, 1209, 1116, 1162, 1134, 1134, 1208, 1214, 1187, 1177, 1169, 1115, 1155, 1230, 1175, 1119, 1109, 1193, 1122, 1170, 1171, 1167, 1156, 1166, 1113, 1129,
      453, 758, 983, 1071, 1153, 1103, 1055, 1188, 1204, 1196, 1110, 1227, 1173, 1116, 1149, 1208, 1243, 1210, 1084, 1186, 1150, 1150, 1201, 1104, 1198, 1122, 1218, 1109, 1120, 1137, 1170, 1163, 1149, 1157,
      493, 802, 913, 1072, 1134, 1134, 1121, 1083, 1078, 1145, 1211, 1167, 1055, 1184, 1186, 1130, 1110, 1114, 1162, 1123, 1095, 1224, 1147, 1163, 1092, 1114, 1154, 1150, 1151, 1105, 1234, 1190, 1108, 1129,
      461, 736, 994, 1052, 1078, 1171, 1138, 1191, 1074, 1118, 1144, 1140, 1203, 1072, 1151, 1084, 1082, 1133, 1179, 1148, 1132, 1135, 1232, 1178, 1203, 1179, 1181, 1129, 1187, 1096, 1087, 1122, 1115, 1139,
      482, 831, 955, 996, 1078, 1169, 1131, 1078, 1039, 1074, 1133, 1160, 1130, 1114, 1139, 1179, 1189, 1153, 1238, 1144, 1253, 1093, 1183, 1147, 1187, 1134, 1196, 1138, 1210, 1076, 1118, 1154, 1166, 1151,
      488, 809, 839, 1047, 1088, 1132, 1089, 1135, 1120, 1047, 1153, 1083, 1130, 1150, 1154, 1231, 1081, 1010, 1225, 1130, 1103, 1101, 1122, 1114, 1026, 1090, 1129, 1112, 1156, 1207, 1145, 1038, 1240, 1139,
      452, 780, 972, 1062, 1134, 1111, 1143, 1113, 1134, 1177, 1070, 1191, 1130, 1087, 1239, 1122, 1132, 1215, 1142, 1139, 1186, 1118, 1185, 971, 1114, 1171, 1136, 1132, 1173, 1072, 1111, 1219, 1165, 1091,
      439, 720, 961, 1059, 1134, 1075, 1097, 1182, 1180, 1183, 1084, 1056, 1030, 1184, 1240, 1110, 1133, 1166, 1072, 1124, 1130, 1133, 1170, 1084, 1134, 1150, 1101, 1156, 1166, 1149, 1196, 1090, 1122, 1161,
      511, 816, 1024, 1136, 1037, 1101, 1170, 1151, 1159, 1137, 1165, 1121, 1098, 1129, 1135, 1146, 1019, 1163, 1096, 1106, 1165, 1167, 1030, 1119, 1124, 1191, 1148, 1144, 1126, 1130, 1093, 1155, 1079, 1090
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
  std::cout << "  >>>>>> Tau2 gen (15) <<<<<"<<"<<< status: "<<iGen->daughter(1)->status()<<"<<<pt:  "<<iGen->daughter(1)->pt()<<"<<<eta:  "<<iGen->daughter(1)->eta()<<"<<<phi:  "<<iGen->daughter(1)->phi()<<"<<<mass:  "<<iGen->daughter(0)->daughter(1)->mass() << std::endl;
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
