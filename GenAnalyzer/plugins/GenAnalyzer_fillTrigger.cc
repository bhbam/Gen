#include "Gen/GenAnalyzer/interface/GenAnalyzer.h"

TH1D *H_accept_trigger;
float V_accept_trigger_1;
float V_accept_trigger_2;
float V_accept_trigger_3;
float V_accept_trigger_4;

//-----------------------now do what ever initialization is needed

void GenAnalyzer::branchesTrigger(TTree* tree, edm::Service<TFileService> &fs)
{
  H_accept_trigger     = fs->make<TH1D>("h_accept_trigger"   , "accept_trigger;accept_trigger;Events"                 ,  0,  2, 1);
  tree->Branch("eta2p1_v4",  &V_accept_trigger_1);
  tree->Branch("eta2p1_PFJet60_v4",  &V_accept_trigger_2);
  tree->Branch("eta2p1_PFJet75_v4",  &V_accept_trigger_3);
  tree->Branch("eta2p1_OneProng_M5to80_v2",  &V_accept_trigger_4);
}

// ---------------------- Fill tree with trigger info  ------------
void GenAnalyzer::fillTrigger(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
  V_accept_trigger_1 = -1111.1111;
  V_accept_trigger_2 = -1111.1111;
  V_accept_trigger_3 = -1111.1111;
  V_accept_trigger_4 = -1111.1111;
  if (isDebug) {std::cout << " >>>>>> Checking TriggerResults" << std::endl;}
  // Study of trigger bit
  edm::Handle<edm::TriggerResults> hltresults;
  iEvent.getByToken(triggerResultsToken_,hltresults);

  if (!hltresults.isValid()) {
    if(print_trigger) std::cout << "!!! Error in getting TriggerResults product from Event !!!" << std::endl;
  }
  int hltAccept_1 = 0;
  int hltAccept_2 = 0;
  int hltAccept_3 = 0;
  int hltAccept_4 = 0;
  edm::TriggerNames const& triggerNames = iEvent.triggerNames(*hltresults);
  std::string used_trgName_1 = "HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1_v4";
  std::string used_trgName_2 = "HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60_v4";
  std::string used_trgName_3 = "HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75_v4";
  std::string used_trgName_4 = "HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_OneProng_M5to80_v2";
  // std::string used_trgName = "HLT_*";

  std::vector< std::vector<std::string>::const_iterator > trgMatches_1 = edm::regexMatch( triggerNames.triggerNames(), used_trgName_1 );
  std::vector< std::vector<std::string>::const_iterator > trgMatches_2 = edm::regexMatch( triggerNames.triggerNames(), used_trgName_2 );
  std::vector< std::vector<std::string>::const_iterator > trgMatches_3 = edm::regexMatch( triggerNames.triggerNames(), used_trgName_3 );
  std::vector< std::vector<std::string>::const_iterator > trgMatches_4 = edm::regexMatch( triggerNames.triggerNames(), used_trgName_4 );

  if ( !trgMatches_1.empty() ) {
    if (isDebug){std::cout << " Number of matches trugger with string : "<< used_trgName_1 <<"---" << trgMatches_1.size() << std::endl;}
  for ( auto const& iT_1 : trgMatches_1 ) {
	    if (print_trigger){std::cout << "["<<triggerNames.triggerIndex(*iT_1)<<"]:"<< *iT_1 << std::endl;}
      if ( hltresults->accept(triggerNames.triggerIndex(*iT_1)) ){
  	    hltAccept_1 = hltAccept_1+1;
        if (print_trigger){std::cout << " name["<<triggerNames.triggerIndex(*iT_1)<<"]:"<< *iT_1 << " -> " << hltresults->accept(triggerNames.triggerIndex(*iT_1)) << std::endl;}
	      }
      }
    }
  if (hltAccept_1 > 0)
  {
    if (isDebug) {std::cout << "*************** HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1_v4:" << hltAccept_1 << std::endl;}
    V_accept_trigger_1 = hltAccept_1;
  }


    if ( !trgMatches_2.empty() ) {
      if (isDebug){std::cout << " Number of matches trugger with string : "<< used_trgName_2 <<"---" << trgMatches_2.size() << std::endl;}
    for ( auto const& iT_2 : trgMatches_2 ) {
  	    if (print_trigger){std::cout << "["<<triggerNames.triggerIndex(*iT_2)<<"]:"<< *iT_2 << std::endl;}
        if ( hltresults->accept(triggerNames.triggerIndex(*iT_2)) ){
    	    hltAccept_2 = hltAccept_2+1;
          if (print_trigger){std::cout << " name["<<triggerNames.triggerIndex(*iT_2)<<"]:"<< *iT_2 << " -> " << hltresults->accept(triggerNames.triggerIndex(*iT_2)) << std::endl;}
  	      }
        }
      }
    if (hltAccept_2 > 0)
    {
      if (isDebug) {std::cout << "*************** HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60_v4:" << hltAccept_2 << std::endl;}
      V_accept_trigger_2 = hltAccept_2;
    }


    if ( !trgMatches_3.empty() ) {
      if (isDebug){std::cout << " Number of matches trugger with string : "<< used_trgName_3 <<"---" << trgMatches_3.size() << std::endl;}
    for ( auto const& iT_3 : trgMatches_3 ) {
        if (print_trigger){std::cout << "["<<triggerNames.triggerIndex(*iT_3)<<"]:"<< *iT_3 << std::endl;}
        if ( hltresults->accept(triggerNames.triggerIndex(*iT_3)) ){
          hltAccept_3 = hltAccept_3+1;
          if (print_trigger){std::cout << " name["<<triggerNames.triggerIndex(*iT_3)<<"]:"<< *iT_3 << " -> " << hltresults->accept(triggerNames.triggerIndex(*iT_3)) << std::endl;}
          }
        }
      }
    if (hltAccept_3 > 0)
    {
      if (isDebug) {std::cout << "*************** HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75_v4:" << hltAccept_3 << std::endl;}
      V_accept_trigger_3 = hltAccept_3;
    }


    if ( !trgMatches_4.empty() ) {
      if (isDebug){std::cout << " Number of matches trugger with string : "<< used_trgName_4 <<"---" << trgMatches_4.size() << std::endl;}
    for ( auto const& iT_4 : trgMatches_4 ) {
        if (print_trigger){std::cout << "["<<triggerNames.triggerIndex(*iT_4)<<"]:"<< *iT_4 << std::endl;}
        if ( hltresults->accept(triggerNames.triggerIndex(*iT_4)) ){
          hltAccept_4 = hltAccept_4+1;
          if (print_trigger){std::cout << " name["<<triggerNames.triggerIndex(*iT_4)<<"]:"<< *iT_4 << " -> " << hltresults->accept(triggerNames.triggerIndex(*iT_4)) << std::endl;}
          }
        }
      }
    if (hltAccept_4 > 0)
    {
      if (isDebug) {std::cout << "*************** HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_OneProng_M5to80_v2:" << hltAccept_4 << std::endl;}
      V_accept_trigger_4 = hltAccept_4;
    }


  if ((hltAccept_4 > 0) || (hltAccept_4 > 0) || (hltAccept_4 > 0) || (hltAccept_4 > 0))
  {
    H_accept_trigger->Fill(1);
  }
}
