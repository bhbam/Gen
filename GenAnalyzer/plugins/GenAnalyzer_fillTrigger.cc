#include "Gen/GenAnalyzer/interface/GenAnalyzer.h"
#include "FWCore/Utilities/interface/RegexMatch.h"
#include "FWCore/Utilities/interface/RegexMatch.h"

TH1D *H_accept_trigger;
float V_accept_trigger;

//-----------------------now do what ever initialization is needed

void GenAnalyzer::branchesTrigger(TTree* tree, edm::Service<TFileService> &fs)
{
  H_accept_trigger     = fs->make<TH1D>("h_accept_trigger"   , "accept_trigger;accept_trigger;Events"                 ,  100,  0, 500);
  tree->Branch("trigger",  &V_accept_trigger);
  // RHTree->Branch("accepted_trigger_name",  &V_accept_trigger_name_);
}

// ---------------------- Fill tree with trigger info  ------------
void GenAnalyzer::fillTrigger(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
  V_accept_trigger = -1111.1111;
  if (isDebug) {std::cout << " >>>>>> Checking TriggerResults" << std::endl;}
  // Study of trigger bit
  edm::Handle<edm::TriggerResults> hltresults;
  iEvent.getByToken(triggerResultsToken_,hltresults);

  if (!hltresults.isValid()) {
    if(print_trigger) std::cout << "!!! Error in getting TriggerResults product from Event !!!" << std::endl;
  }
  int hltAccept = 0;
  edm::TriggerNames const& triggerNames = iEvent.triggerNames(*hltresults);
  std::string used_trgName = "HLT_DoubleMediumDeepTauPFTauHPS*";
  // std::string used_trgName = "HLT_*";

  std::vector< std::vector<std::string>::const_iterator > trgMatches = edm::regexMatch( triggerNames.triggerNames(), used_trgName );

  if ( !trgMatches.empty() ) {
    if (isDebug){std::cout << " Number of matches trugger with string : "<< used_trgName <<"---" << trgMatches.size() << std::endl;}

    for ( auto const& iT : trgMatches ) {
	    if (print_trigger){std::cout << "["<<triggerNames.triggerIndex(*iT)<<"]:"<< *iT << std::endl;}
      if ( hltresults->accept(triggerNames.triggerIndex(*iT)) ){
  	    hltAccept = hltAccept+1;
        if (print_trigger){std::cout << " name["<<triggerNames.triggerIndex(*iT)<<"]:"<< *iT << " -> " << hltresults->accept(triggerNames.triggerIndex(*iT)) << std::endl;}
	      }
      }
    }


  if (hltAccept > 0){
    if (isDebug) {std::cout << "*************** hltAccept:" << hltAccept << std::endl;}
    V_accept_trigger = hltAccept;
  }

H_accept_trigger->Fill(V_accept_trigger);

}
