import FWCore.ParameterSet.Config as cms

from RecoMET.METProducers.METSignificanceParams_cfi import METSignificanceParams

fevt = cms.EDAnalyzer('GenAnalyzer'

,isDebug       = cms.bool(True)
,print_trigger = cms.bool(False)
# ,selected_trgName  = cms.string("HLT_DoubleMediumDeepTauPFTauHPS*")

,genParticles  = cms.InputTag('genParticles',"","")
,ak8GenJets    = cms.InputTag('ak8GenJets',"","")
,genMetTrue    = cms.InputTag('genMetTrue',"","")
,hltresults    = cms.InputTag('TriggerResults', "", "HLT")

)
