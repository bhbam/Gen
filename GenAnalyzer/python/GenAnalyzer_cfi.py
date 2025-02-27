import FWCore.ParameterSet.Config as cms

from RecoMET.METProducers.METSignificanceParams_cfi import METSignificanceParams

fevt = cms.EDAnalyzer('GenAnalyzer'

, isDebug                        = cms.bool(False)
, print_trigger                  = cms.bool(False)

,genParticles    = cms.InputTag('genParticles',"","")
,ak8GenJets    = cms.InputTag('ak8GenJets',"","")
,genMetTrue    = cms.InputTag('genMetTrue',"","")
,hltresults = cms.InputTag('TriggerResults', "", "HLT")

)
