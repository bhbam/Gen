import FWCore.ParameterSet.Config as cms

from RecoMET.METProducers.METSignificanceParams_cfi import METSignificanceParams

fevt = cms.EDAnalyzer('GenAnalyzer'

,isDebug       = cms.bool(False)
,print_trigger = cms.bool(False)
,minJetPt  = cms.double(20.)
,maxJetEta = cms.double(2.4)
,genParticles        = cms.InputTag('genParticles',"","")
,ak8GenJets          = cms.InputTag('ak8GenJets',"","")
,genMetTrue          = cms.InputTag('genMetTrue',"","")
,hltresults          = cms.InputTag('TriggerResults', "", "HLT")
,ak4PFJetCollection  = cms.InputTag('ak4PFJets')
# ,ak4PFJetCollection  = cms.InputTag('ak8PFJetsPuppi')
,tauCollection       = cms.InputTag("hpsPFTauProducer")
,metCollection       = cms.InputTag("pfMet")

)
