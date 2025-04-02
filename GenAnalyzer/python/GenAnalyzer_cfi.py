import FWCore.ParameterSet.Config as cms

from RecoMET.METProducers.METSignificanceParams_cfi import METSignificanceParams

fevt = cms.EDAnalyzer('GenAnalyzer'

, isDebug                        = cms.bool(False)


,genParticles    = cms.InputTag('genParticles',"","")


)
