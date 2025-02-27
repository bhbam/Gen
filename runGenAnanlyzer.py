import os
cfg='GenAnalyzer/python/conFig_cfg.py'
# inputFiles_='file:root://cmsxrootd.fnal.gov//store/group/lpcml/bbbam/MCGeneration/gen_HToAATo4Tau_Hadronic_tauDR0p4_M10_ctau0To3_eta0To2p4_pythia8_2018UL/crab_gen_HToAATo4Tau_Hadronic_tauDR0p4_M10_ctau0To3_eta0To2p4_pythia8_2018UL/230510_053835/0000/GEN_HToAAToTauTau_M10_2018UL_90.root'
# inputFiles_='file:/uscms/home/bbbam/nobackup/CMSSW_10_6_20/src/MCProduction/E2E-HToAATo4Tau/GEN_HToAAToTauTau_M13_2018UL.root'
# inputFiles_='file:/uscms/home/bbbam/nobackup/analysis/MCGeneration/CMSSW_10_6_20/src/MCProduction/E2E-HToAATo4Tau/Unboosted_GEN_HToAAToTauTau_M8_2018UL.root'
# inputFiles_='file:root://cmsxrootd.fnal.gov//store/user/bhbam/MCGeneration/Unboosted_gen_HToAATo4Tau_Hadronic_tauDR0p4_M8_ctau0To3_eta0To2p4_pythia8_2018UL_lessPerFile/crab_Unboosted_gen_HToAATo4Tau_Hadronic_tauDR0p4_M8_ctau0To3_eta0To2p4_pythia8_2018UL_lessPerFile/240103_211550/0000/Unboosted_GEN_HToAAToTauTau_M8_2018UL_15.root'
# inputFiles_='file:/uscms/home/bbbam/nobackup/analysis/MCGeneration/CMSSW_10_6_20/src/MCProduction_with_trigger/GEN_SIM_DIGI_HToAATo4Tau_M5.root'
inputFiles_='file:root://cmsxrootd.fnal.gov//store/group/lpcml/rchudasa/MCGenerationRun3/HToAATo4Tau_hadronic_tauDecay_M3p7_Run3_2023/3p7_AODSIM_newBigProd/250113_144409/0000/step3_AODSIM_97.root'
inputFiles_='file:/uscms/home/bbbam/nobackup/analysis_run3/MCGeneration/CMSSW_13_0_17/src/AOD_HToAATo4Tau_M14.root'
maxEvents_=10
skipEvents_=0#
outputFile_='GenInfo_only_H_AA_4Tau_M3p7.root'
cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
print(cmd)
os.system(cmd)
