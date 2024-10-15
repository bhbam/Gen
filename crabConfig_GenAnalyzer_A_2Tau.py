from CRABClient.UserUtilities import config#, getUsernameFromSiteDB
config = config()
# See parameter defintions here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile#CRAB_configuration_parameters
# To submit to crab:
# crab submit -c crabConfig_data.py
# To check job status:
# crab status -d <config.General.workArea>/<config.General.requestName># To resubmit jobs:
# crab resubmit -d <config.General.workArea>/<config.General.requestName>

# Local job directory will be created in:
# <config.General.workArea>/<config.General.requestName>
config.General.workArea = 'crab_Ato2Tau'
config.General.requestName = 'ATo2Tau_Hadronic_M3p6To20_mass_regression_GenIfo_only'
config.General.transferOutputs = True
config.General.transferLogs = False

# CMS cfg file goes here:
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'GenAnalyzer/python/conFig_cfg.py' # analyzer cfg file
# config.JobType.maxMemoryMB = 2000
config.Data.inputDBS = 'phys03'
config.JobType.allowUndistributedCMSSW = True
# Define input and units per job here:

config.Data.inputDataset = "/GEN_SIM_ATo2Tau_m3p6To20_pt30To200/bhbam-crab_aToTauTau_Hadronic_m3p6To20_pythia8_GEN_SIM-f88935475b466ffb7e3f550798b0e55c/USER"
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10 # units: as defined by config.Data.splitting
config.Data.totalUnits = -1 # -1: all inputs. total jobs submitted = totalUnits / unitsPerJob. cap of 10k jobs per submission
config.Data.publication = False

# Output files will be stored in config.Site.storageSite at directory:
# <config.Data.outLFNDirBase>/<config.Data.outputPrimaryDataset>/<config.Data.outputDatasetTag>/
config.Site.storageSite = 'T3_US_FNALLPC'
config.Data.outLFNDirBase = '/store/user/bhbam/Ntuples_run3' # add your username as subdirectory
# config.Data.outLFNDirBase = '/store/group/lpcml/bbbam/Ntuples_signal_with_trigger' # add your username as subdirectory
# config.Data.outputPrimaryDataset = 'Samples'
config.Data.outputDatasetTag = config.General.requestName
