import os
import glob
import numpy as np
from concurrent.futures import ProcessPoolExecutor

# Configuration file and input directory
cfg = "GenAnalyzer/python/conFig_cfg.py"
local = "/eos/uscms/store/group/lpcml/bbbam/MCGeneration_run3/GEN_SIM_ATo2Tau_m3p6To18_pt30To300_v2/crab_aToTauTau_Hadronic_m3p6To18_pt30To300_pythia8_GEN_SIM_v2/241030_043831/00*/*.root"
# local = "/eos/uscms/store/group/lpcml/bbbam/MCGeneration_run3/GEN_SIM_ATo2Tau_m1p2To3p6_pt30To300_v2/crab_aToTauTau_Hadronic_m1p2To3p6_pt30To300_pythia8_GEN_SIM_v2/241030_050102/000*/*.root"

# Sorting input files
inputFiles = np.sort(glob.glob(local))

# Function to execute a single cmsRun command
def run_cmsRun(input_root):
    inputFiles_ = "file:" + input_root
    tag = (input_root.split("_")[-1]).split(".")[0]
    maxEvents_ = -1
    skipEvents_ = 0
    outputFile_ = f"GenInfo_only_A_2Tau_m3p6To18_pT30To300_v2_{tag}.root"
    cmd = f"cmsRun {cfg} inputFiles={inputFiles_} maxEvents={maxEvents_} skipEvents={skipEvents_} outputFile={outputFile_}"
    print(cmd)
    os.system(cmd)

# Parallel execution using ProcessPoolExecutor
if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=8) as executor:  # Using 8 cores
        executor.map(run_cmsRun, inputFiles)
