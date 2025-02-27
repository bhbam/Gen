import os
import glob
import numpy as np
from concurrent.futures import ProcessPoolExecutor

# Configuration file and input directory
cfg = "GenAnalyzer/python/conFig_cfg.py"

inputFiles = []

for i in range(1,101):
    # file = f"root://xrootd.unl.edu//store/group/phys_diffraction/rchudasa/MCGeneration/HToAATo4Tau_M3p7_Run3_2023/3p7_AODSIM/240911_121009/0000/step3_AODSIM_{i}.root"
    file = f"root://xrootd.unl.edu//store/group/phys_diffraction/rchudasa/MCGeneration/HToAATo4Tau_M14_Run3_2023/14_AODSIM/240911_120940/0000/step3_AODSIM_{i}.root"
    inputFiles.append(file)
# print(inputFiles)

# Function to execute a single cmsRun command
def run_cmsRun(input_root):
    inputFiles_ = "file:" + input_root
    tag = (input_root.split("_")[-1]).split(".")[0]
    maxEvents_ = -1
    skipEvents_ = 0
    outputFile_ = f"GenInfo_H_A_4Tau_M14_{tag}.root"
    cmd = f"cmsRun {cfg} inputFiles={inputFiles_} maxEvents={maxEvents_} skipEvents={skipEvents_} outputFile={outputFile_}"
    print(cmd)
    os.system(cmd)

# Parallel execution using ProcessPoolExecutor
if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=8) as executor:  # Using 8 cores
        executor.map(run_cmsRun, inputFiles)
