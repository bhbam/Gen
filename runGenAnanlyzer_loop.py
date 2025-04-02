import os
import glob
import numpy as np
from concurrent.futures import ProcessPoolExecutor


Mass = '14'

# Configuration file and input directory
cfg = "GenAnalyzer/python/conFig_cfg.py"

inputFiles = []
for i in range(1,101):

    file ={
        "3p7":f"root://xrootd.unl.edu//store/group/lpcml/rchudasa/MCGenerationRun3/HToAATo4Tau_hadronic_tauDecay_M3p7_Run3_2023/3p7_AODSIM_newBigProd/250113_144409/0000/step3_AODSIM_{i}.root"
        ,"4":f"root://xrootd.unl.edu//store/group/lpcml/rchudasa/MCGenerationRun3/HToAATo4Tau_hadronic_tauDecay_M4_Run3_2023/4_AODSIM_newBigProd/250113_145716/0000/step3_AODSIM_{i}.root"
        ,"5":f"root://xrootd.unl.edu//store/group/lpcml/rchudasa/MCGenerationRun3/HToAATo4Tau_hadronic_tauDecay_M5_Run3_2023/5_AODSIM_newBigProd/250113_145836/0000/step3_AODSIM_{i}.root"
        ,"6":f"root://xrootd.unl.edu//store/group/lpcml/rchudasa/MCGenerationRun3/HToAATo4Tau_hadronic_tauDecay_M6_Run3_2023/6_AODSIM_newBigProd/250113_150335/0000/step3_AODSIM_{i}.root"
        ,"8":f"root://xrootd.unl.edu//store/group/lpcml/rchudasa/MCGenerationRun3/HToAATo4Tau_hadronic_tauDecay_M8_Run3_2023/signal_Mass_8_AODSIM_multiThreads/250111_024345/0000/step3_AODSIM_M14_{i}.root"
        ,"10":f"root://xrootd.unl.edu//store/group/lpcml/rchudasa/MCGenerationRun3/HToAATo4Tau_hadronic_tauDecay_M10_Run3_2023/signal_Mass_10_AODSIM_multiThreads/250111_024902/0000/step3_AODSIM_M14_{i}.root"
        ,"12":f"root://xrootd.unl.edu//store/group/lpcml/rchudasa/MCGenerationRun3/HToAATo4Tau_hadronic_tauDecay_M12_Run3_2023/signal_Mass_12_AODSIM_multiThreads/250111_025148/0000/step3_AODSIM_M14_{i}.root"
        ,"14":f"root://xrootd.unl.edu//store/group/lpcml/rchudasa/MCGenerationRun3/HToAATo4Tau_hadronic_tauDecay_M14_Run3_2023/signal_Mass_14_AODSIM_multiThreads/250111_025512/0000/step3_AODSIM_M14_{i}.root"
    }.get(Mass, None)
    inputFiles.append(file)
# print(inputFiles)

# Function to execute a single cmsRun command
def run_cmsRun(input_root):
    inputFiles_ = "file:" + input_root
    tag = (input_root.split("_")[-1]).split(".")[0]
    maxEvents_ = -1
    skipEvents_ = 0
    outputFile_ = f"Gen_reco_Info_with_multi_trigger_H_A_4Tau_M{Mass}_{tag}.root"
    cmd = f"cmsRun {cfg} inputFiles={inputFiles_} maxEvents={maxEvents_} skipEvents={skipEvents_} outputFile={outputFile_}"
    print(cmd)
    os.system(cmd)

# Parallel execution using ProcessPoolExecutor
if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=6) as executor:  # Using 8 cores
        executor.map(run_cmsRun, inputFiles)
# To run
# nohup python3 runGenAnanlyzer_loop.py > gen_M3p7.log 2>&1 &
# to check ps aux | grep runGenAnanlyzer_loop.py
# to pkill -9 -f runGenAnanlyzer_loop.py
