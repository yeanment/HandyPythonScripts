#! /bin/bash
python tecisoZBilger_v2.py --timeStart 0. --timeEnd -1  --writeTimeStepCombine True --case  /mnt/j/CFDiffH2/TTT/CFDiffH2-L_0.05-UF_01.00-Zst_0.50-D4.0E-4-T2.0E-4/postProcessing/isoZBilgerSample --isoSdSurface T --isoSdSurfaceValue 400 500 600 700 800 898.0 900 1000 1100 1144.8 1200 1300 1400 1500 1600  --keepDataOption 1 --averageWidth 11 --savgolFilterWidth 25
python tecisoZBilger_v2.py --timeStart 0. --timeEnd -1  --writeTimeStepCombine True --case  /mnt/j/CFDiffH2/TTT/CFDiffH2-L_0.05-UF_01.00-Zst_0.50-D4.0E-4-T2.0E-4/postProcessing/isoZBilgerSample --isoSdSurface T --isoSdSurfaceValue 400 500 600 700 800 898.0 900 1000 1100 1144.8 1200 1300 1400 1500 1600  --keepDataOption 1 --averageWidth 33 --savgolFilterWidth 25

# python tecisoZBilger_v2.py --timeStart 0. --timeEnd -1  --writeTimeStepCombine True --case  /mnt/j/CFDiffH2/TTT/CFDiffH2-L_0.05-UF_01.00-Zst_0.50-D4.0E-4-T2.0E-4/postProcessing/isoZBilgerSample --isoSdSurface H2O --isoSdSurfaceValue 0.01 0.02 0.03 0.04 0.05 0.05335953 0.06 0.07 0.08 0.08488209 0.09 0.10 0.11 0.012 0.13  --keepDataOption 1 --averageWidth 16 --savgolFilterWidth 25
# python tecisoZBilger_v2.py --timeStart 0. --timeEnd -1  --writeTimeStepCombine True --case  /mnt/j/CFDiffH2/TTT/CFDiffH2-L_0.05-UF_01.00-Zst_0.50-D4.0E-4-T2.0E-4/postProcessing/isoZBilgerSample --isoSdSurface H2O --isoSdSurfaceValue 0.01 0.02 0.03 0.04 0.05 0.05335953 0.06 0.07 0.08 0.08488209 0.09 0.10 0.11 0.012 0.13  --keepDataOption 1 --averageWidth 32 --savgolFilterWidth 25


