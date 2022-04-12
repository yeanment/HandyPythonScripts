#! /bin/bash
# conda activate ofpost
python tecisoZBilger_v2a.py --timeStart 0. --timeEnd -1  --writeTimeStepCombine True --case /mnt/e/CFDIgn-H2-Air/XH20.25-N2-Air/rhoU2_L20_R30_UO5.00_R200_T200_C0.0_PLOT/S1.01461E10/postProcessing/sampleDict/ --isoSdSurface H2O --isoSdSurfaceValue 0.02 0.04 0.06 0.08 0.1 --keepDataOption 1 --averageWidth 11 --savgolFilterWidth 25
python tecisoZBilger_v2a.py --timeStart 0. --timeEnd -1  --writeTimeStepCombine True --case /mnt/e/CFDIgn-H2-Air/XH20.25-N2-Air/rhoU2_L20_R30_UO5.00_R200_T200_C0.0_PLOT/S9.2509E09/postProcessing/sampleDict/ --isoSdSurface H2O --isoSdSurfaceValue 0.02 0.04 0.06 0.08 0.1 --keepDataOption 1 --averageWidth 11 --savgolFilterWidth 25
python tecisoZBilger_v2a.py --timeStart 0. --timeEnd -1  --writeTimeStepCombine True --case /mnt/e/CFDIgn-H2-Air/XH20.25-N2-Air/rhoU2_L20_R30_UO5.00_R200_T200_C0.0_PLOT/S9.5493E09/postProcessing/sampleDict/ --isoSdSurface H2O --isoSdSurfaceValue 0.02 0.04 0.06 0.08 0.1 --keepDataOption 1 --averageWidth 11 --savgolFilterWidth 25

# python tecisoZBilger_v2a.py --timeStart 0. --timeEnd -1  --writeTimeStepCombine True --case /mnt/e/CFDiff-H2-Air/XH20.50-N2-Air/rhoU2_UO0.50_L20_R30_R0.2_T0.2_C0.0/S3.0E10/postProcessing/sampleDict/ --isoSdSurface T --isoSdSurfaceValue  800 1000 1600  --keepDataOption 1 --averageWidth 25 --savgolFilterWidth 25
