#! /bin/bash
# conda activate ofpost
python tecisoZBilger_v2a.py --timeStart 0. --timeEnd -1  --writeTimeStepCombine True --case /mnt/e/CFDiff-H2-Air/XH20.50-N2-Air/rhoU2_UO0.50_L20_R30_R0.2_T0.2_C0.0/S3.0E10/postProcessing/sampleDict/ --isoSdSurface T --isoSdSurfaceValue  800 1000 1600  --keepDataOption 1 --averageWidth 11 --savgolFilterWidth 25
python tecisoZBilger_v2a.py --timeStart 0. --timeEnd -1  --writeTimeStepCombine True --case /mnt/e/CFDiff-H2-Air/XH20.50-N2-Air/rhoU2_UO0.50_L20_R30_R0.2_T0.2_C0.0/S3.0E10/postProcessing/sampleDict/ --isoSdSurface T --isoSdSurfaceValue  800 1000 1600  --keepDataOption 1 --averageWidth 25 --savgolFilterWidth 25
