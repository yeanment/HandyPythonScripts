#!/bin/bash
# python tecisoZBilger_v1.py --timeStart 0.00 --timeEnd -1  --writeTimeStepCombine False --case  /BIGDATA2/pku_zchen_1/xsm/OpenFOAM/xsm-v1712/run/CFDiff-H2Air/CFDH2-L_15mm-ag_100.00-phiS_0.50-UFUO-R_0.2mm-T_0.2ms-W10-E0.220494539mJ/postProcessing/isoZBilgerSample --isoSdSurface T --isoSdSurfaceValue 400 500 600 700 800 900 1000 1100 1200 1300 1400 --isoZBilgerValue 0.6666666666 --rfPosMaximal False --keepDataOption 1

python tecisoZBilger_v1.py --timeStart 0.00 --timeEnd -1  --writeTimeStepCombine False --case  /BIGDATA2/pku_zchen_1/xsm/OpenFOAM/xsm-v1712/run/CFDiff-H2Air/CFDH2-L_15mm-ag_100.00-phiS_0.50-UFUO-R_0.2mm-T_0.2ms-W10-E0.220494539mJ/postProcessing/isoZBilgerSample --isoSdSurface H2O --isoSdSurfaceValue 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 --isoZBilgerValue 0.6666666666 --rfPosMaximal False --keepDataOption 1

