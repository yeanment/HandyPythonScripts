#!/bin/bash
python tecgetrfcsv_v3.py --timeStart 0  --case /mnt/e/CFDIgn-H2-Air/XH20.25-HE-Air/rhoU2_L20_R30_UO0.50_R200_T200_C0.0/S3.0E10/postProcessing --isoSdSurface H2O --isoSdSurfaceValue 0.02 0.04 0.06 0.08 0.10 --savgolFilterWidth 1 --enableSavgolFilter False --leftMaxPos False --writeTimeStepCombine True


python tecgetrfcsv_v3.py --timeStart 0  --case /mnt/e/CFDIgn-H2-Air/XH20.25-HE-Air/rhoU2_L20_R30_UO0.50_R200_T200_C0.0/S3.0E10/postProcessing --isoSdSurface H2O --isoSdSurfaceValue 0.02 0.04 0.06 0.08 0.10 --savgolFilterWidth 1 --enableSavgolFilter False --leftMaxPos True --writeTimeStepCombine True



# python tecgetrfcsv_v3.py --timeStart 0  --case /mnt/e/CFDIgn-H2-Air/XH20.25-N2-Air/rhoU2_L20_R30_UO5.00_R200_T200_C0.0_PLOT/S9.2509E09/postProcessing --isoSdSurface H2O --isoSdSurfaceValue 0.02 0.04 0.06 0.08 0.10 --savgolFilterWidth 1 --enableSavgolFilter False --leftMaxPos False --writeTimeStepCombine True

# python tecgetrfcsv_v3.py --timeStart 0  --case /mnt/e/CFDIgn-H2-Air/XH20.25-N2-Air/rhoU2_L20_R30_UO5.00_R200_T200_C0.0_PLOT/S9.2509E09/postProcessing --isoSdSurface H2O --isoSdSurfaceValue 0.02 0.04 0.06 0.08 0.10 --savgolFilterWidth 1 --enableSavgolFilter False --leftMaxPos True --writeTimeStepCombine True



