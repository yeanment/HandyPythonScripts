#!/bin/bash
# python CFDiffusionFlame/initialParaFuelDiluXAir.py \
#        --Fuel H2 --Diluent N2 --XFuel 0.25 --TransportModel Mix \
#        > data/CFDIgnFuelDiluXAirTMP/log.H2N2_XFF0.25_Mix 2>& 1
python CFDiffusionFlame/initialParaFuelDiluXAir.py \
       --Fuel H2 --Diluent HE --XFuel 0.25 --TransportModel Mix \
       > data/CFDIgnFuelDiluXAirTMP/log.H2HE_XFF0.25_Mix 2>& 1
# python CFDiffusionFlame/initialParaFuelDiluXAir.py \
#        --Fuel H2 --Diluent AR --XFuel 0.25 --TransportModel Mix \
#        > data/CFDIgnFuelDiluXAirTMP/log.H2AR_XFF0.25_Mix 2>& 1

# python CFDiffusionFlame/initialParaFuelDiluXAir.py \
#        --Fuel H2 --Diluent N2 --XFuel 0.25 --TransportModel Multi \
#        > data/CFDIgnFuelDiluXAirTMP/log.H2N2_XFF0.25_Multi 2>& 1
python CFDiffusionFlame/initialParaFuelDiluXAir.py \
       --Fuel H2 --Diluent HE --XFuel 0.25 --TransportModel Multi \
       > data/CFDIgnFuelDiluXAirTMP/log.H2HE_XFF0.25_Multi 2>& 1
# python CFDiffusionFlame/initialParaFuelDiluXAir.py \
#         --Fuel H2 --Diluent AR --XFuel 0.25 --TransportModel Multi \
#         > data/CFDIgnFuelDiluXAirTMP/log.H2AR_XFF0.25_Multi 2>& 1

# python CFDiffusionFlame/initialParaFuelDiluXAir.py \
#        --Fuel H2 --Diluent N2 --XFuel 0.25 --TransportModel Multi --EnableSoret True \
#        > data/CFDIgnFuelDiluXAirTMP/log.H2N2_XFF0.25_MultiSoret 2>& 1
python CFDiffusionFlame/initialParaFuelDiluXAir.py \
       --Fuel H2 --Diluent HE --XFuel 0.25 --TransportModel Multi --EnableSoret True \
       > data/CFDIgnFuelDiluXAirTMP/log.H2HE_XFF0.25_MultiSoret 2>& 1
# python CFDiffusionFlame/initialParaFuelDiluXAir.py \
#         --Fuel H2 --Diluent AR --XFuel 0.25 --TransportModel Multi --EnableSoret True \
#         > data/CFDIgnFuelDiluXAirTMP/log.H2AR_XFF0.25_MultiSoret 2>& 1

# python CFDiffusionFlame/initialParaFuelDiluXAir.py \
#        --Fuel H2 --Diluent N2 --XFuel 0.25 --TransportModel UnityLewis \
#        > data/CFDIgnFuelDiluXAirTMP/log.H2N2_XFF0.25_Unity 2>& 1
python CFDiffusionFlame/initialParaFuelDiluXAir.py \
       --Fuel H2 --Diluent HE --XFuel 0.25 --TransportModel UnityLewis \
       > data/CFDIgnFuelDiluXAirTMP/log.H2HE_XFF0.25_Unity 2>& 1
# python CFDiffusionFlame/initialParaFuelDiluXAir.py \
#        --Fuel H2 --Diluent AR --XFuel 0.25 --TransportModel UnityLewis \
#        > data/CFDIgnFuelDiluXAirTMP/log.H2AR_XFF0.25_Unity 2>& 1