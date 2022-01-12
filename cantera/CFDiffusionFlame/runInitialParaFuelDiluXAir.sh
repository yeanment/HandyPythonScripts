#!/bin/bash
python CFDiffusionFlame/initialParaFuelDiluXAir.py \
       --Diluent N2 --XFuel 0.25 --TransportModel Mix \
       > data/EdgeFlameFuelDiluXAir/log.H2N2_XFF0.25_Mix 2>& 1
python CFDiffusionFlame/initialParaFuelDiluXAir.py \
       --Fuel H2 --Diluent HE --XFuel 0.25 --TransportModel Mix \
       > data/EdgeFlameFuelDiluXAir/log.H2HE_XFF0.25_Mix 2>& 1
python CFDiffusionFlame/initialParaFuelDiluXAir.py \
       --Fuel H2 --Diluent AR --XFuel 0.25 --TransportModel Mix \
       > data/EdgeFlameFuelDiluXAir/log.H2AR_XFF0.25_Mix 2>& 1
python CFDiffusionFlame/initialParaFuelDiluXAir.py \
       --Fuel H2 --Diluent N2 --XFuel 0.25 --TransportModel Multi \
       > data/EdgeFlameFuelDiluXAir/log.H2N2_XFF0.25_Multi 2>& 1
python CFDiffusionFlame/initialParaFuelDiluXAir.py \
       --Fuel H2 --Diluent HE --XFuel 0.25 --TransportModel Multi \
       > data/EdgeFlameFuelDiluXAir/log.H2HE_XFF0.25_Multi 2>& 1
python CFDiffusionFlame/initialParaFuelDiluXAir.py \
        --Fuel H2 --Diluent AR --XFuel 0.25 --TransportModel Multi \
        > data/EdgeFlameFuelDiluXAir/log.H2AR_XFF0.25_Multi 2>& 1

