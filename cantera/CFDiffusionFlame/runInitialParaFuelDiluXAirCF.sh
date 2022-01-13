#!/bin/bash
python CFDiffusionFlame/initialParaFuelDiluXAirCF.py \
       --Diluent N2 --XFuel 0.25 --TransportModel Mix \
       > data/EdgeFlameFuelDiluXAir/log.H2N2_XFF0.25_MixCF 2>& 1
# python CFDiffusionFlame/initialParaFuelDiluXAirCF.py \
#        --Fuel H2 --Diluent HE --XFuel 0.25 --TransportModel Mix \
#        > data/EdgeFlameFuelDiluXAir/log.H2HE_XFF0.25_MixCF 2>& 1
# python CFDiffusionFlame/initialParaFuelDiluXAirCF.py \
#        --Fuel H2 --Diluent AR --XFuel 0.25 --TransportModel Mix \
#        > data/EdgeFlameFuelDiluXAir/log.H2AR_XFF0.25_MixCF 2>& 1
# python CFDiffusionFlame/initialParaFuelDiluXAirCF.py \
#        --Fuel H2 --Diluent N2 --XFuel 0.25 --TransportModel Multi \
#        > data/EdgeFlameFuelDiluXAir/log.H2N2_XFF0.25_MultiCF 2>& 1
# python CFDiffusionFlame/initialParaFuelDiluXAirCF.py \
#        --Fuel H2 --Diluent HE --XFuel 0.25 --TransportModel Multi \
#        > data/EdgeFlameFuelDiluXAir/log.H2HE_XFF0.25_MultiCF 2>& 1
# python CFDiffusionFlame/initialParaFuelDiluXAirCF.py \
#         --Fuel H2 --Diluent AR --XFuel 0.25 --TransportModel Multi \
#         > data/EdgeFlameFuelDiluXAir/log.H2AR_XFF0.25_MultiCF 2>& 1

