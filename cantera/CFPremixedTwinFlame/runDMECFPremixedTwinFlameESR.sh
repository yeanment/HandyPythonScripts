#!/bin/bash
python CFPremixedTwinFlame/CFPremixedTwinFlameESR.py \
       --Mech chem/DME_LuChen2015/DME_LuChen2015.xml \
       --Fuel CH3OCH3 --EquivalenceRatio 1.0 --TransportModel Mix \
       --Tin 300 \
       > data/CFPremixedTwinFlame/log.DMEAir_phi1.0_Mix_Tu300 2>& 1
python CFPremixedTwinFlame/CFPremixedTwinFlameESR.py \
       --Mech chem/DME_LuChen2015/DME_LuChen2015.xml \
       --Fuel CH3OCH3 --EquivalenceRatio 1.0 --TransportModel Mix \
       --Tin 600 \
       > data/CFPremixedTwinFlame/log.DMEAir_phi1.0_Mix_Tu600 2>& 1
