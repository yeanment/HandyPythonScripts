#!/bin/bash
python CFPremixedTwinFlame/CFPremixedTwinFlameESR.py \
       --Mech chem/DME_LuChen2015/DME_LuChen2015.xml \
       --Fuel CH3OCH3 --EquivalenceRatio 1.0 --TransportModel Mix \
       > data/CFPremixedTwinFlame/log.DMEAir_phi1.0_Mix 2>& 1
