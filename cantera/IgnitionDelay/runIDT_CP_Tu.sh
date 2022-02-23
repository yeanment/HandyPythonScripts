#!/bin/bash
python IgnitionDelay/IDT_CP_Tu.py \
       --Mech chem/DME_LuChen2015/DME_LuChen2015.cti \
       --Fuel CH3OCH3 --EquivalenceRatio 1.0 \
       > data/IgnitionDelay/log.IDT_DMEAir_phi1.0 2>& 1