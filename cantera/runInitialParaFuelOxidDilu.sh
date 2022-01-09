#!/bin/bash
# conda activate cantera

python edgeflame/initialParaFuelOxidDilu.py --Diluent AR --DilutionRatio 7.5 --Zst 0.1 --TransportModel UnityLewis >> data/EdgeFlameFuelOxidDiluZst/H2-O2-AR-qDilu7.5-Zst0.1.log
python edgeflame/initialParaFuelOxidDilu.py --Diluent AR --DilutionRatio 7.5 --Zst 0.2 --TransportModel UnityLewis >> data/EdgeFlameFuelOxidDiluZst/H2-O2-AR-qDilu7.5-Zst0.2.log
python edgeflame/initialParaFuelOxidDilu.py --Diluent AR --DilutionRatio 7.5 --Zst 0.3 --TransportModel UnityLewis >> data/EdgeFlameFuelOxidDiluZst/H2-O2-AR-qDilu7.5-Zst0.3.log
python edgeflame/initialParaFuelOxidDilu.py --Diluent AR --DilutionRatio 7.5 --Zst 0.4 --TransportModel UnityLewis >> data/EdgeFlameFuelOxidDiluZst/H2-O2-AR-qDilu7.5-Zst0.4.log
python edgeflame/initialParaFuelOxidDilu.py --Diluent AR --DilutionRatio 7.5 --Zst 0.5 --TransportModel UnityLewis >> data/EdgeFlameFuelOxidDiluZst/H2-O2-AR-qDilu7.5-Zst0.5.log
python edgeflame/initialParaFuelOxidDilu.py --Diluent AR --DilutionRatio 7.5 --Zst 0.6 --TransportModel UnityLewis >> data/EdgeFlameFuelOxidDiluZst/H2-O2-AR-qDilu7.5-Zst0.6.log
python edgeflame/initialParaFuelOxidDilu.py --Diluent AR --DilutionRatio 7.5 --Zst 0.7 --TransportModel UnityLewis >> data/EdgeFlameFuelOxidDiluZst/H2-O2-AR-qDilu7.5-Zst0.7.log
python edgeflame/initialParaFuelOxidDilu.py --Diluent AR --DilutionRatio 7.5 --Zst 0.8 --TransportModel UnityLewis >> data/EdgeFlameFuelOxidDiluZst/H2-O2-AR-qDilu7.5-Zst0.8.log
python edgeflame/initialParaFuelOxidDilu.py --Diluent AR --DilutionRatio 7.5 --Zst 0.9 --TransportModel UnityLewis >> data/EdgeFlameFuelOxidDiluZst/H2-O2-AR-qDilu7.5-Zst0.9.log
