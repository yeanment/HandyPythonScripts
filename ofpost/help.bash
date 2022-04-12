cat run1.log | grep -E "T min/max(\s)*=(\s)*([0-9]*(\.?)[0-9]*)/([0-9]*(\.?)[0-9]*)" | sed 's/T min\/max\s*=//g' | sed 's/\//,/g' > out.log

paste -d, time.csv out.csv > ot.csv


ls H2Ignition_768 | grep -E -v "0|processor" | xargs -i cp -rp ./H2Ignition_768/{} ./H2Ignition_768Plot/


cat log.run2 | grep -E "T min/max(\s)*=(\s)*([0-9]*(\.?)[0-9]*)/([0-9]*(\.?)[0-9]*)" | sed 's/T min\/max\s*=//g' | sed 's/\//,/g' > Tminmax.log

cat log.run2 | grep -E "T min/max(\s)*=(\s)*([0-9]*(\.?)[0-9]*)/([0-9]*(\.?)[0-9]*)" | sed 's/T min\/max\s*=//g' | sed 's/\//,/g' > Tminmax.log


cat log.run2 | grep -E "Flame([^(]*\()(\s)*([0-9]*(\.?)[0-9]*)(\s)*([0-9]*(\.?)[0-9]*)(\s)*([0-9]*(\.?)[0-9]*)(\s)*\).*" | sed 's/Flame\([^(]*(\)\(\s\)*//g' | sed 's/).*//g' | sed 's/\s/,/g' > position.log

ls ../H2Ignition_768 | grep -E "^0\.0(.{1,3})$" | xargs -i cp -rp ../H2Ignition_768/{} .
find . -type f -regex regEXP -exec rm -rf {} \;
find . -maxdepth 2 -type d -regextype posix-egrep -regex ".*/0\.0.{3,9}$" -exec rm -r {} \;

ls processor0 | grep -E -v "^0\.01$|constant|^0\.0(.{4,6})$" | xargs -i reconstructPar -time {}

mapFields -case . -consistent -sourceTime 0.1 -mapMethod cellPointInterpolate ../../Init_noreaction/H2Ig_V0.25_L20.0


ls | grep -E -v "pro.*|postP.*|*.out|log.*|*.dat"

tail -n 100 log.run | grep -E "Time|Flame|delta|min"


qmgr -c "set queue batch resources_default.walltime=3600"

ls | grep -E  "proc.*|postP.*|0.00.*|log.run.*" | xargs -i mv ./{} ./S8.6125E10/

find ./ -regex '.*lineY2.*' -type f |xargs -i rm -rf {}  &

pigz -dc archive.tar.gz | tar xf -

sed -i 's/^    symmetryPlane/    symPlane/g'  0/*

sed -i 's/libcanteraFlamePosDebug/libcanteraFlamePos/g'  system/controlDict*
sed -i 's/libs ("# liboptimizedChemistry_DME.so");/\/\/ libs ("liboptimizedChemistry_DME.so");/g'  system/controlDict

sed -i 's/ebiReacFoamDLBZIgnExp6XSM0/ebiReacFoamDLBZLeIgnExp6XSM0/g'  system/controlDict*


sed -i 's/writeInterval 1e-5;/writeInterval 2e-5;/g'  system/controlDict*
sed -i 's/\.\*sdT/"\.\*sdT"/g'  system/controlDict*

sed -i 's/value           700;/value           400;/g'  system/controlDict*

find . -type f -exec  sed -i 's/\t/    /g' {} +

  
yhbatch -n 24 -N  1 -p bigdata -J xsmReco ./

find . -mindepth 2 -maxdepth 2 -type d -regextype posix-egrep -regex  ".*/processor[0-9]{1,2}/0.00[0-9]{1,1}[0-46-9]$" -exec rm -r {} \; &
 
find . -mindepth 2 -maxdepth 2 -type d -regextype posix-egrep -regex  ".*/processor[0]{1,2}/0\.00[0-4]{1,1}[0-9]{1,9}$" -exec rm -r {} \; &
find . -mindepth 1 -maxdepth 1 -type d -regextype posix-egrep -regex  ".*/0\.01[0-9]{1,1}[1357]{1,9}$" -exec rm -r {} \;
find . -mindepth 3 -maxdepth 3 -type d -regextype posix-egrep -regex  ".*/5e-05$" -exec rm -r {} \;
find . -mindepth 3 -maxdepth 3 -type d -regextype posix-egrep -regex  ".*/0\.[0-9]{4,4}[5]{1,1}$" -exec rm -r {} \;

find . -mindepth 2 -maxdepth 5 -type f -regextype posix-egrep -regex ".*cellCpu.*" > log.find -exec rm -r {} \; &
find . -mindepth 2 -maxdepth 5 -type f -regextype posix-egrep -regex ".*referenceMap.*" > log.find -exec rm -r {} \; &
find . -mindepth 2 -maxdepth 5 -type f -regextype posix-egrep -regex ".*zMix.*" > log.find -exec rm -r {} \; &
find . -mindepth 2 -maxdepth 5 -type f -regextype posix-egrep -regex ".*ZBilger.*" > log.find -exec rm -r {} \; &
find . -mindepth 2 -maxdepth 3 -type f -regextype posix-egrep -regex ".*_sdT.*" > log.find -exec rm -r {} \; &
find . -mindepth 2 -maxdepth 3 -type f -regextype posix-egrep -regex ".*_sdH2O.*" > log.find -exec rm -r {} \; &

ls | grep -E "^0\.[0-9][0-46-9]$" | xargs -i rm -r ./{}

find -type d -empty



rm -r processor0* processor1* &
rm -r processor2* processor3* &
rm -r processor4* processor5* &
rm -r processor6* processor7* &
rm -r processor8* processor9* &


rm -r processor0* & 
rm -r processor1* &
rm -r processor2* & 
rm -r processor3* &
rm -r processor4* & 
rm -r processor5* &
rm -r processor6* & 
rm -r processor7* &
rm -r processor8* &
rm -r processor9* &

tar zcvfp lineSample.tgz ./lineSample >> /dev/null &
tar zcvfp lineSampleSdH2O.tgz ./lineSampleSdH2O >> /dev/null &
tar zcvfp lineSampleSdT.tgz ./lineSampleSdT/ >> /dev/null &

ls | grep "lineY.-iso" | xargs -i tar zcvfp {}.tgz ./{} >> /dev/null &

rm -r lineSample lineSampleSdH2O lineSampleSdT &


ls | grep -E -v "S.*" | xargs -i cp -r ./{} ../../

0.220494539


ln -s ~/xsm/OpenFOAM/localOpenFOAM5.x localOpenFOAM5.x



postProcess -funcs "(sampleDictIsoT sampleDictIsoH2O)" -time 0.00000001: -fields \
 "(p U Qdot T rho H O OH HO2 H2O2 H2 O2 H2O 
  Ks_sdT W_sdT cur_sdT grad_sdT igradU_sdT n_sdT nngradU_sdT 
  sd_conv_sdT sd_corr_sdT sd_diff_sdT sd_rr_sdT sd_sdT sd_unsteady_sdT 
  Ks_sdH2O W_sdH2O cur_sdH2O grad_sdH2O igradU_sdH2O n_sdH2O nngradU_sdH2O 
  sd_conv_sdH2O sd_corr_sdH2O sd_diff_sdH2O sd_rr_sdH2O sd_sdH2O sd_unsteady_sdH2O)"


ebiReacFoamDLBIgnExp6XSM0 -postProcess  -time '0.000001:' -dict ../controlDict.postEnergyBudget  
postProcess -dict ../controlDict.postSampleIsoSurface -time 0.00000001: -fields \
"(p U Qdot T rho H O OH HO2 H2O2 H2 O2 H2O N2 rho
 Ks_sdT W_sdT cur_sdT grad_sdT igradU_sdT n_sdT nngradU_sdT
  sd_conv_sdT sd_corr_sdT sd_diff_sdT sd_rr_sdT sd_sdT sd_unsteady_sdT
   Ks_sdH2O W_sdH2O cur_sdH2O grad_sdH2O igradU_sdH2O n_sdH2O nngradU_sdH2O
    sd_conv_sdH2O sd_corr_sdH2O sd_diff_sdH2O sd_rr_sdH2O sd_sdH2O sd_unsteady_sdH2O)"

 
 
ls |  awk '{if (strtonum($0) > 0.008000001) print($0);}'  | xargs -i rm -r ./{} 


git filter-branch --commit-filter 'if [ "$GIT_COMMITTER_EMAIL" = "yeanment@outlook.com" ];
  then git commit-tree -S "$@";
  else git commit-tree "$@";
  fi' HEAD



gcc alternative
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 7
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-7 7
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 9
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-9 9

sudo update-alternatives --config gcc
sudo update-alternatives --config g++

sudo update-alternatives --set gcc /usr/bin/gcc-9
sudo update-alternatives --set g++ /usr/bin/g++-9




bundle install
bundle exec jekyll serve
