#!/bin/bash

# cat log.run2 | grep -E "^Time (\s)*=(\s)*([0-9]*(\.?)[0-9]*)" | sed 's/^Time\(\s\)*=\s*//g' > time.log
cat log.run2 | grep -E "^Time (\s)*=(\s)*" | sed 's/^Time\(\s\)*=\s*//g' > time.log
cat log.run2 | grep -E "T min/max(\s)*=(\s)*" | sed 's/T min\/max\s*=//g' | sed 's/\//,/g' > Tminmax.log
cat log.run2 | grep -E "rho min/max(\s)*=" | sed 's/rho min\/max\s*=//g' | sed 's/\//,/g' > rhominmax.log
cat log.run2 | grep -E "p min/max(\s)*=" | sed 's/p min\/max\s*=//g' | sed 's/\//,/g' > pminmax.log
cat log.run2 | grep -E "q min/max(\s)*=" | sed 's/q min\/max\s*=//g' | sed 's/\//,/g' > qminmax.log
cat log.run2 | grep -E "Flame([^(]*\()" | sed 's/Flame\([^(]*(\)\(\s\)*//g' | sed 's/).*//g' | sed 's/\s/,/g' > position.log

paste -d, time.log Tminmax.log rhominmax.log pminmax.log qminmax.log position.log> timevar.dat
# sed -i '1s/^/time, Tmin, Tmax, rhomin, rhomax, qmin, qmax, xf, yf, zf\n/' timevar.dat
rm time.log Tminmax.log rhominmax.log qminmax.log position.log pminmax.log