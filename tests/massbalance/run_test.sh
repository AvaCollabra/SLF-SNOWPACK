#!/bin/bash
testvalue=0.001
#Do mass balance check
bash massbalancecheck.sh ../res1exp/output/MST96_res.met 2>&1 >/dev/null | \
#If word ERROR is detected, make it first word on the line
sed -r 's/(.*)(ERROR)(.*)/\2 \1\2\3/' | \
#Now check if the limits are not exceeded. If so, write out ERROR in the end.
awk '{if(NR==1) {if(substr($1,1,5)=="ERROR") print "ERROR"}; if((NR==3 && $NF*$NF>'${testvalue}'*'${testvalue}') || (NR==4 && $NF*$NF>'${testvalue}'*'${testvalue}')) {error=1}; print $0} END {if(error==1) {print "ERROR"}}'
