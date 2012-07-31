#!/bin/bash
TESTDELTA=10
bash massbalancecheck.sh ../res1exp/output/MST96_res.met 2>&1 >/dev/null | sed -r 's/.*(ERROR).*/\1/' | awk '(NR==1) {print ($NF>'${TESTDELTA}' || NF==1)?"ERROR":0}'
