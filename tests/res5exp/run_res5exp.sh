#!/bin/bash
../../bin/snowpack -c io_res5exp.ini -e 1996-06-17T00:00

PREC="1e-3"
#north slopes
numdiff -r ${PREC} output_ref/MST961_tst.haz output/MST961_tst.haz | grep "+++"
numdiff -r ${PREC} output_ref/MST961_tst.sno output/MST961_tst.sno | grep "+++"

#south slopes
numdiff -r ${PREC} output_ref/MST963_tst.haz output/MST963_tst.haz | grep "+++"
numdiff -r ${PREC} output_ref/MST963_tst.sno output/MST963_tst.sno | grep "+++"

#north slopes
rm -f output_ref/MST961_tst.met
bunzip2 -k output_ref/MST961_tst.met.bz2
sed -i '11d' output_ref/MST961_tst.met; sed -i '11d' output/MST961_tst.met
numdiff -r ${PREC} --speed-large-files output_ref/MST961_tst.met output/MST961_tst.met | grep "+++"
rm -f output_ref/MST961_tst.met

rm -f output_ref/MST961_tst.pro
bunzip2 -k output_ref/MST961_tst.pro.bz2
sed -i '10d' output_ref/MST961_tst.pro; sed -i '10d' output/MST961_tst.pro
numdiff -r ${PREC} --speed-large-files output_ref/MST961_tst.pro output/MST961_tst.pro | grep "+++"
rm -f output_ref/MST961_tst.pro

#south slopes
rm -f output_ref/MST963_tst.met
bunzip2 -k output_ref/MST963_tst.met.bz2
sed -i '11d' output_ref/MST963_tst.met; sed -i '11d' output/MST963_tst.met
numdiff -r ${PREC} --speed-large-files output_ref/MST963_tst.met output/MST963_tst.met | grep "+++"
rm -f output_ref/MST963_tst.met

rm -f output_ref/MST963_tst.pro
bunzip2 -k output_ref/MST963_tst.pro.bz2
sed -i '10d' output_ref/MST963_tst.pro; sed -i '10d' output/MST963_tst.pro
numdiff -r ${PREC} --speed-large-files output_ref/MST963_tst.pro output/MST963_tst.pro | grep "+++"
rm -f output_ref/MST963_tst.pro

