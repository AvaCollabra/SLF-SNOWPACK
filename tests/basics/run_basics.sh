#!/bin/bash
#This test assumes that run_res1exp has sucfcessfully completed, so it can analyze its results

TMP_REF="/tmp/run_basics_$$_ref"
TMP_NEW="/tmp/run_basics_$$_new"

function compare_result {
	param=$1
	prec=$2
	
	comp=`numdiff ${prec} ${TMP_REF} ${TMP_NEW} | grep "+++"`
	comp_result=`echo "${comp}" | grep -vE "are equal$"`
	printf "Comparing %-50s" "${param}"
	if [ -z "${comp_result}" ]; then
		printf "[OK]\n"
	else
		printf "[fail]\n"
	fi
}

#prepare the reference file
rm -f ../res1exp/output_ref/MST96_res.met
bunzip2 -k ../res1exp/output_ref/MST96_res.met.bz2

printf "*** basic checks:\n"
#check the snow height
../../tools/SnExtract.sh ../res1exp/output_ref/MST96_res.met 30 > ${TMP_REF}
../../tools/SnExtract.sh ../res1exp/output/MST96_res.met 30 > ${TMP_NEW}
compare_result "snow height (HS)" "-a 1."
#check the surface temperature
../../tools/SnExtract.sh ../res1exp/output_ref/MST96_res.met 13 > ${TMP_REF}
../../tools/SnExtract.sh ../res1exp/output/MST96_res.met 13 > ${TMP_NEW}
compare_result "surface temperature (TSS)" "-a 1."
#check the albedo
../../tools/SnExtract.sh ../res1exp/output_ref/MST96_res.met 11 > ${TMP_REF}
../../tools/SnExtract.sh ../res1exp/output/MST96_res.met 11 > ${TMP_NEW}
compare_result "albedo (ALB)" "-a 0.05"



printf "\n**** check the mass balance:\n"
#check the snow water equivalent
../../tools/SnExtract.sh ../res1exp/output_ref/MST96_res.met 36 > ${TMP_REF}
../../tools/SnExtract.sh ../res1exp/output/MST96_res.met 36 > ${TMP_NEW}
compare_result "snow water eq	uiv. (SWE)" "-a .5"
#check the snow rate
../../tools/SnExtract.sh ../res1exp/output_ref/MST96_res.met 29 > ${TMP_REF}
../../tools/SnExtract.sh ../res1exp/output/MST96_res.met 29 > ${TMP_NEW}
compare_result "snow rate" "-a .1"
#check the rain rate
../../tools/SnExtract.sh ../res1exp/output_ref/MST96_res.met 38 > ${TMP_REF}
../../tools/SnExtract.sh ../res1exp/output/MST96_res.met 38 > ${TMP_NEW}
compare_result "rain rate" "-a .1"
#check the snowpack runoff
../../tools/SnExtract.sh ../res1exp/output_ref/MST96_res.met 39 > ${TMP_REF}
../../tools/SnExtract.sh ../res1exp/output/MST96_res.met 39 > ${TMP_NEW}
compare_result "snowpack runoff" "-a .5"



printf "\n**** check the energy balance components:\n"
#check the sensible heat
../../tools/SnExtract.sh ../res1exp/output_ref/MST96_res.met 3 > ${TMP_REF}
../../tools/SnExtract.sh ../res1exp/output/MST96_res.met 3 > ${TMP_NEW}
compare_result "sensible heat" "-r 1e-2"
#check the latent heat
../../tools/SnExtract.sh ../res1exp/output_ref/MST96_res.met 4 > ${TMP_REF}
../../tools/SnExtract.sh ../res1exp/output/MST96_res.met 4 > ${TMP_NEW}
compare_result "latent heat" "-r 1e-2"
#check olwr
../../tools/SnExtract.sh ../res1exp/output_ref/MST96_res.met 5 > ${TMP_REF}
../../tools/SnExtract.sh ../res1exp/output/MST96_res.met 5 > ${TMP_NEW}
compare_result "OLWR" "-r 1e-2"
#check ilwr
../../tools/SnExtract.sh ../res1exp/output_ref/MST96_res.met 6 > ${TMP_REF}
../../tools/SnExtract.sh ../res1exp/output/MST96_res.met 6 > ${TMP_NEW}
compare_result "ILWR" "-r 1e-2"
#check rswr
../../tools/SnExtract.sh ../res1exp/output_ref/MST96_res.met 8 > ${TMP_REF}
../../tools/SnExtract.sh ../res1exp/output/MST96_res.met 8 > ${TMP_NEW}
compare_result "RSWR" "-r 1e-2"
#check iswr
../../tools/SnExtract.sh ../res1exp/output_ref/MST96_res.met 9 > ${TMP_REF}
../../tools/SnExtract.sh ../res1exp/output/MST96_res.met 9 > ${TMP_NEW}
compare_result "ISWR" "-r 1e-2"
#ground flux
../../tools/SnExtract.sh ../res1exp/output_ref/MST96_res.met 16 > ${TMP_REF}
../../tools/SnExtract.sh ../res1exp/output/MST96_res.met 16 > ${TMP_NEW}
compare_result "ground flux" "-r 1e-2"
#rain heat flux
../../tools/SnExtract.sh ../res1exp/output_ref/MST96_res.met 19 > ${TMP_REF}
../../tools/SnExtract.sh ../res1exp/output/MST96_res.met 19 > ${TMP_NEW}
compare_result "rain heat flux" "-r 1e-2"
#surface input heat flux
../../tools/SnExtract.sh ../res1exp/output_ref/MST96_res.met 97 > ${TMP_REF}
../../tools/SnExtract.sh ../res1exp/output/MST96_res.met 97 > ${TMP_NEW}
compare_result "surface input heat flux" "-r 1e-2"



printf "\n**** check the internal energy state:\n"
#internal energy change
../../tools/SnExtract.sh ../res1exp/output_ref/MST96_res.met 96 > ${TMP_REF}
../../tools/SnExtract.sh ../res1exp/output/MST96_res.met 96 > ${TMP_NEW}
compare_result "internal energy change" "-r 1e-2"
#phase change heat flux
../../tools/SnExtract.sh ../res1exp/output_ref/MST96_res.met 102 > ${TMP_REF}
../../tools/SnExtract.sh ../res1exp/output/MST96_res.met 102 > ${TMP_NEW}
compare_result "phase change heat flux" "-r 1e-2"
#check the liquid water content
../../tools/SnExtract.sh ../res1exp/output_ref/MST96_res.met 54 > ${TMP_REF}
../../tools/SnExtract.sh ../res1exp/output/MST96_res.met 54 > ${TMP_NEW}
compare_result "liquid water content (LWC)" "-a 1."

##Cleanup
rm -f ../res1exp/output_ref/MST96_res.met
rm -f ${TMP_REF}
rm -f ${TMP_NEW}
