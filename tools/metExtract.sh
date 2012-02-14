#!/bin/bash

#Purpose: cut out data fields from a *.met file
#

format="old"

if [ $# -lt 2 ]; then
	echo "You need at least two parameters: (dflt | ant | cal | opera) <filename>"
	exit
fi
headerFile=/home/fierz/usr/share/metFields.$1.txt
if [ $# -lt 3 ]; then
	head -21 ${headerFile}
	echo "Enter your choice as a field list for cut, e.g. 3,6-8,89:"
	read fields
else
	fields=$3
fi
if [ $# -eq 4 ]; then
	tz=".utc${4}"
else
	tz=""
fi
fieldNames=`tail -1 ${headerFile} | cut -d',' -f${fields}`
echo "Extract fields ${fieldNames} from file $2"
echo "date${tz}" > tmpH0
tail -1 ${headerFile} | cut -d',' -f${fields} | paste -d',' tmpH0 - > tmpHead

tail -n +17 $2 > tmp0
# Deal with date
cut -d',' -f2 tmp0 > tmpDate
if [ ${format} == "old" ]; then
	cut -d' ' -f1 tmpDate | cut -d'.' -f1 > tmpDay
	cut -d' ' -f1 tmpDate | cut -d'.' -f2 > tmpMonth
	cut -d' ' -f1 tmpDate | cut -d'.' -f3 > tmpYear
	cut -d' ' -f2 tmpDate > tmpTime
	paste -d'-' tmpYear tmpMonth tmpDay > tmpDate
	paste -d' ' tmpDate tmpTime > tmpDatime
fi
# Extract data fields
cut -d',' -f3- tmp0 > tmp2
cut -d',' -f${fields} tmp2 | paste -d',' tmpDatime - > tmpData
filetrunc=`ls ${2} | cut -d'.' -f1`
fileout=${filetrunc}.txt
cat tmpHead tmpData > ${fileout}

rm -f tmp*
