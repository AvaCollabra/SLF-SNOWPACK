#!/bin/bash

# use  nohup *.sh > fout  to redirect output to named file

TOOL="valgrind --tool=callgrind --simulate-cache=yes --dump-instr=yes --trace-jumps=yes"
TOOL="/software/bin/valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes --suppressions=../../tools/valgrindOCCI.supp --log-fd=2"
TOOL="time"
TOOL=""

cp -p init/DAV5* current_snow
rm -f current_snow/DAV5.haz
${TOOL} ./snowpack -c cfgfiles/io_oper.ini -m operational -s DAV5 -e NOW
