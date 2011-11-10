#!/bin/bash

# use  nohup *.sh > fout  to redirect output to named file

TOOL="valgrind --tool=callgrind --simulate-cache=yes --dump-instr=yes --trace-jumps=yes"
TOOL="/software/bin/valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes --suppressions=../../tools/valgrindOCCI.supp --log-fd=2"
TOOL="time"
TOOL=""

cp -p init/DAV2* current_snow
${TOOL} snowpack -c cfgfiles/io_resDAV2.ini -e 2010-06-01T00:00
