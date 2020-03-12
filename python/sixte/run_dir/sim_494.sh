#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 494 
python ../pre-process-esass.py 494 
python ../esass.py 494 
