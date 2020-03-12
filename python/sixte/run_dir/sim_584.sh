#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 584 
python ../pre-process-esass.py 584 
python ../esass.py 584 
