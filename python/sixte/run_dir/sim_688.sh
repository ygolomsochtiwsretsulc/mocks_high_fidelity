#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 688 
python ../pre-process-esass.py 688 
python ../esass.py 688 
