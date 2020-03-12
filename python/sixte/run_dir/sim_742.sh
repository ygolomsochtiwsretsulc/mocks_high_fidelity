#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 742 
python ../pre-process-esass.py 742 
python ../esass.py 742 
