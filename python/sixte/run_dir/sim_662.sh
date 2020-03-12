#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 662 
python ../pre-process-esass.py 662 
python ../esass.py 662 
