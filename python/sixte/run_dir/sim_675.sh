#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 675 
python ../pre-process-esass.py 675 
python ../esass.py 675 
