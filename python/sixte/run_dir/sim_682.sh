#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 682 
python ../pre-process-esass.py 682 
python ../esass.py 682 
