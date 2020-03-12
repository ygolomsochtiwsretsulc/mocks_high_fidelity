#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 466 
python ../pre-process-esass.py 466 
python ../esass.py 466 
