#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 016 
python ../pre-process-esass.py 016 
python ../esass.py 016 
