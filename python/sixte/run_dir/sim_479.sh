#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 479 
python ../pre-process-esass.py 479 
python ../esass.py 479 
