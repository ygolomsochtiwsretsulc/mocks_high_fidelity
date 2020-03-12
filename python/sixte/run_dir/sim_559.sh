#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 559 
python ../pre-process-esass.py 559 
python ../esass.py 559 
