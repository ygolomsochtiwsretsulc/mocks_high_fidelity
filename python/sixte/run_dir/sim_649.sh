#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 649 
python ../pre-process-esass.py 649 
python ../esass.py 649 
