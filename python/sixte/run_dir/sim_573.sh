#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 573 
python ../pre-process-esass.py 573 
python ../esass.py 573 
