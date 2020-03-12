#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 638 
python ../pre-process-esass.py 638 
python ../esass.py 638 
