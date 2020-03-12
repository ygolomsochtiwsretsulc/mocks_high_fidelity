#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 672 
python ../pre-process-esass.py 672 
python ../esass.py 672 
