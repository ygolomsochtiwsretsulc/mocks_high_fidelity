#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 664 
python ../pre-process-esass.py 664 
python ../esass.py 664 
