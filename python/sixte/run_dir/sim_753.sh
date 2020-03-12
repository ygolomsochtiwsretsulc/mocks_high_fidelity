#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 753 
python ../pre-process-esass.py 753 
python ../esass.py 753 
