#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 723 
python ../pre-process-esass.py 723 
python ../esass.py 723 
