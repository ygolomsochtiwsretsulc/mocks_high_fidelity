#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 611 
python ../pre-process-esass.py 611 
python ../esass.py 611 
