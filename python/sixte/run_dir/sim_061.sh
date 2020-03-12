#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 061 
python ../pre-process-esass.py 061 
python ../esass.py 061 
