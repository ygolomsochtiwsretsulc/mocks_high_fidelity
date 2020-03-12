#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 572 
python ../pre-process-esass.py 572 
python ../esass.py 572 
