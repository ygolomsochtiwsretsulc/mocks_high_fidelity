#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 586 
python ../pre-process-esass.py 586 
python ../esass.py 586 
