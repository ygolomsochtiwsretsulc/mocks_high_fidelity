#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 020 
python ../pre-process-esass.py 020 
python ../esass.py 020 
