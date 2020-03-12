#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 609 
python ../pre-process-esass.py 609 
python ../esass.py 609 
