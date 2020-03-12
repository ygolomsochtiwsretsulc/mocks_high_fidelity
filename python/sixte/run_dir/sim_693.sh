#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 693 
python ../pre-process-esass.py 693 
python ../esass.py 693 
