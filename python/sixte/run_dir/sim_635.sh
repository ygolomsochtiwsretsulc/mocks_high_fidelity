#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 635 
python ../pre-process-esass.py 635 
python ../esass.py 635 
