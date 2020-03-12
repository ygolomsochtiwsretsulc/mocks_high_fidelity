#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 648 
python ../pre-process-esass.py 648 
python ../esass.py 648 
