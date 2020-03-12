#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 717 
python ../pre-process-esass.py 717 
python ../esass.py 717 
