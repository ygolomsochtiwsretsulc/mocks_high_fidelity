#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 622 
python ../pre-process-esass.py 622 
python ../esass.py 622 
