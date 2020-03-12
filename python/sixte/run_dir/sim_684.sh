#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 684 
python ../pre-process-esass.py 684 
python ../esass.py 684 
