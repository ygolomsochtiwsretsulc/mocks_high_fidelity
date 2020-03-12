#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 530 
python ../pre-process-esass.py 530 
python ../esass.py 530 
