#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 551 
python ../pre-process-esass.py 551 
python ../esass.py 551 
