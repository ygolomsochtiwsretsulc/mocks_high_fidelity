#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 711 
python ../pre-process-esass.py 711 
python ../esass.py 711 
