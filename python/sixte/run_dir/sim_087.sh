#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 087 
python ../pre-process-esass.py 087 
python ../esass.py 087 
