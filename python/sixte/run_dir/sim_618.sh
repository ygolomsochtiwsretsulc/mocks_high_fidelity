#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 618 
python ../pre-process-esass.py 618 
python ../esass.py 618 
