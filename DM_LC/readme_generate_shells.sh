#!/bin/bash

# computes the gemoetry
python 000_geometry.py 1000.0 MD10

# writes the commands to execute
python 001_process_hlists_WRITE_commands.py

# executes the commands
nohup sh md10_r000_all.sh > md10_r000_all.log & 

