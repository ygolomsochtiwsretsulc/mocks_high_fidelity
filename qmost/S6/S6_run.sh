#!/bin/bash

python2 S6_4FS_input_WI_format.py
rm eROSITA_component_4FS_new_format.fits
rm IR_component_4FS_new_format.fits
rm S6_summary.fits.gz
gzip S6_summary.fits
python2 S6_LSM_create_function.py
python2 S6_SSM_create_function.py
python2 S6_tabulate_NZ.py
python2 S6_generate_sub_survey_params_table.py

