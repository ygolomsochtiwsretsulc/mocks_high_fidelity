#!/bin/bash

python2 S8_4FS_input_WI_format.py
gzip S8_Cosmology_latest_input_cat_4FS_new_format.fits
python2 S8_LSM_create_function.py
python2 S8_SSM_create_function.py
python2 S8_tabulate_NZ.py
python2 S8_generate_sub_survey_params_table.py

