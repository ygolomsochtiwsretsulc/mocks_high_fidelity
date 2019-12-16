#!/bin/bash

python2 S5_4FS_input_WI_format.py
gzip S5_Clusters_latest_input_cat_4FS_new_format.fits
python2 S5_LSM_create_function.py
python2 S5_SSM_create_function.py
python2 S5_tabulate_NZ.py
python2 S5_generate_sub_survey_params_table.py

