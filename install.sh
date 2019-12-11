#!/bin/bash

# run model installation scripts
python install_all_models_step1.py
python install_control_model_step2.py
bash install_control_model_step3.sh
python install_control_model_step4.py
python install_reduced_model.py
python install_manipulation1_model.py
python install_manipulation2_model.py

# run analysis installations scripts
python install_control_model_analysis.py
python install_reduced_model_analysis.py
python install_manipulation1_model_analysis.py
python install_manipulation2_model_analysis.py
