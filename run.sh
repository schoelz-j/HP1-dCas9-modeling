#!/bin/bash

mamba activate hp1_modeling

# Run this code from the directory hp1_modeling

# 1. HP1 only model. Process data, train model, create plots Figure 1,
# and Figure S1. Run simulations for Figures 6-8.
R code/models_and_predictions/HP1-only_Models_DataPrep.R
python code/models_and_predictions/HP1-only-models.py
R code/plots/Figure_1.R
R code/plots/Figure_S1.R

# 2. HP1 plus and No HP1 models.

# 3. Simulations.
