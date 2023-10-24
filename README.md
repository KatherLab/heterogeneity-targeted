# Heterogeneity-targeted

"Heterogeneity-targeted" is a Python-based tool designed for computing heterogeneity metrics, specifically utilizing the outputs (mutation classification: MUT or WT) Marugoto tile-based method to calculate Shannon entropy. The tool extends its functionality by conducting survival analyses, particularly employing COX proportional-hazards (COX PH) analysis, and producing Kaplan-Meier curves. Moreover, it facilitates the generation of composite statuses by combining binary genetic alteration and heterogeneity statuses, leading to the formation of four distinct status combinations. Subsequently, COX PH analysis is applied to each combination, followed by the visualization of Kaplan-Meier curves for the resultant four groups.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Requirements](#requirements)
- [Features](#features)

## Installation

Make sure you have Python version 3.9 or higher installed.


    pip install .


## Usage

    python -m heterogeneity-targeted --target genetic_alteration_name \
                                     --path_tile_preds /absolute/path/to/patient_tile_prediction_csv \
                                     --path_patient_pred /absolute/path/to/patient_patient_prediction_csv \
                                     --path_clini_table /path/to/clini_table \
                                     --output_directory /absolute/path/to/desired_output_dir_name \
                                     --duration_col Fuyears \
                                     --event_col crc_death
             

## Requirements

- Python >= 3.9
- argparse=1.4.0
- lifelines=0.27.8
- matplotlib=3.8.0"
- numpy=1.26.1
- pandas~=2.1
- pathlib=1.0.1   


## Features
1) Shannon Entropy Calculation:

    Utilizes a Marugoto tile-based approach for heterogeneity metric calculation.

2) COX PH Analysis:

    Performs survival analyses using COX proportional-hazards analysis.

3) Kaplan-Meier Curves:

    Plots Kaplan-Meier curves for the calculated metrics.

4) Combined Status Analysis:

    Creates combined forms of genetic alteration status and heterogeneity status, resulting in four different combinations.

5) COX PH Analysis for Combined Status:

    Applies COX PH analysis to each combination of combined statuses.

6) Kaplan-Meier Curves for Combined Status:

    Generates Kaplan-Meier curves for the four resulting groups.
