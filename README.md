
This repository contains Matlab code and required files to generate design alternatives and to calculate their robustness metrics for the robust design of five RoR hydropower plants.
The methodology presented in the paper by  V. Yildiz, S.F. Brown, C. Rougé  "Revisiting small hydropower design in a drought-prone world" submitted to the Water Resources Research.


# MATLAB Library: Amalgam_HYPER and Robustness_Analysis

This MATLAB library consists of two main folders: Amalgam_Hyper for multi-objective optimization (MO) and Robustness_Analysis for robustness analysis under different scenarios. Each folder contains separate setups for five case studies; 

Case Study 1: Besik RoR
Case Study 2: Tepe RoR
Case Study 3: Karacay RoR
Case Study 4: Kaplan RoR
Case Study 5: Buyukdere RoR

Additionally, a MATLAB code is provided for post-processing the robustness results once the analysis is completed.

## Amalgam_HYPER

The Amalgam_HYPER folder contains the multi-objective optimization (MO) scripts and setup for the five case studies. These experimental setups are designed to find optimal solutions by simultaneously considering multiple objectives. Each case study has its own setup along with project's inputs within the folder.

To use the Amalgam_HYPER code:

1. Ensure that you have MATLAB installed on your system.
2. Clone or download this repository to your local machine.
3. Open MATLAB and navigate to the Amalgam_HYPER folder.
4. The required input files, including streamflow records as well as project-based parameters, are already provided in the respective case study folders.
5. Adjust the population size and number of generations as desired. Default values are [population size = 100] and [number of generations = 100].
6. Run the `run_AMALGAM` MATLAB script for the desired case study.
7. The script will execute the optimization algorithm and generate the Pareto set of design alternatives based on the provided inputs.


## Robustness_Analysis

The Robustness_Analysis folder contains the robustness analysis code and setup along with HYPER simulation files  for the five case studies. This code is used to evaluate generated design alternatives robustness under various scenarios. Each case study has its own setup within the folder. Note that scenarios (SOWs) are generated based on the defined ranges in the related paper within each case study setup.

To use the Robustness_Analysis code:

1. Ensure that you have MATLAB installed on your system.
2. Clone or download this repository to your local machine.
3. Open MATLAB and navigate to the Robustness_Analysis folder.
4. The required input files, including historical and synthetic streamflow records as well as project-based parameters, are already provided in the respective case study folders.
4. Run the MATLAB script (main_CASE STUDY NAME)  of the desired case study to perform robustness analysis.
5. The output of the robustness analysis will include a large-sized file containing parameters for both financial and streamflow scenarios and the evaluation of objective functions, energy generation, financial parameters of each alternative, and their robustness results.
Please note that the robustness analysis for each case study was performed on a high-performance computing (HPC) system due to its computational expense. 

## Post-processing the Robustness Results

After running the robustness analysis, you can use the provided MATLAB code file, `Robustness_Metrics.m`,   for post-processing the results. This file is designed to analyze and visualize the robustness results obtained from the Robustness_Analysis code.


## Robustness Results Repository

Since running the case studies are computationally expensive, design alternatives of each case study generated through MO are provided tn this repository as Des_pars_CASE STUDY NAME. Additionally, pre-calculated robustness results for each case study are provided in a separate repository due to their large size. You can access the robustness results repository at [https://drive.google.com/drive/folders/1wA52HRU3jIeiFyXRta67kjQEUmON2D_H?usp=drive_link]. The repository contains the pre-calculated robustness results for each case study, allowing you to review and analyze the robustness without re-running the analysis.

## License

This project is licensed under the MIT License. See the LICENSE file for more details.

## Contact

If you have any questions or need further assistance, please feel free to contact us at [vyildiz1@sheffield.ac.uk].

Happy robustness analysis and post-processing!