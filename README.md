# CARE

Classical MR analyses utilizing summary data may still produce biased causal effect estimates due to the the winner’s curse and pleiotropic issues. To address these two issues and establish valid causal conclusions, we propose a unified robust Mendelian Randomization framework with summary data, which systematically removes the winner’s curse and screens out invalid genetic instruments
with pleiotropic effects. 

In this repo, we provide the following sources.

- Simulations folder: the codes for replicating simulation results in the above manuscript
- RealData folder: the codes for replicating real data results in the above manuscript

We aim to write a separate manuscript that focuses on software and pipeline development. The codes provided here (especially under the RealData folder) can provide a good template for experts using our new methods. The users may need to specify the working directory and install relevant packages to run it smoothly.



## Simulation
- Simulation/simulation_dir_pleio2.R: Generate simulation results for the main setting.
- Simulation/simulation_dir_pleio.R: Generate simulation results for uniform distributed effects for alpha.
- Simulation/simulation_balance_pleio.R: Generate simulation results for balanced horizontal pleiotropy with InSIDE assumption satisfied.
- Simulation/simulation_balance_pleio2.R: Generate simulation results for balanced horizontal pleiotropy wih InSIDE assumption violated.
- Simulation/simulation_nonlinear.R: Generate simulation results for nonlinear X on G without interaction term. 
- Simuation/simulation_nonlinear2.R: Generate simulation results for nonlinear X on G with interaction term. 
- Simulation/simulation_nonlinear3.R: Generate simulation results for nonlinear Y on X without interaction term.
- Simulation/simulation_dir_pleio3.R/simulation_dir_pleio7.R: Generate simulation results for GBIC.
- Simulation/simulation_dir_pleio4.R: Different sample sizes.
- Simulation/simulation_dir_pleio5.R: Third sample under the main setting. P cutoff = 5e-5 for all methods
- Simulation/simulation_dir_pleio7.R: Different eta under the main setting.
- Simulation/simulalation_dir_pleio8.R: Generate simulation results for different sample sizes of SNPs.
- simulation_addition/simulation_dir_pleio9.R: Comparison of l0, and two l1 algorithms under the main setting.
- Simulation/simulation_dir_pleio10.R: P cutoff = 5e-5 for other methods.
- Support files: CARE_support.R, CARE_support2.R, CARE_support3.R, mr_raps_own.R, mr_lasso_own.R, cML_support.R, CARE_support_measurement_with_CD2.cpp, CARE_support_measurement_with_CD2_2.cpp, CARE_support_measurement_overlap.cpp, cML_support.cpp
- /pbs: generate lsf files used for parallel computing.

### Figure 2 (Main setting)
- resAna/summarize_simRes_final.R: Generate summary results.
- resAna/simulation_figure2.R: Reproduce results in figure 2.
- summaryRes/: This folder contains simulation results under the main setting. Due to file upload size limitations, supplementary results are not included here.

### Supplementary figures

- Simulation/resAna/summarize_simRes_final.R: Generate summary results.
- Simulation/resAna/summarize_simRes_final2.R: Generate summary results for GBIC.
- Simulation/resAna/summarize_simRes_final_eta.R: Generate summary results for different eta.
- simulation_addition/resAna/summarize_simRes_final.R: Generate summary results for comparison of l0 and two l1 algorithms under the main setting.
- Simulation/resAna/simulation_supplementary_figure.R: Reproduce results in Supplementary Figure 1.
- Simulation/resAna/simulation_figure2.R: Reproduce results in Supplementary Figure 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16, 17, 18, 19, 28.
- Simulation/resAna/simulation_figure3.R: Reproduce results in Supplementary Figure 20, 21, 22, 23, 24, 25, 26, 27.
- Simulation/resAna/simulation_figure4.R: Reproduce results in Supplementary Figure 29, 30.
- Simulation/resAna/simulation_supplementray_eta_figure.R: Reproduce results in Supplementary Figure 14.
- Simulation/resAna/simulation_supplementray_gic_figure.R: Reproduce results in Supplementary Figure 15.
- simulation_addition/resAna/simulation_figure.R: Reproduce results in Supplementary Figure 29, 30
- Simulation/resAna/summarize_simRes_computational_time_setting: Generate summary results for running time.
- Simulation/resAna/simulation_supplementary_figure_runtime.R: Reproduce results in Supplementary Figure 4.



## Real data 
### Negative Control

#### Prepare data

- Download data from IEU-GWAS: https://gwas.mrcieu.ac.uk/. The exposures we used are listed in exposure_list.csv.
- step0b_convert_data.R: Convert raw summary data to txt files.
- step1_process_GWAS.R: Pre-processing and quality control codes for GWAS summary data. This is a very useful in-house code to pre-process GWAS data. If you use this for other purposes, please cite this paper.

#### Analysis 

- MR_pipeline.R: Main codes for applying CARE and other MR methods for analysis.

#### Generate results

- step5_resAna_all.R: Reproduce results in Figure 3 and Supplementary Figure 14.
- smmry_res/: Due to file upload size limitations, this folder only contains 20 results generated from step5_resAna_all.R.

#### Support files

- CARE_support.R, cML_support.R, mr_raps_own.R, CARE_support.cpp, CARE_support_measurement_with_CD2.cpp, Sfig1_qq_plot_inflated_type1.R.
- exposure_list.csv: All the exposures used in this study.

### COVID-19 Severity

#### Download outcome data
- http://covid19hg.org/results/r7/

#### Prepare data
- step0a_download_IEUGWAS.R: Download exposures from IEU-GWAS.
- step0b_convert_data.R: Convert raw data to txt files.
- step1_process_GWAS.R: Pre-processing and quality control codes for GWAS summary data. This is a very useful in-house code to pre-process GWAS data. If you use this for other purposes, please cite this paper.

#### Analysis
- resAna/resAna_step1_all_summarized_COVID19.R: Generate summary results for analysis use.
- COVID19_summarized.RData contains the results from resAna_step1_all_summarized_COVID19.R.

#### Generate results
- resAna/COVID_histgram.R: Reproduce results in Figure 4.
- resAna/COVID_DataPlots.R: Reproduce results in Figure 5.

## Disclaimer

The codes are provided "as is" and the author disclaims all warranties with regard to these codes including all implied warranties of merchantability and fitness. In no event shall the author be liable for any special, direct, indirect, or consequential damages or any damages whatsoever resulting from loss of use, data or profits, whether in an action of contract, negligence or other tortious action, arising out of or in connection with the use or performance of these codes. 

