# GE-PRD
Title: Gene-environment Interaction Identification via Penalized Robust Divergence

Authors: Mingyang Ren, Sanguo Zhang, Shuangge Ma, Qingzhao Zhang

Code maintainer: Mingyang Ren (renmingyang17@mails.ucas.ac.cn)

R version: 4.0.0 (2020-04-24) ( platform: x86_64-w64-mingw32/x64 (64-bit) )
R package versions:
glmnet: 4.0    Matrix: 1.2-18   	grpreg: 3.3.0  	Rcpp: 1.0.4.6

R and cpp files:
The following .R and .cpp files are included in the submitted supplementary files, and the folders in which these files are placed are described in “Structure of code submission”.
(Main functions)
1.function.R: This file includes main functions of proposed methods for joint and marginal analysis, which are used to support numerical simulation studies and real data analysis. Detailed descriptions of these functions can be found in the notes to this R file.
2.Four .cpp files in total named in the “Rcpp_fun_...” or “Rcpp_fun_...” form: These files are base Rcpp functions that implement the proposed methods, which are needed to be loaded before using function.R.

(Codes for real data analysis)
3.TNBC_preprocess.R: This file includes R codes performing the preprocessing of the Triple-Negative Breast Cancer data.
4.TNBC_joint.R: This file includes R codes performing joint analysis of the Triple-Negative Breast Cancer data.
5.TNBC_marginal.R: This file includes R codes performing marginal analysis of the Triple-Negative Breast Cancer data.

(Codes for simulation studies)
6.gen_data.R: This file includes functions for generating simulated data and evaluating performances of competitors, which are used to support numerical simulation studies. Detailed descriptions of these functions can be found in the notes to this R file.
7.main_joint_simulations.R: This file is used to perform joint analysis for logistic regression and linear regression in simulation studies.
8.main_marginal_simulations.R: This file is used to perform marginal analysis for logistic regression and linear regression in simulation studies.

Data:
There is a data file used to perform real data analysis in the submitted supplementary files, and the folders in which these files are placed are described in “Structure of code submission”.
TNBC.name.list.RData: This file contains the name lists of genes and samples in the Triple Negative Breast Cancer working data, which can be obtained by running codes available at “http://web.tecnico.ulisboa.pt/susanavinga/TNBC/Lopes_et_al_BMCBioinformatics.Rmd” (Lopes et al., 2018).

How to conduct real data analysis:
The files “TNBC_preprocess.R”, “TNBC_joint.R”, “TNBC_marginal.R” can be used to conduct real data analysis. 

For Triple-Negative Breast Cancer data analysis, below we describe the procedure.
1.Put the data file “TNBC.name.list.RData” in the working directory, such as “case_study”.
2.Run all code in “TNBC_preprocess.R”. Then we can obtain the working data saved as "TNBC.RData".
(joint analysis)
3.Using “TNBC_joint.R”, load data “TNBC.RData”, and source the function.R and the relevant .cpp files (Rcpp_fun_mdpdlogistic.cpp, Rcpp_fun_gammalinear.cpp).
4.Run all codes in “TNBC_joint.R”. Then we can obtain the results which can be saved in .csv files.
(marginal analysis)
5.Using “TNBC_marginal.R”, load data “TNBC.RData”, and source the function.R and the relevant .cpp files (Rcpp_fun_mdpdlogistic.cpp, Rcpp_fun_gammalinear.cpp).
6.Run all codes in “TNBC_marginal.R”. Then we can obtain the results which can be saved in .csv files.


How to conduct simulations:
The files main_joint_simulations.R and main_marginal_simulations.R are used to conduct simulation studies.

Take main_joint_simulations.R as an example. Below we describe how to conduct joint analysis. The marginal analysis can be conducted in the same manner.
1.Load R libraries glmnet, grpreg, and Rcpp. 
2.Source the following files: 
function.R,		              gen_data.R,	
Rcpp_fun_mdpdlogistic.cpp,		Rcpp_fun_gammalogistic.cpp,		
Rcpp_fun_gammalinear.cpp,   	Rcpp_fun_mdpdlinear.cpp.
Run all the code in main_joint_simulations.R. Then we can obtain the simulation results.

References
1.Lopes, M. B., Veríssimo, A., Carrasquinha, E., Casimiro, S., Beerenwinkel, N., & Vinga, S. (2018). Ensemble outlier detection and gene selection in triple-negative breast cancer data. BMC bioinformatics, 19(1), 168.
