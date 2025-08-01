# RelativePhaseAnalysis
This repository provides core code for Relative Phase Analysis (RPA).

The materials here are produced by the "Brain States and Transitions Laboratory": https://moonbrainlab.org
All questions and inquiries are to be directed to the Principal Investigator, Dr. Joon-Young Moon.
Please refer to the above link for the contact address.

The codes are written in MATLAB, and this repository includes the necessary .m files and .mat files to perform RPA.
The RPA will be applied to a preprocessed, band-pass filtered multi-channel EEG/ECoG time series.
We will assume that the EEG/ECoG time series under investigation is already preprocessed.

The main script to run is "run_RPA.m"
The names of the directories and files in "run_RPA.m" must be modified accordingly to the user's environment.

The order of analysis included in the "run_RPA.m" is as follows for the pipeline #1, 
where the centroids are constructed from the given data itself: 

1. make_movie_rel_phase.m
This script performs RPA and produce Relative Phase Dynamics Movie.

2. cluster_movie_frames_eval.m
This script performs K-means clustering across the frames of the Relative Phase Dynamics Movie.

3. cluster_movie_frames_eval_PCA.m
This script performs PCA across the frames of the Relative Phase Dynamics Movie. 


The order of analysis included in the "run_RPA.m" is as follows for the pipeline #2, 
where the centroids are pre-constructed: 

1. make_movie_rel_phase.m
This script performs RPA and produce Relative Phase Dynamics Movie.

2. cluster_movie_frames_eval_reg.m
This script performs linear regression, using the pre-constructed "universal centroids" as the regressors, across the frames of the Relative Phase Dynamics Movie.
It is equivalent to performing PCA across the frames of the Relative Phase Dynamics Movie, using the "universal centroids" as the eigenvectors. 



Other codes necessary to perform the analysis are the following:

1. cal_rel_phase_v_final.m
This function calculates relative phase from the given time series.

2. change_cluster_idx_v_final.m
This function changes the order of clusters automatically.

3. moving_time_window.m
This function performs moving time window averaging.

4. cal_regression_clustering.m
This function calculates the cluster indices of the given time series.

6. cal_transition_prop_v_final.m
This function calculates the transition matrix and dwelling time of the given time series.


Directories necessary to perform the analysis are the following:

1. /topoplot_script
This directory contains all the data to construct Relative Phase Dynamics Moive.
This directory also contains function "topoplot_general_test.m", which is to be utilized in other functions and scripts.

2. /comb_centroids_v1
This directory contains pre-constructed centroids, which can be utilized in the RPA. 
