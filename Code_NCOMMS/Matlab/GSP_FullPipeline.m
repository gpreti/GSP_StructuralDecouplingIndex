
%% Nature Communications submission NCOMMS-19-13468 (METHODS Pipeline)
clear all

%% 0 load variables (structural connectome + functional timecourses), compute & decompose Graph Laplacian, project RS fMRI signals into structural harmonics to find cutoff and split high and low frequencies
findpath=which('GSP_FullPipeline');
findpath=findpath(1:end-18);
mypath=findpath; %set code folder path

GSP_Laplacian 

%% 1 analyze graph signal
GSanalysis

%% 2 generate SC-informed surrogates
GSrandomozation_create_SCinformed_surrogates

%% 3 compute structural decoupling index (SDI) for empirical data and for SC-informed surrogates + test significance of empirical against surrogates
SDI_SCinformed_surrogates

%% 4 generate SC-ignorant surrogates
GSrandomozation_create_SCignorant_surrogates

%% 5 compute SDI for SC-ignorant surrogates
SDI_SCignorant_surrogates

%% 6 Functional Connectivity analysis
GS_FC

%% 7 Create maps for Neurosynth meta-analysis (which will be run in Python)
Neurosynth_inputs
