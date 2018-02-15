%% Run RSA ROI pipeline
%
% Alex Feb 2015

% Set stimulus selections to use
load('/work/imaging3/MEG/AC_MEG_object_processing/PCA_2011/Youngs/MEG/regressors/trialorder_to_catcode_order.txt');
set = {trialorder_to_catcode_order'};
%set = {[1:302]};

%% Task switches
sens_extract = 0;    % Extract sensor timecourses
extract      = 0;    % Extract source ROI timecourses for vertices
RSA_do       = 0;    % Calculate data RDMs and correlate with models
ROI_noise    = 0;    % Calculate noice ceilings for each region
rfx_stats    = 0;    % Do permutation stats
rfx_plots    = 1;    % Generate plots


% other tasks
createRDMs_only = 0;  % Create MEG ROI RDMs
rsa_calc_only = 0;    % Create RSA timecourses for models


%% Add path to RSA tools and get settings
addpath('/work/imaging3/MEG/AC_MEG_object_processing/OscRSA_2016/scripts/RS_spm12');
addpath('/home/alex/software/EEGLAB/eeglab13_4_4b/');
%eeglab;


%% Run tasks
for s = 1:length(set)  
    option = optionsfile(s);
    
    if sens_extract; roiRSA_sensor_extract(option,s); end
    if extract; roiRSA_source_extract(option,s); end            
    if RSA_do; roiRSA_create_test_MEG_RDMs(option,set{s},s); end;    
        if createRDMs_only; roiRSA_create_MEG_RDMs(option,set{s},s); end;
        if rsa_calc_only; roiRSA_timecourses(option,s); end;
    if ROI_noise; roiRSA_noise_ceilings(option,s); end
    if rfx_stats; roiRSA_permstats(option); end
    if rfx_plots; roiRSA_rfxplots(option,s,1); end    
%   if rfx_plots; roiRSA_rfxplots_clusterTF_best(option,s,1,c); end    
end

clear all

% Could also add:
%   Data RDM temporal correlations
%   Noise ceilings
%   TFCE