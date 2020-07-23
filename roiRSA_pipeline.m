%% Run RSA pipeline
%
% Alex Feb2015/Sept2015/March2016/Dec2016 etc

% Set stimulus selections to use
set = {[1:51]};

rs=rng; %set random seed for stats replication

%% Task switches
ecog_extract = 1;    % Extract from 3D ecog matrix
sens_extract = 0;    % Extract sensor timecourses
extract      = 0;    % Extract source ROI timecourses for vertices
RSA_do       = 1;    % Calculate data RDMs and correlate with models
RSA_do_gffx  = 0;    % RSA based on RDMs averaged over subjects
ROI_noise    = 0;    % Calculate noice ceilings for each region
rfx_stats    = 1;    % Do permutation stats
rfx_plots    = 1;    % Generate plots


%% Add path to RSA tools and get settings
addpath(genpath('/work/imaging3/MEG/AC_ECOG_objects/scripts/RSA_MEEG_tools-master'));
addpath('/home/alex/software/EEGLAB/eeglab13_4_4b/');
%eeglab;


%% Run tasks
for s = 1:length(set)  
    option = optionsfile(s);
    
    if ecog_extract; roiRSA_ecog_extract(option,s); end
    if sens_extract; roiRSA_sensor_extract(option,s); end
    if extract; roiRSA_source_extract(option,s); end            

    if RSA_do; roiRSA_create_test_MEG_RDMs(option,set{s},s); end;     
    if RSA_do_gffx; roiRSA_create_test_MEG_RDMs_groupFFX(option); end    
    if ROI_noise; roiRSA_noise_ceilings(option,s); end
    
    if rfx_stats; roiRSA_permstats(option,rs); end   % for parametric rfx
%   if rfx_stats; roiRSA_permstats_wilcoxon(option,rs); end
    if rfx_plots; roiRSA_rfxplots(option,s,1); end    

end

%clear