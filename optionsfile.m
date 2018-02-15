function option = optionsfile(s,part);

% contains all the options and settings for the RSA ROI analysis
% Changed for SPM8/12 data
%
% Alex Feb2015/Sept2015/March2016/Dec2016

if ~exist('part','var')
    part = 0;
end

if part == 0 || part == 1;
% Directory, subject numbers, subjects to exclude
option.s = s;
option.datadir = '/work/imaging3/MEG/AC_MEG_object_processing/OscRSA_2016/MEG/Subjects/';   % root of where to find subject data
% for RFX dir, see below
option.sub_beg = 'meg11_';
option.subs = {'0233' '0247' '0248' '0251' '0255' '0256' '0258' '0259' '0261' '0262' '0263' '0264' '0265' '0266' '0274'};
option.val = 4;  % SPM inversion to use
option.trialfront = 'dabMefffbasic_'; % file to extract trial data from (can be different to above)
option.subdir1 = '/BD_mf22_corr06_notrans/';  % subfolder where sub data is
option.subdir = [option.subdir1 'ROI_data_Mval' num2str(option.val) '/'];
option.sub_out = [1 3 6];  % index of subjects to exclude from stats from list above


% File parameters reading EEGlab files and output filenames
option.front = 'wmabMefffbasic_';
option.lastbit = '_ssst.mat';
option.trials = ['reject_trials_sub_'];
option.rsafront = 'RS_timecourse_';  % output name after rsa
option.noisefront = 'RS_NC_';
end


if part == 0 || part == 1
% Model RDMs
option.var_dir = '/work/imaging3/MEG/AC_MEG_object_processing/OscRSA_2016/MEG/RSA_ROI/RSA_models/'; % directory where RS models are saved
option.models = ['RSA_sim_modelRDMs_DNN_comb_final'];


% Region definitions
option.usemask = 0;  % use mask (1) or coordinate (0)
option.mask100 = 0;  % use only top 100 vertices in mask (based on zscore of trial data)
option.zscore = 1; % zscore extracted data before saving
option.maskdir = '/work/imaging3/MEG/AC_MEG_object_processing/OscRSA_2016/MEG/ROIs/';
%option.maskname = {'BA17_18.hdr' 'HO_compiled_pVTC_left.hdr' 'HO_compiled_pVTC_right.hdr' 'Holdstock_PRc_2plus_LH.hdr' 'Holdstock_PRc_2plus_RH.hdr'};  % name of mask.nii, should be 1mm spacing
option.maskname = {'HO_compiled_EVC.hdr' 'HO_compiled_pVTC_left.hdr' 'HO_compiled_pVTC_right.hdr' 'HO_compiled_ATL_left.hdr' 'HO_compiled_ATL_right.hdr'};
%option.masknic = {'EVC_' 'LpVTC' 'RpVTC' 'LATL' 'RATL'};
option.masknic = {'OccipCoord_' 'LpVTCCoord_' 'RpVTCCoord_' 'LATLCoord_' 'RATLCoord_'};
option.ROI_coords = {[-10 -94 -16] [-50 -52 -20] [52 -56 -16] [-30 -6 -40] [30 -4 -42]};
option.srad = 20; % radius in mm if using coordinate

% Deal with time confounds
option.remtime = 0; % remove linear and quadratic effects of time from data RDMs before RSA testing
option.timevector = [1:576];  % represents the time-onsets of the trials (can be ordinal, or seconds etc)
option.ndrifts = 2; % number of polynomial drift models to include


% Timing settings
option.srate = 2;        % Sampling rate (in ms; 1000/sample rate in Hz)
option.epoch_length = 2000;    % total length of the epoched data in ms (not including 0 ms)
option.baseline_length = 1500;  % total length of the pre-stimulus period in ms (not including 0 ms)
option.bl_ms = [-200 0];    % time in ms for baseline correction
option.epoch = [1:option.epoch_length/option.srate];
option.baseline = option.baseline_length/option.srate;  % baseline length in samples


% RSA settings
option.tw = round(60/option.srate); % size of sliding time-window (ms/samp_rate), cleanest if tw an even number
option.doavg = 0; % use averaged temporal signal over time window (spatial pattern, 1), or spatial-temporal pattern (0)
option.s = s;
option.rfxdir = ['/work/imaging3/MEG/AC_MEG_object_processing/OscRSA_2016/MEG/RSA_ROI/RSA_timecourses_coord_' option.models '_' num2str(option.s) '/']; % where to save RFX data
option.dist1 = 'correlation';  % distance to use for first level RSA
option.dist = 'Spearman';  % distance to use for second level RSA
option.saveRDMs = 1; % save data RDMs (can be large files)
option.parmodels = 1; % run brain-model correlations in one step, or nmod steps (turn off if models have different pattern of NaNs, but much slower)
option.partial = 1; % run partial correlations (1 - only works in parmodels mode right now)
end


if part == 0 || part == 1
% RFX stats
option.tails = 1;
option.alpha = 0.01/option.tails; % one-tailed significance level for height threshold
option.perms = 10000;
option.stats_epoch = [0 800]; % start and stop of epoch for stats in ms
option.anaepoch = [(option.stats_epoch(1)/option.srate)+option.baseline+1 (option.stats_epoch(2)/option.srate)+option.baseline+1]; % epoch for stats in samples


% TF settings
% (please check extraction script, line 38, to check TF extraction settings)
option.doTF = 1; % Do time-freq version (need both TF and TFbands if TF not been extracted before)
option.doTFbands = 0; % Do time-freq bands version
option.doTWs = 0;     % Use data averaged over timewindows for TFbands analysis
option.doPhase = 1;   % use phase instead of power
option.doPhasePower = 0; % use phase and power
option.timew = {[-500 0] [0 500]};   % TWs to average data over
option.fsc = [4 95]; % frequency limits for TF (also used for bands)
option.nfs = 50; % number of frequencies to extract for log spaced fsc
option.fs = {[4 9] [9 15] [16 30] [31 95]}; % start and stop of bands to do (in Hz)
option.tf_pre = 'TF_';
%option.tf_pre2 = 'rptf_';  % if present it will use both this and above together
option.cycles = 5;
option.tfstep = 20; % time step for TF estimation (ms)
option.nonlocked = 0; % remove evoked signal from the analysis at extraction?

% Plot dimesions
option.start = -option.baseline_length;
option.stop = (option.epoch_length-option.baseline_length);
option.bl_samps = abs(option.bl_ms(1))/option.srate;
option.ylimits = [0 0.04];
option.tlimits = [-6 6];
option.range = [option.baseline-option.bl_samps+1 (option.stats_epoch(2)/option.srate)+option.baseline]; % displays from the start of the period BC was applied, to the end of the period stats were calculated

option.jump = 100;  % points in the x-axis lable in ms


% These shouldn't need changing
option.midname = ['source_'];

end