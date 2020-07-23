function roiRSA_create_test_MEG_RDMs_groupFFX(option)

warning off

%matlabpool(5);

%% Begin
% Set empty output matrix
load([option.var_dir option.models '.mat']);  % Models
mod_names = fieldnames(Models);
nmods = size(mod_names,1);

for mask = 1:length(option.masknic)
    
    sprintf('......Averaging subjects in Region %s......', num2str(mask))
    
    for sub = 1:length(option.subs)
        
        cd([option.datadir option.sub_beg option.subs{sub} option.subdir]);
        tmpout = zeros(length(option.masknic),nmods,option.epoch_length/option.srate+1);
        
        % Load list of trials to reject
        reject = load([option.datadir option.sub_beg option.subs{sub} option.subdir option.trials option.subs{sub} '.mat']);
        
        %% Load RDM files
        if option.tw > 1
            if mod(option.tw,2)  % must be even number
                option.tw = option.tw+1;
            end
        end
        if option.doavg
            dat = load([option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.srate) 'ms_sTW.mat']);
        else
            dat = load([option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.srate) 'ms_sTW.mat']);
        end
        
        if sub==1; ROI_RDMs = nan([length(option.subs),size(dat.ROI_RDMs)]); end
        
        ROI_RDMs(sub,:,:,:) = dat.ROI_RDMs;
        clear dat
        
    end
    
    % Group average RDM
    ROI_RDMs = squeeze(nanmean(ROI_RDMs,1));
    
    % Save
    outname = [option.rfxdir option.masknic{mask} option.midname 'group_RDMs_' num2str(option.tw*option.srate) 'ms_sTW.mat'];
    save(outname,'ROI_RDMs','option');
    
    
    %% RSA analysis
    % model RDMs
    for m = 1:nmods
        model_tmp = Models.(mod_names{m});
        model(:,m) = vectorizeRDM(model_tmp)'; clear model_tmp
    end
    
    vs = length(ROI_RDMs(1,1,:,1));
    meg_data = zeros((((vs*vs)-vs)/2),length(ROI_RDMs(:,1,1,1)),'single');
    for time = 1:length(ROI_RDMs(:,1,1))
        meg_data_tmp = squeeze(ROI_RDMs(time,1:vs,1:vs));  % get MEG RDM for this timepoint
        meg_data(:,time) = single(vectorizeRDM(meg_data_tmp)'); clear meg_data_tmp
    end
    
    if option.parmodels  % run all models together
        % Remove NaNs
        x = [];
        x = [model, meg_data]; clear meg_data % combine models and data in same matrix
        x(any(isnan(x(:,[1:nmods round(size(x,2)/2)])),2),:) = [];
        
        % correlate
        sprintf('......model correlations......')
        y = x(:,1:nmods); y = double(y);  % seperate for memory saving
        x(:,1:nmods) = []; x = double(x);
        
        if option.partial
            % partial correlation controlling for effect of all models on one model
            tot_mods = 1:nmods;
            for i = tot_mods
                % Get model sets
                exc = tot_mods; inc = tot_mods(i); exc(i) = [];
                pc = partialcorr(x,y(:,inc),y(:,exc),'type', option.dist, 'rows', 'pairwise')';
                r(i,:) = pc;
            end
        else
            r = corr(y,x,'type',option.dist,'rows','pairwise');
        end
        r = single(r);
    else % run each model seperately
        for m = 1:nmods
            % Remove NaNs
            x = [];
            x = [model(:,m), meg_data]; % combine models and data in same matrix
            x(any(isnan(x(:,[1 round(size(x,2)/2)])),2),:) = [];
                        
            % correlate
            tmpr = corr(double(x(:,1)),double(x(:,2:end)),'type',option.dist,'rows','pairwise');
            r(m,:) = single(tmpr); clear tmpr
        end
        clear meg_data
    end
    
    % Save RSA timecourse
    for m = 1:size(mod_names,1);
        rsa_out.(mod_names{m}) = r(m,:);  % store RSA timecourses
    end    
    outname = [option.rfxdir option.rsafront option.masknic{mask} '_groupFFX_' option.midname num2str(option.tw*option.srate) 'ms_sTW.mat'];
    save(outname, 'rsa_out', 'mod_names','option'); clear rsa_out r        
    
end
