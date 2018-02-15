function roiRSA_create_test_MEG_RDMs(option,trials,s)

% Creates MEG-based ROI RDMs from MEG source localised ROI data
%
% Alex (10/2014) Feb 2015

%warning off

matlabpool(5);

%% Begin
% Set empty output matrix
optionz = option; % for parfor
load([option.var_dir option.models '.mat']);  % Models
mod_names = fieldnames(Models);
nmods = size(mod_names,1);

parfor sub = 1:length(optionz.subs)
    option = optionz;
    
    if option.doTFbands
        %   start = (1000/min(option.fsc))*(round(option.cycles/2));
        tmpout = zeros(length(option.masknic),nmods,size(option.fs,2),(option.epoch_length)/option.tfstep+1);
    elseif option.doTF
        %   start = (1000/min(option.fsc))*(round(option.cycles/2));
        tmpout = zeros(length(option.masknic),nmods,option.nfs,(option.epoch_length)/option.tfstep+1);
    else
        tmpout = zeros(length(option.masknic),nmods,option.epoch_length/option.srate+1);
    end
    Models = [];
    
    % Load list of trials to reject
    reject = load([option.datadir option.sub_beg option.subs{sub} option.subdir option.trials option.subs{sub} '.mat']);
    
    cd([option.datadir option.sub_beg option.subs{sub} option.subdir]);
    
    for mask = 1:length(option.masknic)
        
        sprintf('......Subject %s, Region %s......', option.subs{sub},num2str(mask))
        
        % Need for files extracted with tw as one value, and now want to
        % use another without re-extracting
        options = optionsfile(s,1);
        option = setfield(option,'masknic',options.masknic);
        option = setfield(option,'ROI_coords',options.ROI_coords);
        option = setfield(option,'tw',options.tw);
        option = setfield(option,'doTFbands',options.doTFbands);
        option = setfield(option,'doTF',options.doTF);
        option = setfield(option,'fs',options.fs);
        option = setfield(option,'tw',options.tw);
        option = setfield(option,'doavg',options.doavg); options = [];
        
        %% Setup matfile for outputs
        if option.doTF && option.doTFbands==0
            option.tw = round((option.tw*option.srate)/option.tfstep);
            if option.tw > 1
                if mod(option.tw,2)  % must be even number
                    option.tw = option.tw+1;
                end
            end
            if option.doavg
                if option.doPhase || option.doPhasePower
                    if exist([option.tf_pre 'phz_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'file'); delete([option.tf_pre 'phz_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']); end
                    rdms = matfile([option.tf_pre 'phz_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'Writable',true);
                else
                    if exist([option.tf_pre option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'file'); delete([option.tf_pre option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']); end
                    rdms = matfile([option.tf_pre option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'Writable',true);
                end
            else
                if option.doPhase || option.doPhasePower
                    if exist([option.tf_pre 'phz_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'file'); delete([option.tf_pre 'phz_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']); end
                    rdms = matfile([option.tf_pre 'phz_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'Writable',true);
                else
                    if exist([option.tf_pre option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'file'); delete([option.tf_pre option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']); end
                    rdms = matfile([option.tf_pre option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'Writable',true);
                end
            end
        elseif option.doTFbands
            if option.doTWs
                option.tw = 1;
            else
                option.tw = round((option.tw*option.srate)/option.tfstep);
            end
            if option.tw > 1
                if mod(option.tw,2)  % must be even number
                    option.tw = option.tw+1;
                end
            end
            if option.doavg
                if option.doPhase || option.doPhasePower
                    if exist([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'file'); delete([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']); end
                    rdms = matfile([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'Writable',true);
                else
                    if exist([option.tf_pre 'bands_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'file'); delete([option.tf_pre 'bands_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']); end
                    rdms = matfile([option.tf_pre 'bands_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'Writable',true);
                end
            else
                if option.doPhase || option.doPhasePower
                    if option.doTWs
                        if exist([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'TW_RDMs.mat'],'file'); delete([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'TW_RDMs.mat']); end
                        rdms = matfile([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'TW_RDMs.mat'],'Writable',true);
                    else
                        if exist([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'file'); delete([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']); end
                        rdms = matfile([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'Writable',true);
                    end
                else
                    if option.doTWs
                        if exist([option.tf_pre 'bands_' option.masknic{mask} option.midname 'TW_RDMs.mat'],'file'); delete([option.tf_pre 'bands_' option.masknic{mask} option.midname 'TW_RDMs.mat']); end
                        rdms = matfile([option.tf_pre 'bands_' option.masknic{mask} option.midname 'TW_RDMs.mat'],'Writable',true);
                    else
                        if exist([option.tf_pre 'bands_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'file'); delete([option.tf_pre 'bands_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']); end
                        rdms = matfile([option.tf_pre 'bands_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'Writable',true);
                    end
                end
            end
        else
            if option.tw > 1
                if mod(option.tw,2)  % must be even number
                    option.tw = option.tw+1;
                end
            end
            if option.doavg
                if exist([option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.srate) 'ms_sTW.mat'],'file'); delete([option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.srate) 'ms_sTW.mat']); end
                rdms = matfile([option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.srate) 'ms_sTW.mat'],'Writable',true);
            else
                if exist([option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.srate) 'ms_sTW.mat'],'file'); delete([option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.srate) 'ms_sTW.mat']); end
                rdms = matfile([option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.srate) 'ms_sTW.mat'],'Writable',true);
            end
        end
        
        %% Read ROI data
        if option.doTF || option.doTFbands
            if option.doPhase || option.doPhasePower
                infile = [option.tf_pre 'phz_' option.masknic{mask} option.midname option.front option.subs{sub} '.mat'];
                D = load(infile);
            else
                infile = [option.tf_pre option.masknic{mask} option.midname option.front option.subs{sub} '.mat'];
                D = load(infile);
            end
        else
            infile = [option.masknic{mask} option.midname option.front option.subs{sub} '.mat'];
            D = load(infile);
        end
        
        % Need for files extracted with one value, and now want to
        % use another without re-extracting
        options = optionsfile(s,0);
        option = setfield(option,'masknic',options.masknic);
        option = setfield(option,'ROI_coords',options.ROI_coords);
        option = setfield(option,'subs',options.subs);
        option = setfield(option,'partial',options.partial);
        option = setfield(option,'tw',options.tw);
        option = setfield(option,'models',options.models);
        option = setfield(option,'doTFbands',options.doTFbands);
        option = setfield(option,'doTF',options.doTF);
        option = setfield(option,'doavg',options.doavg);
        option = setfield(option,'fs',options.fs);
        option = setfield(option,'saveRDMs',options.fs);
        options = [];
        
        % Create RDMs over time
        %% TF version
        if option.doTF || option.doTFbands
            if option.doTF
                if option.doTWs
                    option.tw = 1;  % tw needs to be 1 if we'll do TW averaging later
                else
                    option.tw = round((option.tw*option.srate)/option.tfstep);
                    if option.tw > 1
                        if mod(option.tw,2)  % must be even number
                            option.tw = option.tw+1;
                        end
                    end
                end
                
                % Read all data first
                ao = D.D(:,:,:,trials);
                at = length(D.D(1,1,:,1)); % times
                af = length(D.D(1,:,1,1)); % freqs
                D = [];
                
                %                 % Vectorise phase/power data
                %                 if ~isreal(ao)
                %                     ao = [zscore(angle(ao),0,3); zscore(abs(ao),0,3)];
                %                     %ao = [angle(ao); abs(ao)]; % the spatial dimension is effectively doubled
                %                 end
                
                if option.tw > 1
                    for time = (option.tw/2+1):at-(option.tw/2)
                        if option.doavg
                            a = mean(ao(:,:,(time-option.tw/2):(time+option.tw/2),:),3);  % vertices x freq x time x trials
                        else
                            a = ao(:,:,(time-option.tw/2):(time+option.tw/2),:);  % vertices x freq x time x trials
                        end
                        a = reshape(permute(a,[2 4 1 3]),af,length(trials),[]); % freqs x trials x [time x vertices]
                        
                        % RDMs for each freq
                        for f = 1:af
                            
                            if ~isreal(a) % i.e. power/phase distances
                                st_rdm = complexdist(squeeze(a(f,:,:)));
                            else  % normal distances
                                st_rdm = squareform(pdist(squeeze(a(f,:,:)),option.dist1));
                            end
                            
                            % set incorrect/rejected as NaNs
                            st_rdm(find(reject.reject(trials)),:) = NaN;
                            st_rdm(:,find(reject.reject(trials))) = NaN;
                            
                            % collect data
                            temp_o(1,1,:,:) = single(st_rdm);
                            rdms.ROI_RDMs(time,f,1:size(st_rdm,1),1:size(st_rdm,2)) = temp_o;
                            st_rdm = []; temp_o = [];
                        end
                        a = [];
                    end
                    rdms.ROI_RDMs(at-(option.tw/2):at,1:af,1:length(rdms.ROI_RDMs(1,1,:,1)),1:length(rdms.ROI_RDMs(1,1,:,1))) = NaN;
                else
                    for time = 1:at
                        a = ao(:,:,time,:);  % vertices x freqs x time x trials
                        a = reshape(permute(a,[2 4 1 3]),af,length(trials),[]); % freqs, trials x [time x vertices]
                        
                        % RDMs for each freq
                        for f = 1:af
                            
                            if ~isreal(ao) % i.e. power/phase distances
                                st_rdm = complexdist(squeeze(a(f,:,:)));
                            else  % normal distances
                                st_rdm = squareform(pdist(squeeze(a(f,:,:)),option.dist1));
                            end
                            
                            % set incorrect/rejected as NaNs
                            st_rdm(find(reject.reject(trials)),:) = NaN;
                            st_rdm(:,find(reject.reject(trials))) = NaN;
                            
                            % collect data
                            temp_o(1,1,:,:) = single(st_rdm);
                            rdms.ROI_RDMs(time,f,1:size(st_rdm,1),1:size(st_rdm,2)) = temp_o;
                            st_rdm = []; temp_o = [];
                        end % f
                        a = [];
                    end % time
                end
                ao = [];
            else
                % Load TF data from previous
                if option.doTWs
                    option.tw = 1;  % tw needs to be 1 if we'll do TW averaging later
                else
                    option.tw = round((option.tw*option.srate)/option.tfstep);
                    if option.tw > 1
                        if mod(option.tw,2)  % must be even number
                            option.tw = option.tw+1;
                        end
                    end
                end
                if option.doavg
                    if option.doPhase
                        rdms = matfile([option.tf_pre 'phz_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                    else
                        rdms = matfile([option.tf_pre option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                    end
                else
                    if option.doPhase
                        rdms = matfile([option.tf_pre 'phz_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                    else
                        rdms = matfile([option.tf_pre option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                    end
                end
            end
            
            %% Time/frequency averaing of similarities
            if option.doTFbands
                
                % Create empty matrices
                if option.doTWs
                    ROI_RDMs_TW = zeros(length(option.timew),length(option.fs),length(rdms.ROI_RDMs(1,1,:,1)),length(rdms.ROI_RDMs(1,1,:,1)));
                else
                    ROI_RDMs_TW = zeros(length(rdms.ROI_RDMs(:,1,1,1)),length(option.fs),length(rdms.ROI_RDMs(1,1,:,1)),length(rdms.ROI_RDMs(1,1,:,1)));
                end
                
                % Fband averaging
                for band = 1:length(option.fs)
                    
                    f = option.fs{band};
                    f = find(((option.fsc>=f(1))&(option.fsc<=f(end)))); % index of freqs in band
                    
                    if option.doTWs % Time window averaging
                        for tw = 1:length(option.timew)
                            t = find(((option.tfepoch>=option.timew{tw}(1))&(option.tfepoch<=option.timew{tw}(end))));
                            ROI_RDMs_TW(tw,band,1:length(rdms.ROI_RDMs(1,1,:,1)),1:length(rdms.ROI_RDMs(1,1,1,:))) = mean(mean(rdms.ROI_RDMs(t,f,1:length(rdms.ROI_RDMs(1,1,:,1)),1:length(rdms.ROI_RDMs(1,1,1,:))),1),2);
                        end
                    else
                        ROI_RDMs_TW(1:length(rdms.ROI_RDMs(:,1,1,1)),band,1:length(rdms.ROI_RDMs(1,1,:,1)),1:length(rdms.ROI_RDMs(1,1,:,1))) = mean(rdms.ROI_RDMs(1:length(rdms.ROI_RDMs(:,1,1,1)),f,1:length(rdms.ROI_RDMs(1,1,:,1)),1:length(rdms.ROI_RDMs(1,1,1,:))),2);
                    end
                end
                
                % Delete old matfile and make new TFbands one (ineligant way to get around replacing data in matfile)
                if option.doavg
                    if option.doPhase || option.doPhasePower
                        delete([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                        rdms = matfile([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'Writable',true);
                    else
                        delete([option.tf_pre 'bands_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                        rdms = matfile([option.tf_pre 'bands_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'Writable',true);
                    end
                else
                    if option.doPhase || option.doPhasePower
                        if option.doTWs
                            delete([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'TW_RDMs.mat']);
                            rdms = matfile([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'TW_RDMs.mat'],'Writable',true);
                        else
                            delete([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                            rdms = matfile([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'Writable',true);
                        end
                    else
                        if option.doTWs
                            delete([option.tf_pre 'bands_' option.masknic{mask} option.midname 'TW_RDMs.mat']);
                            rdms = matfile([option.tf_pre 'bands_' option.masknic{mask} option.midname 'TW_RDMs.mat'],'Writable',true);
                        else
                            delete([option.tf_pre 'bands_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                            rdms = matfile([option.tf_pre 'bands_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'Writable',true);
                        end
                    end
                end
                
                rdms.ROI_RDMs(1:length(ROI_RDMs_TW(:,1,1,1)),1:band,1:length(ROI_RDMs_TW(1,1,:,1)),1:length(ROI_RDMs_TW(1,1,1,:))) = ROI_RDMs_TW(1:length(ROI_RDMs_TW(:,1,1,1)),1:band,1:length(ROI_RDMs_TW(1,1,:,1)),1:length(ROI_RDMs_TW(1,1,1,:)));
                clear ROI_RDMs_TW
            end
            
            %% Do normal
        else
            
            if option.tw > 1
                if mod(option.tw,2)  % must be even number
                    option.tw = option.tw+1;
                end
            end
            
            % Read all data first
            ao = D.D(:,:,trials);
            at = length(D.D(1,:,1)); % times
            D = [];
            
            if option.tw > 1
                for time = (option.tw/2+1):at-(option.tw/2)
                    if option.doavg
                        a = mean(ao(:,(time-option.tw/2):(time+option.tw/2),:),2);  % vertices x time x trials
                    else
                        a = ao(:,(time-option.tw/2):(time+option.tw/2),:);  % vertices x time x trials
                    end
                    a = reshape(permute(a,[3 1 2]),length(trials),[]); % trials x [time x vertices]
                    st_rdm = squareform(pdist(a,option.dist1));
                    
                    % set incorrect/rejected as NaNs
                    st_rdm(find(reject.reject(trials)),:) = NaN;
                    st_rdm(:,find(reject.reject(trials))) = NaN;
                    
                    % collect data
                    temp_o(1,:,:) = single(st_rdm);
                    rdms.ROI_RDMs(time,1:size(st_rdm,1),1:size(st_rdm,2)) = temp_o;
                    a=[]; st_rdm=[]; temp_o=[];
                end
                rdms.ROI_RDMs(at-(option.tw/2):at,1:length(rdms.ROI_RDMs(1,1,:,1)),1:length(rdms.ROI_RDMs(1,1,:,1))) = NaN;
            else
                for time = 1:at
                    a = ao(:,time,:);  % vertices x time x trials
                    a = reshape(permute(a,[3 1 2]),length(trials),[]); % trials x [time x vertices]
                    st_rdm = squareform(pdist(a,option.dist1));
                    
                    % set incorrect/rejected as NaNs
                    st_rdm(find(reject.reject(trials)),:) = NaN;
                    st_rdm(:,find(reject.reject(trials))) = NaN;
                    
                    % collect data
                    temp_o(1,:,:) = single(st_rdm);
                    rdms.ROI_RDMs(time,1:size(st_rdm,1),1:size(st_rdm,2)) = temp_o;
                    a=[]; st_rdm=[]; temp_o=[];
                end
            end
        end % doTF
        
        ao=[];
        
        options = optionsfile(s,1);
        option = setfield(option,'saveRDMs',options.saveRDMs); options =[];
        
        %% Run RSA
        [r mod_names] = RSA_timecourses(option,s,sub,rdms);
        tmpout(mask,:,:,:) = r;        
        rdms = [];
        
        if option.saveRDMs
            rdms.option = option;
            %% Add info to saved matfiles
            if option.doTF && option.doTFbands==0
                rdms.freqs = option.fsc;
            elseif option.doTFbands
                rdms.freqs = option.fs;
            end
        elseif option.saveRDMs == 0
            % delete matfile
            if option.doTF && option.doTFbands==0
                if option.doavg
                    if option.doPhase || option.doPhasePower
                        delete([option.tf_pre 'phz_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                    else
                        delete([option.tf_pre option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                    end
                else
                    if option.doPhase || option.doPhasePower
                        delete([option.tf_pre 'phz_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                    else
                        delete([option.tf_pre option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                    end
                end
            elseif option.doTFbands
                if option.doavg
                    if option.doPhase || option.doPhasePower
                        delete([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                    else
                        delete([option.tf_pre 'bands_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                    end
                else
                    if option.doPhase || option.doPhasePower
                        if option.doTWs
                            delete([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'TW_RDMs.mat']);
                        else
                            delete([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                        end
                    else
                        if option.doTWs
                            delete([option.tf_pre 'bands_' option.masknic{mask} option.midname 'TW_RDMs.mat']);
                        else
                            delete([option.tf_pre 'bands_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                        end
                    end
                end
            else
                if option.doavg
                    delete([option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.srate) 'ms_sTW.mat']);
                else
                    delete([option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.srate) 'ms_sTW.mat']);
                end
            end
        end
        
    end % ROI
    parsave('temp.mat',tmpout); tmpout = [];
end % sub

mod_names = fieldnames(Models);
option = optionz;

% Compile subject data
for sub = 1:length(optionz.subs)
    load([option.datadir option.sub_beg option.subs{sub} option.subdir 'temp.mat'])
    subtmpout(sub,:,:,:,:) = tmpout;
end
tmpout = subtmpout; clear subtempout
tmpout = permute(tmpout,[2 1 3 4 5]);

% Save RSA timecourses
for mask = 1:length(option.masknic)
    
    for m = 1:size(mod_names,1);
        rsa_out.(mod_names{m}) = squeeze(tmpout(mask,:,m,:,:));  % store RSA timecourses
    end
    
    % Added incase you've changed output directory since first data extraction
    options = optionsfile(s,1);
    option = setfield(option,'rfxdir',options.rfxdir);
    clear options
    rfxdir = option.rfxdir;
    if ~exist(rfxdir)
        mkdir(rfxdir)
    end
    
    % Save output
    if option.doTF && option.doTFbands==0
        if option.doPhase || option.doPhasePower
            outname = [option.rfxdir option.tf_pre 'phz_' option.rsafront option.masknic{mask} option.midname num2str(option.tw*option.tfstep) 'ms_sTW.mat'];
            save(outname, 'rsa_out', 'mod_names','option'); clear rsa_out
        else
            outname = [option.rfxdir option.tf_pre option.rsafront option.masknic{mask} option.midname num2str(option.tw*option.tfstep) 'ms_sTW.mat'];
            save(outname, 'rsa_out', 'mod_names','option'); clear rsa_out
        end
    elseif option.doTFbands
        if option.doTWs
            if option.doPhase || option.doPhasePower
                outname = [option.rfxdir option.tf_pre 'bands_phz_' option.rsafront option.masknic{mask} option.midname '.mat'];
                save(outname, 'rsa_out', 'mod_names','option'); clear rsa_out
            else
                outname = [option.rfxdir option.tf_pre 'bands_' option.rsafront option.masknic{mask} option.midname '.mat'];
                save(outname, 'rsa_out', 'mod_names','option'); clear rsa_out
            end
        else
            if option.doPhase || option.doPhasePower
                outname = [option.rfxdir option.tf_pre 'bands_phz_' option.rsafront option.masknic{mask} option.midname num2str(option.tw*option.tfstep) 'ms_sTW.mat'];
                save(outname, 'rsa_out', 'mod_names','option'); clear rsa_out
            else
                outname = [option.rfxdir option.tf_pre 'bands_' option.rsafront option.masknic{mask} option.midname num2str(option.tw*option.tfstep) 'ms_sTW.mat'];
                save(outname, 'rsa_out', 'mod_names','option'); clear rsa_out
            end
        end
    else
        outname = [option.rfxdir option.rsafront option.masknic{mask} option.midname num2str(option.tw*option.srate) 'ms_sTW.mat'];
        save(outname, 'rsa_out', 'mod_names','option'); clear rsa_out
    end
end % mask


function [r mod_names] = RSA_timecourses(option,s,sub,rdms)

% Correlates data RDMs with model RDMs and returns RSA timecourses
% AC (11/2012)(10/2014) Feb 2016

% Load subject-specific model RDMs
load([option.var_dir option.models '.mat']);  % Models
mod_names = fieldnames(Models);
nmods = size(mod_names,1);
for m = 1:nmods
    model_tmp = Models.(mod_names{m});
    model(:,m) = vectorizeRDM(model_tmp)'; clear model_tmp
end

sprintf('......Subject %s, RSA timecourses......', option.subs{sub})

if or(option.doTF,option.doTFbands)
    if option.parmodels  % run all models together
        vs = length(rdms.ROI_RDMs(1,1,:,1));
        meg_data = zeros((((vs*vs)-vs)/2),length(rdms.ROI_RDMs(1,:,1,1)),length(rdms.ROI_RDMs(:,1,1,1)),'single');
        for f = 1:length(rdms.ROI_RDMs(1,:,1,1))
            for time = 1:length(rdms.ROI_RDMs(:,1,1,1))
                meg_data_tmp = squeeze(rdms.ROI_RDMs(time,f,1:vs,1:vs));  % get MEG RDM for this timepoint
                meg_data(:,f,time) = single(vectorizeRDM(meg_data_tmp)'); clear meg_data_tmp
            end
        end % f
        meg_data = reshape(meg_data,size(meg_data,1),[]);
        
        % Remove effect of time
        if option.remtime
            warning off
            sprintf('......Removing effect of time......')
            times=option.timevector;
            drift_reg = [];
            for n = 1:option.ndrifts
                drift_reg = [drift_reg, (pdist(times').^n)']; % Drift models
            end
            % Remove NaNs
            x = [];
            x = [model, meg_data]; clear meg_data % combine models and data in same matrix
            x(any(isnan(x(:,[1:nmods round(size(x,2)/2)])),2),:) = [];
            for i = nmods+1:size(x,2)
                [~,~,r] = regress(x(:,i),[ones(size(x,1),1) drift_reg]);
                x(:,i) = r;
            end % i
            clear drift_reg r i
            warning on
        else
            % Remove NaNs
            x = [];
            x = [model, meg_data]; clear meg_data % combine models and data in same matrix
            x(any(isnan(x(:,[1:nmods round(size(x,2)/2)])),2),:) = [];
        end
        
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
        r = reshape(r,nmods,length(rdms.ROI_RDMs(1,:,1,1)),length(rdms.ROI_RDMs(:,1,1,1)));
    else % run each model seperately
        for m = 1:nmods
            vs = length(rdms.ROI_RDMs(1,1,:,1));
            meg_data = zeros((((vs*vs)-vs)/2),length(rdms.ROI_RDMs(1,:,1,1)),length(rdms.ROI_RDMs(:,1,1,1)),'single');
            for f = 1:length(rdms.ROI_RDMs(1,:,1,1))
                for time = 1:length(rdms.ROI_RDMs(:,1,1,1))
                    meg_data_tmp = squeeze(rdms.ROI_RDMs(time,f,1:vs,1:vs));  % get MEG RDM for this timepoint
                    meg_data(:,f,time) = single(vectorizeRDM(meg_data_tmp)'); clear meg_data_tmp
                end
            end % f
            meg_data = reshape(meg_data,size(meg_data,1),[]);
            
            % Remove effect of time
            if option.remtime
                warning off
                sprintf('......Removing effect of time......')
                times=option.timevector;
                drift_reg = [];
                for n = 1:option.ndrifts
                    drift_reg = [drift_reg, (pdist(times').^n)']; % Drift models
                end
                % Remove NaNs
                x = [];
                x = [model(:,m), meg_data];
                drift_reg(any(isnan(x(:,[1 round(size(x,2)/2)])),2),:) = [];
                x(any(isnan(x(:,[1 round(size(x,2)/2)])),2),:) = [];
                for i = 2:size(x,2)
                    [~,~,r] = regress(x(:,i),[ones(size(x,1),1) drift_reg]);
                    x(:,i) = r;
                end % i
                clear drift_reg r i
                warning on
            else
                % Remove NaNs
                x = [];
                x = [model(:,m), meg_data]; % combine models and data in same matrix
                x(any(isnan(x(:,[1 round(size(x,2)/2)])),2),:) = []; % only check model and first time for NaNs
            end
            
            % correlate
            sprintf('......model correlations......')
            y=double(x(:,1)); % seperate model and data for memory saving
            x(:,1) = []; x=double(x);
            tmpr = corr(y,x,'type',option.dist,'rows','pairwise');
            r(m,:,:) = single(reshape(tmpr,length(rdms.ROI_RDMs(1,:,1,1)),length(rdms.ROI_RDMs(:,1,1,1)))); clear tmpr x y
        end
        clear meg_data
    end
else
    vs = length(rdms.ROI_RDMs(1,1,:,1));
    meg_data = zeros((((vs*vs)-vs)/2),length(rdms.ROI_RDMs(:,1,1,1)),'single');
    for time = 1:length(rdms.ROI_RDMs(:,1,1))
        meg_data_tmp = squeeze(rdms.ROI_RDMs(time,1:vs,1:vs));  % get MEG RDM for this timepoint
        meg_data(:,time) = single(vectorizeRDM(meg_data_tmp)'); clear meg_data_tmp
    end
    
    if option.parmodels  % run all models together
        
        % Remove effect of time
        if option.remtime
            warning off
            sprintf('......Removing effect of time......')
            times=option.timevector;
            drift_reg = [];
            for n = 1:option.ndrifts
                % Drift models
                drift_reg = [drift_reg, (pdist(times').^n)'];
            end
            % Remove NaNs
            x = [];
            x = [model, meg_data];  clear meg_data
            drift_reg(any(isnan(x(:,[1:nmods round(size(x,2)/2)])),2),:) = [];
            x(any(isnan(x(:,[1:nmods round(size(x,2)/2)])),2),:) = [];
            for i = nmods+1:size(x,2)
                [~,~,r] = regress(x(:,i),[ones(size(x,1),1) drift_reg]);
                x(:,i) = r;
            end % i
            clear drift_reg r i
            warning on
        else
            % Remove NaNs
            x = [];
            x = [model, meg_data]; clear meg_data % combine models and data in same matrix
            x(any(isnan(x(:,[1:nmods round(size(x,2)/2)])),2),:) = [];
        end
        
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
            % Remove effect of time
            if option.remtime
                warning off
                sprintf('......Removing effect of time......')
                times=option.timevector;
                drift_reg = [];
                for n = 1:option.ndrifts
                    % Drift models
                    drift_reg = [drift_reg, (pdist(times').^n)'];
                end
                % Remove NaNs
                x = [];
                x = [model(:,m), meg_data];
                drift_reg(any(isnan(x(:,[1:nmods round(size(x,2)/2)])),2),:) = [];
                x(any(isnan(x(:,[1 round(size(x,2)/2)])),2),:) = [];
                for i = 2:size(x,2)
                    [~,~,r] = regress(x(:,i),[ones(size(x,1),1) drift_reg]);
                    x(:,i) = r;
                end % i
                clear drift_reg r i
                warning on
            else
                % Remove NaNs
                x = [];
                x = [model(:,m), meg_data]; % combine models and data in same matrix
                x(any(isnan(x(:,[1 round(size(x,2)/2)])),2),:) = [];
            end
            
            % correlate
            tmpr = corr(double(x(:,1)),double(x(:,2:end)),'type',option.dist,'rows','pairwise');
            r(m,:) = single(tmpr); clear tmpr
        end
        clear meg_data
    end
end


function [st_rdm] = complexdist(a)
% Distances based on complex numbers. In essense, this is the distance
% between the two end-points of the vectors on the real and imag axis
% Input a is a complex valued matrix

% calc degrees from origin of each vector
degs = 180./(pi./(angle(a)));

% Calculate distances between end-points of vectors
st_rdm = single(zeros(size(a,2),((size(a,1))^2-size(a,1))/2));
ind=0;
for ii = 1:size(a,1)
    for jj = ii:size(a,1)
        if ii == jj
            ii=ii;
        else
            ind = ind+1;
            st_rdm(:,ind) = abs(sqrt(sum(real(a([ii,jj],:)).^2)-2.*real(a(ii,:)).*real(a(jj,:)).*cosd(degs(ii,:)-degs(jj,:))));
        end
    end
end
% sum vector distances over time (gives something like the area between the two sets of vectors)
st_rdm = squeeze(sum(st_rdm,1));
st_rdm = squareform(st_rdm);
