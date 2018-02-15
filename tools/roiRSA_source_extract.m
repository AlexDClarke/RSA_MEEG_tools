function roiRSA_source_extract(option,s)

%% Extracts timecourses from vertices within an ROI mask image
%
% Outputs SPM ROI file containing the trial timecourses for each vertex in
% the ROI in the same trial order as the original data (including rejected
% trials)
%
% Alex Feb 2015

spm eeg

val = option.val;

%% Get coordinates of vertices in ROIs
for mask = 1:length(option.masknic)
    
    % Get either single coordinate for list of coordinates
    if option.usemask
        P = ([option.maskdir option.maskname{mask}]);
        [M, XYZa] = spm_read_vols(spm_vol(P));
        XYZ = XYZa(:,find(M))'; clear P M XYZa V
        option.midname = 'mask_';
    else
        XYZ = option.ROI_coords{mask};
        option.midname = 'coord_';
    end        
    
    %% Extract ROI data
    for sub = 1:length(option.subs);
        
        option = optionsfile(s); % need to reset as update below
        mkdir([option.datadir option.sub_beg option.subs{sub} option.subdir]);
        cd([option.datadir option.sub_beg option.subs{sub} option.subdir]);
                      
        % Only extract vertex data if not already done
        if exist([option.masknic{mask} option.midname option.front option.subs{sub} option.lastbit],'file')
        sprintf('Subject %s, Region %s already extracted, will not do again', option.subs{sub},num2str(mask))
        % Load vertex data
        Ds = spm_eeg_load([option.masknic{mask} option.midname option.front option.subs{sub} option.lastbit]);
        else        
        theData = [option.datadir option.sub_beg option.subs{sub} option.subdir1 option.front option.subs{sub} option.lastbit];
        
        % Extract Virtual Channel(s)
        D = spm_eeg_load(theData);
        
        % Get vertices in ROI
        vert  = D.inv{1}.mesh.tess_mni.vert;
        
        if option.usemask
            vert = int16(vert);
            member = ismember(vert,XYZ, 'rows');
            XYZs = double(vert(member,:)); member=[]; vert=[];
        else
            dist = sqrt(sum([vert(:,1) - XYZ(1,1), ...
                vert(:,2) - XYZ(1,2), ...
                vert(:,3) - XYZ(1,3)].^2, 2));
            XYZs = vert(find(dist < option.srad),:);
        end
        
%       epoch = find(abs(D.time - -option.baseline/1000) < 1/10000):find(abs(D.time - (option.epoch_length-option.baseline)/1000) < 1/10000);
        
        D.val = option.val;
        D.inv{val}.source.fname = [option.datadir option.sub_beg option.subs{sub} option.subdir option.masknic{mask} option.midname option.front option.subs{sub} option.lastbit];
        D.inv{val}.source.type  = 'trials';
        D.inv{val}.source.rad = 0;  % This should be zero if using only want 1 vertex
        D.inv{val}.source.XYZ = XYZs;
        D.inv{val}.source.It = [-option.baseline_length option.epoch_length]; % TW to get from ROIs (ms)
        % SPM file you want to extract the data from
        trialfile = [option.datadir option.sub_beg option.subs{sub} option.subdir1 option.trialfront option.subs{sub} option.lastbit];
        
        numClusters = size(XYZs,1);
        for r = 1:numClusters
            D.inv{val}.source.label{r} = sprintf('ROI %d %d %d chan%d',XYZs(r,1), XYZs(r,2), XYZs(r,3), r);  % Insert your names if you want
        end
        
        Ds = spm_eeg_inv_extract_keep_order_ac(D,trialfile);
        D=[];
        end  % extract data again?
                
        if option.nonlocked
            clear Ds
            S.D = [option.datadir option.sub_beg option.subs{sub} option.subdir option.masknic{mask} option.midname option.front option.subs{sub} option.lastbit];
            Ds = spm_eeg_remove_evoked(S)
        end        
        
        %% Select only top 100 verticies based on average zscore
        if option.mask100
            zdata = []; topind = [];
            % Get mean zscore for stats period
            etimes = find(Ds.time < (option.stats_epoch(2))/1000 & Ds.time > (option.stats_epoch(1)/1000))-1;
            btimes = find(Ds.time < (option.bl_ms(2))/1000 & Ds.time > (option.bl_ms(1)/1000))-1;
            for i = 1:size(Ds,1)
                zdata(i) = abs(mean2((squeeze(Ds(i,etimes,:)) - repmat(squeeze(mean(Ds(i,btimes,:),2)),[1 size(etimes,2)])')./repmat(squeeze(std(Ds(i,btimes,:),[],2)),[1 size(etimes,2)])'));
            end
            % Find top 100
            zdata = [zdata; [1:size(Ds,1)]];
            zdata = sortrows(zdata',-1);
            zdata(isnan(zdata(:,1)),:) = [];  % remove rows of NaNs (e.g. a flat timecourse)
            topind = sort(zdata(1:100,2));
        else
            topind = [1:size(Ds,1)]';  % use all
        end               
        
        %% Process extracted ROI data                
        if option.doTF || option.doTFbands
            addpath(genpath('/home/alex/software/EEGLAB/eeglab13_4_4b/'));
            eeglab
            
            % create matfile for partial saving on each loop (remove any that already exist)
            if option.doPhase || option.doPhasePower
                if exist([option.tf_pre 'phz_' option.masknic{mask} option.midname option.front option.subs{sub} '.mat'],'file');
                    delete([option.tf_pre 'phz_' option.masknic{mask} option.midname option.front option.subs{sub} '.mat']);
                end
                tfrs = matfile([option.tf_pre 'phz_' option.masknic{mask} option.midname option.front option.subs{sub}],'Writable',true);
            else
                if exist([option.tf_pre option.masknic{mask} option.midname option.front option.subs{sub} '.mat'],'file');
                    delete([option.tf_pre option.masknic{mask} option.midname option.front option.subs{sub} '.mat']);
                end
                tfrs = matfile([option.tf_pre option.masknic{mask} option.midname option.front option.subs{sub}],'Writable',true);
            end

            %% TF version
            % Work out when TF limits should be
            start = (1000/min(option.fsc))*(round(option.cycles/2));
            stop = (size(Ds,2)*(1000/Ds.fsample)) - (1000/min(option.fsc))*(round(option.cycles/2));
        
            % TF loop for each vertex in ROI
            ei = 0;
            for e = topind'
                ei=ei+1;
                                                               
                [roi_data,freqs,times] = timefreq(squeeze(Ds(e,:,:)),Ds.fsample,'cycles',[option.cycles option.cycles*3],'freqs',[option.fsc(1) option.fsc(end)],'nfreqs',option.nfs,'freqscale','log','timesout',[start:option.tfstep:stop]);   % currently set to match prev analysis I did
                
                % Update time info
                if ei == 1
                    D.times = (round(times/10)*10) - abs(Ds.timeonset*1000);  %round(times - abs(Ds.timeonset*1000),-1); % correct times output, round to nearest 10
                end
                
                if option.doPhase
                    % Get phases
                    roi_data  = angle(roi_data);
                
                elseif option.doPhasePower
                    roi_data  = roi_data;  % keep the complex data                                       
                else
                    % Convert to power
                    if option.cycles == 0
                        roi_data = 2/0.375*roi_data/max(pow2(nextpow2(size(D.data,2))-3),4);
                    else
                        roi_data  = roi_data.*conj(roi_data);
                    end
                    roi_data = 10*log10(roi_data);
                                       
                    % Baseline correct (zscore) each freq
                    btimes = find(D.times < option.bl_ms(2) & D.times > option.bl_ms(1));
                    for i = 1:length(roi_data(:,1,1))
                        roi_data(i,:,:) = (squeeze(roi_data(i,:,:)) - repmat(squeeze(mean(roi_data(i,btimes,:),2)),[1 size(roi_data,2)])')./repmat(squeeze(std(roi_data(i,btimes,:),[],2)),[1 size(roi_data,2)])';
                    end
                end
                
                temp_o(1,:,:,:) = single(roi_data);
                tfrs.D(ei,1:length(roi_data(:,1,1)),1:length(roi_data(1,:,1)),1:length(roi_data(1,1,:))) = temp_o;
                clear roi_data btimes temp_o
                
            end % electrode loop
            
            times = D.times;
                        
            % Save data
            % Update time info with new time resolution
            option.epoch_length = max(times) + abs(min(times));
            option.baseline_length = abs(min(times));
            option.tfepoch = times;
            option.baseline = option.baseline_length/option.tfstep;
            option.fsc = freqs;  % just in case different
            
             % Save list of trials to reject
            reject = zeros(length(Ds.conditions),1);
            reject(Ds.badtrials) = 1;
            save([option.datadir option.sub_beg option.subs{sub} option.subdir option.trials option.subs{sub} '.mat'],'reject');
            
            tfrs.option = option;                        
            Ds=[]; S=[];
            
        else
            %% Standard data
                 
            if option.zscore
                % Zscore data using baseline
                S.D = Ds.fname;
                S.timewin = option.bl_ms;
                Ds = spm_eeg_zscore(S);
            end
            
            % Save list of trials to reject
            reject = zeros(length(Ds.conditions),1);
            reject(Ds.badtrials) = 1;
            save([option.datadir option.sub_beg option.subs{sub} option.subdir option.trials option.subs{sub} '.mat'],'reject');
            
            % Save ROI data for RSA processing
            D = Ds(topind,:,:);            
            outfile = [option.masknic{mask} option.midname option.front option.subs{sub}];
            save(outfile,'D','option');
            
            D=[]; Ds=[]; S=[];
        end
    end % sub
    
end % mask
