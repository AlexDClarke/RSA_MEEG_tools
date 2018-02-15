function roiRSA_sensor_extract(option,s)

%% Extracts timecourses SPM sensor file
%

spm eeg

%% Get coordinates of vertices in ROIs
for mask = 1:length(option.masknic)
    
    option.midname = 'sensor_';
        
    %% Extract ROI data
    for sub = 1:length(option.subs);
        
        option = optionsfile(s); % need to reset as update below
        mkdir([option.datadir option.sub_beg option.subs{sub} option.subdir]);
        cd([option.datadir option.sub_beg option.subs{sub} option.subdir]);
                      
        % Only extract vertex data if not already done
        if exist([option.masknic{mask} option.midname option.trialfront option.subs{sub} option.lastbit],'file')
        sprintf('Subject %s, Region %s already extracted, will not do again', option.subs{sub},num2str(mask))
        % Load vertex data
        Ds = spm_eeg_load([option.masknic{mask} option.midname option.trialfront option.subs{sub} option.lastbit]);
        else
        theData = [option.datadir option.sub_beg option.subs{sub} option.subdir1 option.trialfront option.subs{sub} option.lastbit];
        
        % Open data file
        Ds = spm_eeg_load(theData);
                
        end  % extract data again?
                
        % Select data from sensor type
        if strcmp(option.masknic{mask},'EEG')
            topind = setxor(find(strcmp(Ds.chantype,'EEG')),Ds.badchannels);
        else
            topind = find(strcmp(Ds.chantype,option.masknic{mask}));
        end

        %% Process extracted data                
        if option.doTF || option.doTFbands
            addpath(genpath('/home/alex/software/EEGLAB/eeglab13_4_4b/'));
            eeglab
            
            % create matfile for partial saving on each loop (remove any that already exist)
            if option.doPhase
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
            for e = topind
                ei=ei+1;
                                                               
                [roi_data,freqs,times] = timefreq(squeeze(Ds(e,:,:)),Ds.fsample,'cycles',[option.cycles option.cycles*3],'freqs',[option.fsc(1) option.fsc(end)],'nfreqs',option.nfs,'freqscale','log','timesout',[start:option.tfstep:stop]);   % currently set to match prev analysis I did
                
                % Update time info
                if ei == 1
                    D.times = (round(times/10)*10) - abs(Ds.timeonset*1000);  %round(times - abs(Ds.timeonset*1000),-1); % correct times output, round to nearest 10
                end
                
                if option.doPhase
                    % Get phases
                    roi_data  = angle(roi_data);
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
                        
            % Zscore data using baseline
            S.D = theData;
            S.time = option.bl_ms;
            Ds = spm_eeg_zscore(S);
            
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
