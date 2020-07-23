function roiRSA_ecog_extract_alt(option,s)

%% Extracts timecourses for RSA from 3D ECOG file
%

%spm eeg

m=load('/work/imaging3/MEG/AC_ECOG_objects/Localizations/electrode_roi_locations_mono.mat');

%% Get coordinates of vertices in ROIs
for mask = 1:length(option.masknic)
    
    mask_elecs = m.roi_elecs(:,strcmp(m.ROIs,option.masknic{mask}));
    
    option.midname = 'ecog_';
        
    %% Extract ROI data
    for sub = 1:length(option.subs);
        
%        % Is there data for this mask/subject? If not skip to next
%        if isempty(mask_elecs{sub});
%            sprintf('......Subject %s has no electrodes in region %s......', option.subs{sub},option.masknic{mask})
%            continue;
%        end;
        
        option = optionsfile(s); % need to reset as update below
        if ~exist([option.datadir option.sub_beg option.subs{sub} option.subdir],'dir');
            mkdir([option.datadir option.sub_beg option.subs{sub} option.subdir]);
        end
        cd([option.datadir option.sub_beg option.subs{sub} option.subdir]);
        
        theData = [option.datadir option.sub_beg option.subs{sub} option.subdir1 option.trialfront option.subs{sub} option.lastbit];
        % Open data file
        load(theData);
        Ds = concepts_all;
        clear concepts_all
        
        % specify data for ROI
        topind = mask_elecs{find(strcmp(m.subs,option.subs{sub}))};
%       topind = mask_elecs{sub};
        
        if isempty(mask_elecs{find(strcmp(m.subs,option.subs{sub}))});
            sprintf('......Subject %s has no good electrodes in region %s......', option.subs{sub},option.masknic{mask})
            break;
        end;
               
        
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

            %% Get precomputed TF
        
            % Save list of trials to reject
            reject = squeeze(isnan(Ds(1,1,1,:)));
            save([option.datadir option.sub_beg option.subs{sub} option.subdir option.trials option.subs{sub} '.mat'],'reject');            

            % Replace trial NaNs with zeros
            Ds(:,:,:,isnan(Ds(1,1,1,:))) = 0;
            
            % TF loop for each vertex in ROI
            ei = 0;
            for e = topind
                ei=ei+1;
                
                tfrs.D(ei,1:length(Ds(1,:,1,1)),1:length(Ds(1,1,:,1)),1:length(Ds(1,1,1,:))) = single(Ds(e,:,:,:));
                       
            end % electrode loop
            
            times = o.times;
                        
            % Save data
            % Update time info with new time resolution
            option.epoch_length = max(times) + abs(min(times));
            option.baseline_length = abs(min(times));
            option.tfepoch = times;
            option.baseline = option.baseline_length/option.tfstep;
            option.fsc = o.freqs;  % just in case different
                        
            tfrs.option = option;                        
            Ds=[];
            
        else
            %% Standard data
            
            % Limit to electrodes we want
            Ds = Ds(topind,:,:); 
            
            % Baseline correct (zscore)
            btimes = find(o.times < option.bl_ms(2) & o.times > option.bl_ms(1));
            for i = 1:length(Ds(:,1,1))
                Ds(i,:,:) = (squeeze(Ds(i,:,:)) - repmat(squeeze(mean(Ds(i,btimes,:),2)),[1 size(Ds,2)])')./repmat(squeeze(std(Ds(i,btimes,:),[],2)),[1 size(Ds,2)])';
            end

            % Save list of trials to reject
            reject = squeeze(isnan(Ds(1,1,:)));  % find NaN trials that are missing
%           reject = zeros(length(Ds.conditions),1);
%           reject(Ds.badtrials) = 1;
            save([option.datadir option.sub_beg option.subs{sub} option.subdir option.trials option.subs{sub} '.mat'],'reject');
            
            % Save ROI data for RSA processing
            D = Ds;
            D(:,:,isnan(Ds(1,1,:))) = 0; % convert trial nans to zeros
            outfile = [option.masknic{mask} option.midname option.front option.subs{sub}];
            save(outfile,'D','option');
            
            D=[]; Ds=[]; S=[];
        end
    end % sub
    
end % mask
