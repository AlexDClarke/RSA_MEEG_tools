function roiRSA_noise_ceilings(option,s)

% Correlates each subjects data RDMs with group_average RDM to get noice
% ceilings
%
% AC Feb 2017

% Begin
for mask = 1:length(option.masknic)
    
    masknam = [option.masknic{mask}];           
    
    %% Collect subject matfile names
    for sub = 1:length(option.subs)
        
        cd([option.datadir option.sub_beg option.subs{sub} option.subdir]);
                
        % Load data RDMs
        if option.doTF
            if sub == 1
            option.tw = round((option.tw*option.srate)/option.tfstep);
            end
            if option.tw > 1
                if mod(option.tw,2)  % must be even number
                    option.tw = option.tw+1;
                end
            end
            if option.doavg
                if option.doPhase
                    rdms.(['sub_' num2str(sub)]) = matfile([option.tf_pre 'phz_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                else
                    rdms.(['sub_' num2str(sub)]) = matfile([option.tf_pre option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                end
            else
                if option.doPhase
                    rdms.(['sub_' num2str(sub)]) = matfile([option.tf_pre 'phz_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                else
                    rdms.(['sub_' num2str(sub)]) = matfile([option.tf_pre option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                end
            end
        elseif option.doTFbands
            if sub == 1
            if option.doTWs
                option.tw = 1;
            else
                option.tw = round((option.tw*option.srate)/option.tfstep);
            end
            end
            if option.tw > 1
                if mod(option.tw,2)  % must be even number
                    option.tw = option.tw+1;
                end
            end
            if option.doavg
                if option.doPhase
                    rdms.(['sub_' num2str(sub)]) = matfile([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                else
                    rdms.(['sub_' num2str(sub)]) = matfile([option.tf_pre 'bands_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                end
            else
                if option.doPhase
                    if option.doTWs
                        rdms.(['sub_' num2str(sub)]) = matfile([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'TW_RDMs.mat']);
                    else
                        rdms.(['sub_' num2str(sub)]) = matfile([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                    end
                else
                    if option.doTWs
                        rdms.(['sub_' num2str(sub)]) = matfile([option.tf_pre 'bands_' option.masknic{mask} option.midname 'TW_RDMs.mat']);
                    else
                        rdms.(['sub_' num2str(sub)]) = matfile([option.tf_pre 'bands_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
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
                rdms.(['sub_' num2str(sub)]) = matfile([option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.srate) 'ms_sTW.mat']);
            else
                rdms.(['sub_' num2str(sub)]) = matfile([option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.srate) 'ms_sTW.mat']);
            end
        end
        
    end  % sub
     
    %% Calculate Noice ceilings
    for time = 1:length(rdms.(['sub_1']).ROI_RDMs(:,1,1))
        
        vs = length(rdms.(['sub_' num2str(sub)]).ROI_RDMs(1,1,:,1));
        meg_data = zeros(1,(((vs*vs)-vs)/2),length(option.subs),'single');
        for sub = 1:length(option.subs)
            meg_data_tmp = squeeze(rdms.(['sub_' num2str(sub)]).ROI_RDMs(time,1:vs,1:vs));  % get MEG RDM for this timepoint
            meg_data(1,:,sub) = single(vectorizeRDM(meg_data_tmp)'); clear meg_data_tmp
        end
                
        meg_data(:,:,option.sub_out) = []; % exclude rejected subjects
        
        % Use adapted RSA toolbox function
        [ceiling_upperBound, ceiling_lowerBound, RDMcorrs_LOO, bestFitRDM] = ceilingAvgRDMcorr_ac(meg_data,option.dist);

        % Store results
        nc_out.upper(time) = ceiling_upperBound;
        nc_out.lower(time) = ceiling_lowerBound;
        nc_out.sub_lower(time,:) = RDMcorrs_LOO;
        
    end
    
    %% Plots  
    rfxdir = option.rfxdir;
    if ~exist(rfxdir)
        mkdir(rfxdir)
    end        
    cd(rfxdir);
    
    % Bounds
    figure;
    Color = [1 0 0];
    hold on
    plot(zeros(1,length(nc_out.upper)),'LineWidth',1,'Color',[0 0 0])
    plot(nc_out.upper,'LineWidth',2);
    plot(nc_out.lower,'LineWidth',2);    
    set(gca,'XTick',[1:option.jump/option.srate:option.epoch_length/option.srate]);
    set(gca,'XTickLabel',[round(option.start):option.jump:round(option.stop)]);
    set(get(gca,'XLabel'),'String','Time (ms)', 'fontweight','b')
    set(get(gca,'YLabel'),'String','Noice ceilings', 'fontweight','b')
    xlim(option.range)
    hold off
    set(gcf,'PaperPositionMode','auto')
    print('-dtiff',[masknam '_noise_ceiling_bounds']);
    close
    
    % subject lower bounds
    figure;
    Color = [1 0 0];
    hold on
    plot(zeros(1,length(nc_out.upper)),'LineWidth',1,'Color',[0 0 0])
    plot(nc_out.sub_lower,'LineWidth',2);
    set(gca,'XTick',[1:option.jump/option.srate:option.epoch_length/option.srate]);
    set(gca,'XTickLabel',[round(option.start):option.jump:round(option.stop)]);
    set(get(gca,'XLabel'),'String','Time (ms)', 'fontweight','b')
    set(get(gca,'YLabel'),'String','Noice ceilings', 'fontweight','b')
    xlim(option.range)
    hold off
    set(gcf,'PaperPositionMode','auto')
    print('-dtiff',[masknam '_noise_ceiling_subjects']);
    close
    
    %% Save outputs                   
    if option.doTF                    
        if option.doPhase
            outname = [option.rfxdir option.tf_pre 'phz_' option.noisefront option.masknic{mask} option.midname num2str(option.tw*option.tfstep) 'ms_sTW.mat'];
            save(outname, 'nc_out', 'option'); clear nc_out
        else
            outname = [option.rfxdir option.tf_pre option.noisefront option.masknic{mask} option.midname num2str(option.tw*option.tfstep) 'ms_sTW.mat'];
            save(outname, 'nc_out', 'option'); clear nc_out
        end
    elseif option.doTFbands        
        if option.doTWs
            if option.doPhase
                outname = [option.rfxdir option.tf_pre 'bands_phz_' option.noisefront option.masknic{mask} option.midname '.mat'];
                save(outname, 'nc_out', 'option'); clear nc_out
            else
                outname = [option.rfxdir option.tf_pre 'bands_' option.noisefront option.masknic{mask} option.midname '.mat'];
                save(outname, 'nc_out', 'option'); clear nc_out
            end
        else
            if option.doPhase
                outname = [option.rfxdir option.tf_pre 'bands_phz_' option.noisefront option.masknic{mask} option.midname num2str(option.tw*option.tfstep) 'ms_sTW.mat'];
                save(outname, 'nc_out', 'option'); clear nc_out
            else
                outname = [option.rfxdir option.tf_pre 'bands_' option.noisefront option.masknic{mask} option.midname num2str(option.tw*option.tfstep) 'ms_sTW.mat'];
                save(outname, 'nc_out', 'option'); clear nc_out
            end
        end
    else        
        outname = [option.rfxdir option.noisefront option.masknic{mask} option.midname num2str(option.tw*option.srate) 'ms_sTW.mat'];
        save(outname, 'nc_out', 'option'); clear nc_out
    end
    
end % masknic
