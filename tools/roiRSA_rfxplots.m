function roiRSA_rfxplots(option,s,resetoption)

% Creates timecourse plots for RFX analysis
% Works for timecourses, TFR, and TF band plots
%
% Alex Feb 2015

cd(option.rfxdir);
load([option.var_dir option.models '.mat']);  % Load to get names
alpha = option.alpha;

if ~or(option.doTF,option.doTFbands)
%% Standard timecourses
for mask = 1:length(option.masknic)
    
    masknam = [option.masknic{mask}];
    if option.tw > 1
    if mod(option.tw,2)  % must be even number
        option.tw = option.tw+1;
    end
    end
    
    % Load data
    load([option.rsafront masknam option.midname num2str(option.tw*option.srate) 'ms_sTW.mat'],'rsa_out','option');
    
    mod_names = fieldnames(Models);
    nmods = size(mod_names,1);
    
    if resetoption == 1
    option = optionsfile(s);
    end
    
    for modl = 1:nmods
        
        % Calculate mean, st, t-stats
        x = rsa_out.(mod_names{modl});
        x(isnan(x)) = 0;
        x(option.sub_out,:) = [];
        
        [h p ci stats] = ttest(x);  % get t-value timecourse
        av = mean(x);
        rsa_std = std(x,[],1);
        
        % Spearmans plot
        figure;
        Color = [1 0 0];
        hold on
        plot(zeros(1,length(av)),'LineWidth',1,'Color',[0 0 0])
        se = rsa_std(1,:)/sqrt(size(x,1));        
        h1 = boundedline([1:length(av)],av,se,'alpha','cmap',Color);        
        % For matlab 2014b onwards
%         h1.LineWidth = 2;
%         ax=gca;
%         ax.XTick = [1:option.jump/option.srate:option.epoch_length/option.srate];
%         ax.XTickLabel = [round(option.start):option.jump:round(option.stop)];        
%         ax.XLabel.String = 'Time (ms)';
%         ax.XLabel.FontWeight = 'b';
%         ax.YLabel.String = 'Similarity';
%         ax.YLabel.FontWeight = 'b';
%         ax.XLim = option.range;        
%         % Pre-2014b
        set(h1,'LineWidth',2);        
        set(gca,'XTick',[1:option.jump/option.srate:option.epoch_length/option.srate]);
        set(gca,'XTickLabel',[round(option.start):option.jump:round(option.stop)]);
        set(get(gca,'XLabel'),'String','Time (ms)', 'fontweight','b')
%       ylim(option.ylimits)
        set(get(gca,'YLabel'),'String','Similarity', 'fontweight','b')
        xlim(option.range)        
        hold off
        set(gcf,'PaperPositionMode','auto')
        print('-dtiff',[masknam mod_names{modl} '_correlation']);
        close
        
        % Stats plot
        figure;
        Color = [0 0.75 0.25];
        hold on
        plot(zeros(1,length(av)),'LineWidth',1,'Color',[0 0 0])
        plot(ones(1,length(av)).*abs(tinv(alpha,size(x,1)-1)),'--','LineWidth',0.75,'Color',[0 0 0])
        plot(stats.tstat,'-','LineWidth',2,'Color',Color)
        % For matlab 2014b onwards
%         ax=gca;
%         ax.XTick = [1:option.jump/option.srate:option.epoch_length/option.srate];
%         ax.XTickLabel = [round(option.start):option.jump:round(option.stop)];        
%         ax.XLabel.String = 'Time (ms)';
%         ax.XLabel.FontWeight = 'b';
%         ax.YLabel.String = 'Tscore';
%         ax.YLabel.FontWeight = 'b';
%         ax.XLim = option.range;         
%         % Pre-2014B
        set(gca,'XTick',[1:option.jump/option.srate:option.epoch_length/option.srate]);
        set(gca,'XTickLabel',[round(option.start):option.jump:round(option.stop)]);        
        set(get(gca,'XLabel'),'String','Time (ms)', 'fontweight','b')
        set(get(gca,'YLabel'),'String','Tscore', 'fontweight','b')
        xlim(option.range)        
        hold off
        set(gcf,'PaperPositionMode','auto')
        print('-dtiff',[masknam mod_names{modl} '_tstats']);
        close                
        
    end % modl
    
end % mask
end

if option.doTF && option.doTFbands==0
    
    if resetoption == 1
        option = optionsfile(s);
    end
    
%% TFR plots
for mask = 1:length(option.masknic)
    option = optionsfile(s);    
    masknam = [option.masknic{mask}];        
    
%     if option.tw > 1
%         if mod(option.tw,2)  % must be even number
%             option.tw = option.tw+1;
%         end
%     end   
    if (option.tw*option.srate)/option.tfstep > 1
        if mod((option.tw*option.srate)/option.tfstep,2)  % must be even number
            tfw = round((option.tw*option.srate)/option.tfstep)+1;
        else
            tfw = round((option.tw*option.srate)/option.tfstep);
        end
    else
        tfw = 1;
    end
    
    % Load data
    if isfield(option,'tf_pre2')
        prenam = [option.tf_pre option.tf_pre2];        
    else
        prenam = [option.tf_pre];
    end
    
    if option.doPhase || option.doPhasePower % phase
        prenam = [prenam 'phz_'];    
    end
        
    if option.tw == 1
        load([prenam option.rsafront masknam option.midname num2str(0*option.tfstep) 'ms_sTW.mat'],'rsa_out','option');
    else
        load([prenam option.rsafront masknam option.midname num2str(tfw*option.tfstep) 'ms_sTW.mat'],'rsa_out','option');
    end
    
    mod_names = fieldnames(Models);
    nmods = size(mod_names,1);
        
    if resetoption == 1
        options = optionsfile(s,1);        
        option = setfield(option,'start',options.start);
        option = setfield(option,'stop',options.stop);
        option = setfield(option,'bl_ms',options.bl_ms);
        option = setfield(option,'stats_epoch',options.stats_epoch);
        option = setfield(option,'ylimits',options.ylimits);
        option = setfield(option,'tlimits',options.tlimits);
        option = setfield(option,'jump',options.jump);
        option = setfield(option,'sub_out',options.sub_out); clear options
    end        
    freqs = option.fsc;
    
    option.range = [min(find(option.tfepoch > option.bl_ms(1))) min(find(option.tfepoch > option.stats_epoch(2)))];
    
    for modl = 1:nmods
        
        % Calculate mean, st, t-stats
        x = rsa_out.(mod_names{modl});
        x(isnan(x)) = 0;
        x(option.sub_out,:,:) = [];
        
        [h p ci stats] = ttest(x);  % get t-value timecourse
        av = squeeze(mean(x));
        
        % Spearmans plot
        figure;
        imagesc(flipdim(av,1));
        hold on
        mark = min(find(abs(option.tfepoch) == min(abs(0-option.tfepoch))));
        line([mark mark],get(gca,'Ylim'),'Color',[0 0 0],'LineWidth',2); hold off
        
        %% Fix time axis
        aa = flipdim(round(freqs),2);
        set(gca,'XTick',[1:option.jump/option.tfstep:(option.epoch_length-1)/option.tfstep]);
        set(gca,'XTickLabel',[-option.baseline_length:option.jump:(option.epoch_length-option.baseline_length-1)]);
        set(get(gca,'XLabel'),'String','Time (ms)', 'fontweight','b')               
        set(gca,'YTick',[1:2:size(freqs,2)]);
        set(gca,'YTickLabel',aa(1:2:size(freqs,2)));
        set(get(gca,'YLabel'),'String','Freq (Hz)', 'fontweight','b')
        xlim(option.range)
        caxis([option.ylimits])
        colorbar; colormap jet
        set(gcf,'PaperPositionMode','auto')
        print('-dtiff',[prenam masknam mod_names{modl} '_correlation']);
        close
        
        % Stats plot        
        figure;
        imagesc(flipdim(squeeze(stats.tstat),1));
        hold on
        mark = min(find(abs(option.tfepoch) == min(abs(0-option.tfepoch))));
        line([mark mark],get(gca,'Ylim'),'Color',[0 0 0],'LineWidth',2); hold off
        
        set(gca,'XTick',[1:option.jump/option.tfstep:(option.epoch_length-1)/option.tfstep]);
        set(gca,'XTickLabel',[-option.baseline_length:option.jump:(option.epoch_length-option.baseline_length-1)]);
        set(get(gca,'XLabel'),'String','Time (ms)', 'fontweight','b')       
        set(gca,'YTick',[1:2:size(freqs,2)]);
        set(gca,'YTickLabel',aa(1:2:size(freqs,2)));
        set(get(gca,'YLabel'),'String','Freq (Hz)', 'fontweight','b')
        xlim(option.range)
        caxis([option.tlimits])
        colorbar; colormap jet
        set(gcf,'PaperPositionMode','auto')
        print('-dtiff',[prenam masknam mod_names{modl} '_tmap']);
        close
        
        % Stats plot 2                
        if option.tails == 2;
            xx = flipdim(or(squeeze(stats.tstat) < (tinv(alpha/2,size(x,1)-1)), ...
            squeeze(stats.tstat) > abs(tinv(alpha/2,size(x,1)-1))),1);           % plots both tails
        else                        
            xx = flipdim(squeeze(stats.tstat) > abs(tinv(alpha,size(x,1)-1)),1);
        end
        
        figure;
        imagesc(flipdim(av,1),'AlphaData',0.5);
        hold on
        imagesc(flipdim(av,1),'AlphaData',1*(xx));
        contour(xx,1,'black','LineWidth',2);
        
        mark = min(find(abs(option.tfepoch) == min(abs(0-option.tfepoch))));                
        set(gca,'XTick',[1:option.jump/option.tfstep:(option.epoch_length-1)/option.tfstep]);
        set(gca,'XTickLabel',[-option.baseline_length:option.jump:(option.epoch_length-option.baseline_length-1)]);
        set(get(gca,'XLabel'),'String','Time (ms)', 'fontweight','b')                       
        set(gca,'YTick',[1:2:size(freqs,2)]);
        set(gca,'YTickLabel',aa(1:2:size(freqs,2)));
        set(get(gca,'YLabel'),'String','Freq (Hz)', 'fontweight','b')
        xlim(option.range)
        caxis([option.ylimits])
        colorbar; colormap jet        
%         for n=1:size(xx,2)
%             for m = 1:size(xx,1)
%                 if(xx(m,n) == 0)
%                     rectangle('Position',[n-.5 m-.5 1 1],'FaceColor',[.4 .4 .4],'LineStyle','none');
%                 end
%             end
%         end
        line([mark mark],get(gca,'Ylim'),'Color',[0 0 0],'LineWidth',2);    
        hold off
        
        set(gcf,'PaperPositionMode','auto')
        print('-dtiff',[prenam masknam mod_names{modl} '_correlation_p' num2str(alpha) '.tif']);
        close             
        
        % Stats plot 3      
        if option.tails == 2;
        xx =squeeze(stats.tstat) .* or(squeeze(stats.tstat) < (tinv(alpha/2,size(x,1)-1)), ...
            squeeze(stats.tstat) > abs(tinv(alpha/2,size(x,1)-1)));           % plots both tails
        else
        xx =squeeze(stats.tstat) .* (squeeze(stats.tstat) > abs(tinv(alpha,size(x,1)-1)));  % only positive effects
        end
        figure;
        imagesc(flipdim(xx,1));
        hold on
        mark = min(find(abs(option.tfepoch) == min(abs(0-option.tfepoch))));
        line([mark mark],get(gca,'Ylim'),'Color',[0 0 0],'LineWidth',2); hold off
        
        set(gca,'XTick',[1:option.jump/option.tfstep:(option.epoch_length-1)/option.tfstep]);
        set(gca,'XTickLabel',[-option.baseline_length:option.jump:(option.epoch_length-option.baseline_length-1)]);
        set(get(gca,'XLabel'),'String','Time (ms)', 'fontweight','b')                       
        set(gca,'YTick',[1:2:size(freqs,2)]);
        set(gca,'YTickLabel',aa(1:2:size(freqs,2)));
        set(get(gca,'YLabel'),'String','Freq (Hz)', 'fontweight','b')
        xlim(option.range)
        caxis([option.tlimits])
        colorbar; colormap jet
        set(gcf,'PaperPositionMode','auto')
        print('-dtiff',[prenam masknam mod_names{modl} '_tmap_p' num2str(alpha) '.tif']);
        close      
        
        
    end % modl
    
end % mask
end

if option.doTFbands
%% TF bands
scrsz = get(0,'ScreenSize');

% if option.tw > 1
%     if mod(option.tw,2)  % must be even number
%         option.tw = option.tw+1;
%     end
% end
if (option.tw*option.srate)/option.tfstep > 1
    if mod((option.tw*option.srate)/option.tfstep,2)  % must be even number
        tfw = round((option.tw*option.srate)/option.tfstep)+1;
    else
        tfw = round((option.tw*option.srate)/option.tfstep);
    end
else
    tfw = 1;
end

for mask = 1:length(option.masknic)

    masknam = [option.masknic{mask}];
        
    
    % Load data
    prenam = [];   
    if isfield(option,'tf_pre2')
        prenam = [option.tf_pre option.tf_pre2 'bands_'];        
    else
        prenam = [option.tf_pre 'bands_'];
    end
    
    if option.doPhase || option.doPhasePower % phase
        prenam = [prenam 'phz_'];    
    end
        
    if option.tw == 1
        load([prenam option.rsafront masknam option.midname num2str(0*option.tfstep) 'ms_sTW.mat'],'rsa_out','option');
    else
        load([prenam option.rsafront masknam option.midname num2str(tfw*option.tfstep) 'ms_sTW.mat'],'rsa_out','option');
    end
    
    mod_names = fieldnames(Models);
    nmods = size(mod_names,1);
    
    if resetoption == 1
        options = optionsfile(s,1);        
        option = setfield(option,'start',options.start);
        option = setfield(option,'stop',options.stop);
        option = setfield(option,'bl_ms',options.bl_ms);
        option = setfield(option,'stats_epoch',options.stats_epoch);
        option = setfield(option,'ylimits',options.ylimits);
        option = setfield(option,'tlimits',options.tlimits);
        option = setfield(option,'jump',options.jump);
        option = setfield(option,'tw',options.tw);
        clear options
    end        
    freqs = option.fs;
    
    for modl = 1:nmods
        
        % Calculate mean, st, t-stats
        x = rsa_out.(mod_names{modl});
        x(isnan(x)) = 0;
        x(option.sub_out,:,:) = [];        
        [h p ci stats] = ttest(x);  % get t-value timecourse
        av = squeeze(mean(x));
        rsa_std = squeeze(std(x,[],1));
        
        % Spearmans plot
        figure1 = figure('Position',[1 10 scrsz(3)/2 scrsz(4)/1.5]);
        for bands = 1:length(freqs(1,:))            
            if floor(length(freqs(1,:))/2) == (length(freqs(1,:))/2)
                subplot((length(freqs(1,:))/2),(length(freqs(1,:))/2),bands);
            else                
                subplot(ceil(length(freqs(1,:))/2),floor(length(freqs(1,:))/2),bands);
            end            
            title([num2str(freqs{bands}(1)) 'Hz to ' num2str(freqs{bands}(2)) 'Hz']);
            hold on
            Color = [1 0 0];
            plot(zeros(1,length(av)),'LineWidth',1,'Color',[0 0 0])                        
            se = rsa_std(bands,:)/sqrt(size(x,1));
            h1 = boundedline([1:length(av(bands,:))],av(bands,:),se,'alpha','cmap',Color);        
%             % For matlab 2014b onwards
%             h1.LineWidth = 2;
%             ax=gca;
%             ax.XTick = [1:option.jump/option.tfstep:option.epoch_length/option.tfstep];
%             ax.XTickLabel = [option.tfepoch(1):option.jump:option.tfepoch(end)];
%             ax.XLabel.String = 'Time (ms)';
%             ax.XLabel.FontWeight = 'b';
%             ax.YLabel.String = 'Similarity';
%             ax.YLabel.FontWeight = 'b';
%             ax.XLim = [1 length(option.tfepoch)];  % should update to use times within a time vector
%             ax.YLim = [-(max(max(abs(av)))+(0.2*max(max(abs(av))))) max(max(abs(av)))+(0.2*max(max(abs(av))))];
            % Pre-2014b                     
            set(h1,'LineWidth',2);            
            set(gca,'XTick',[1:option.jump/option.tfstep:(option.epoch_length-1)/option.tfstep]);
            set(gca,'XTickLabel',[-option.baseline_length:option.jump:(option.epoch_length-option.baseline_length-1)]);
            set(get(gca,'XLabel'),'String','Time (ms)', 'fontweight','b')            
%           ylim([option.ylimits])
            set(get(gca,'YLabel'),'String','Spearman rho', 'fontweight','b')
            xlim([1 length(option.tfepoch)]);
            hold off
        end
%        set(gcf,'PaperPositionMode','auto')
        print('-dtiff',[prenam masknam mod_names{modl} '_correlation']);
        close
                
        % Stats plot        
        figure1 = figure('Position',[1 10 scrsz(3)/2 scrsz(4)/1.5]);
        for bands = 1:length(freqs(1,:))
            if floor(length(freqs(1,:))/2) == (length(freqs(1,:))/2)
                subplot((length(freqs(1,:))/2),(length(freqs(1,:))/2),bands);
            else                
                subplot(ceil(length(freqs(1,:))/2),floor(length(freqs(1,:))/2),bands);
            end              
            title([num2str(freqs{bands}(1)) 'Hz to ' num2str(freqs{bands}(2)) 'Hz']);         
            Color = [0 0.75 0.25];
            hold on
            plot(zeros(1,length(av)),'LineWidth',1,'Color',[0 0 0])
            plot(ones(1,length(av)).*abs(tinv(alpha,size(x,1)-1)),'--','LineWidth',0.75,'Color',[0 0 0])
            plot(ones(1,length(av)).*(tinv(alpha,size(x,1)-1)),'--','LineWidth',0.75,'Color',[0 0 0])            
            plot(squeeze(stats.tstat(1,bands,:)),'-','LineWidth',2,'Color',Color)
%             % For matlab 2014b onwards
%             ax=gca;
%             ax.XTick = [1:option.jump/option.tfstep:option.epoch_length/option.tfstep];
%             ax.XTickLabel = [option.tfepoch(1):option.jump:option.tfepoch(end)];
%             ax.XLabel.String = 'Time (ms)';
%             ax.XLabel.FontWeight = 'b';
%             ax.YLabel.String = 'Tscore';
%             ax.YLabel.FontWeight = 'b';
%             ax.XLim = [1 length(option.tfepoch)];  % should update to use times within a time vector
%             ax.YLim = [-(max(max(abs(stats.tstat(1,:,:))))+(0.2*max(max(abs(stats.tstat(1,:,:)))))) max(max(abs(stats.tstat(1,:,:))))+(0.2*max(max(abs(stats.tstat(1,:,:)))))];
            % Pre-2014b                     
            set(gca,'LineWidth',2);            
            set(gca,'XTick',[1:option.jump/option.tfstep:(option.epoch_length-1)/option.tfstep]);
            set(gca,'XTickLabel',[-option.baseline_length:option.jump:(option.epoch_length-option.baseline_length-1)]);
            set(get(gca,'XLabel'),'String','Time (ms)', 'fontweight','b')                             
%           ylim([option.ylimits])
            set(get(gca,'YLabel'),'String','Tscore', 'fontweight','b')
            xlim([1 length(option.tfepoch)]);           
            hold off
        end
%        set(gca,'PaperPositionMode','auto')
        print('-dtiff',[prenam masknam mod_names{modl} '_tstats']);
        close
                        
    end % modl
    
end % mask
end
