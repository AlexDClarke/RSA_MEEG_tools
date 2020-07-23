function roiRSA_rfxplots_3D_groupFFX(option,s,resetoption)

% Creates timecourse plots for RFX analysis
% Works for '3D' timecourses where models are stacked along 1 dimension
%
% Alex Feb 2015/2019

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
    load([option.rsafront masknam '_groupFFX_' option.midname num2str(option.tw*option.srate) 'ms_sTW.mat'],'rsa_out','option');
    
    mod_names = fieldnames(Models);
    nmods = size(mod_names,1);
    
    if resetoption == 1
    option = optionsfile(s);
    end
    
    for modl = 1:nmods       
        % Calculate mean, st, t-stats
        x = rsa_out.(mod_names{modl});
        x(isnan(x)) = 0;                
        av(modl,:) = x;        
    end
        
    figure; surf(av,'FaceColor','interp',...
        'EdgeColor','none',...
        'FaceLighting','phong')
    axis tight
%   grid off
    view(-25,40)    
    set(gca,'XTick',[1:option.jump*2/option.srate:option.epoch_length/option.srate]);
    set(gca,'XTickLabel',[round(option.start):option.jump*2:round(option.stop)]);
    set(get(gca,'XLabel'),'String','Time (ms)', 'fontweight','b')
    xlim(option.range)
    set(get(gca,'YLabel'),'String','MikeNet Layer', 'fontweight','b')   
    zlim([-0.075 0.15]); %option.ylimits)
    caxis([0 0.1])
    set(get(gca,'zLabel'),'String','Similarity', 'fontweight','b')    
    set(gcf,'PaperPositionMode','auto')
    print('-dtiff',[masknam '_' option.models '_3D_correlation']);
    close
    
%     figure;
%     z = ones(size(tstat))*abs(tinv(option.alpha,size(option.subs,2)-1));
%     surf(z,'FaceColor',[.8 .8 .8],...
%         'EdgeColor','none',...
%         'FaceLighting','phong','FaceAlpha',0.75)    
%     hold on
%     surf(tstat,'FaceColor','interp',...
%         'EdgeColor','none',...
%         'FaceLighting','phong')
%     axis tight
%     view(-25,40)
%     colormap hot
%     set(gca,'XTick',[1:option.jump*2/option.srate:option.epoch_length/option.srate]);
%     set(gca,'XTickLabel',[round(option.start):option.jump*2:round(option.stop)]);
%     set(get(gca,'XLabel'),'String','Time (ms)', 'fontweight','b')
%     xlim(option.range)
%     set(get(gca,'YLabel'),'String','MikeNet Layer', 'fontweight','b')
%     zlim([-2 5])
%     caxis([-2 5])
%     set(get(gca,'zLabel'),'String','T-value', 'fontweight','b')
%     set(gcf,'PaperPositionMode','auto')
%     print('-dtiff',[masknam '_' option.models '_3D_tmap']);
%     close
    
end % mask
end
