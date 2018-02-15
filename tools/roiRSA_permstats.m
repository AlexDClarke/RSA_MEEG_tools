function roiRSA_permstats(option)

% Generates permutation based cluster stats for the RSA timecourses
%
% Alex Feb 2015

% Clear values and set alpha
max_perm_hist = [];
clusters = [];
alpha = option.alpha;
cd(option.rfxdir);
rand_perms = round(rand(length(option.subs)-length(option.sub_out),option.perms));

load([option.var_dir option.models '.mat']);  % Models
mod_names = fieldnames(Models);
nmods = size(mod_names,1);

%% Standard timecourses
if ~or(option.doTF,option.doTFbands)
    
epoch = option.anaepoch;
    
for mask = 1:length(option.masknic)
    
    % Load data
    masknam = [option.masknic{mask}];
    if option.tw > 1
    if mod(option.tw,2)  % must be even number
        option.tw = option.tw+1;
    end
    end
    
    load([option.rsafront masknam option.midname num2str(option.tw*option.srate) 'ms_sTW.mat'],'rsa_out')
    
    for modl = 1:nmods
        
        % Get data
        x = rsa_out.(mod_names{modl});
        x(isnan(x)) = 0; % set nans to zero
        x(option.sub_out,:) = [];
        % put data into shape
        for i = 1:length(x(:,1))
            r(i,1,:) = x(i,:);
        end
        a = zeros(size(r)); clear x i
        
        % Test
        [clust_stat_p, clust_mass_p, base_tmap] = rsa_permutation_cluster_test_2dtfr_func(r,a,1,1,epoch(1),epoch(2),option.perms,1,alpha,rand_perms);
        
        % work out start/stop of clusters
        start = nonzeros([1; diff(find(base_tmap>abs(tinv(alpha,(size(r,1)-1)))))>1] .* find(base_tmap>abs(tinv(alpha,(size(r,1)-1)))));
        stop = nonzeros([diff(find(base_tmap>abs(tinv(alpha,(size(r,1)-1)))))>1; 1] .* find(base_tmap>abs(tinv(alpha,(size(r,1)-1)))));
        if isempty(start)
            start = 0; stop = 0;
        else
            start = (start*option.srate)+1-((option.baseline-epoch(1)+1)*option.srate); % convert to ms
            stop = (stop*option.srate)-((option.baseline-epoch(1)+1)*option.srate); % convert to ms
        end
        
        % compile cluster size and permutation data
        r_clust = [];
        ii=0;
        if and(start,stop) == 0; r_clust = [];
        else
            for i = 1:length(clust_stat_p(:,1))
                if clust_stat_p(i,2) < 0.2  % only record if p < 0.2
                    ii=ii+1;
                    r_clust{ii,1} = mod_names{modl};  % model name
                    r_clust{ii,2} = masknam;       % ROI
                    r_clust{ii,3} = [num2str(start(i)) '-' num2str(stop(i)) ' ms'];       % times
                    r_clust{ii,4} = clust_stat_p(i,1);   % Cluster mass
                    r_clust{ii,5} = clust_stat_p(i,2);   % cluster p
                end
            end
        end
        clusters = [clusters; r_clust];        % cluster sizes with label index
        max_perm_hist = [max_perm_hist; clust_mass_p]; clear r i
        
    end % mod
end % roi

max_perm_hist2 = max(max_perm_hist);
for c = 1:length(clusters(:,1))
    max_p_val{c} = length(find(max_perm_hist2 > clusters{c,4}))/length(max_perm_hist2(1,:));
end

ROI_corrected_stats = [clusters, max_p_val'];
outname = [option.rfxdir 'RSA_timecourse_perm_stats_p' num2str(alpha) '.mat'];
save(outname,'ROI_corrected_stats','option','max_perm_hist');

end % standard

%% TFbands
if option.doTFbands        
    
    option.tw = round((option.tw*option.srate)/option.tfstep);
    if option.tw > 1
        if mod(option.tw,2)  % must be even number
            option.tw = option.tw+1;
        end
    end     
%     if mod((option.tw*option.srate)/option.tfstep,2)  % must be even number
%         tfw = ((option.tw*option.srate)/option.tfstep)+1;
%     else
%         tfw = ((option.tw*option.srate)/option.tfstep);
%     end

for mask = 1:length(option.masknic)
    
    masknam = [option.masknic{mask}];            
    
    % Load data
    masknam = [option.masknic{mask}];    
    if option.doPhase || option.doPhasePower% phase        
            load([option.tf_pre 'bands_phz_' option.rsafront masknam option.midname num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'rsa_out','option');        
    else        
            load([option.tf_pre 'bands_' option.rsafront masknam option.midname num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'rsa_out','option');        
    end
    
    
    % re-check some options
    options = optionsfile(option.s,1);        
    option = setfield(option,'sub_out',options.sub_out);
    option = setfield(option,'tails',options.tails);
    option = setfield(option,'alpha',options.alpha);
    option = setfield(option,'perms',options.perms);
    option = setfield(option,'stats_epoch',options.stats_epoch);
    clear options      
    
    freqs = option.fs;    
    epoch = [find(option.tfepoch == option.stats_epoch(1)) find(option.tfepoch == option.stats_epoch(2))];
    if isempty(epoch)  % use approximate
    epoch = [min(find((option.tfepoch>option.stats_epoch(1))))-1 ...
       max(find((option.tfepoch>option.stats_epoch(1))))];
    end
    
    for modl = 1:nmods
        r=[];
        % Get data
        x = rsa_out.(mod_names{modl});
        x(isnan(x)) = 0; % set nans to zero
        x(option.sub_out,:,:) = [];
        % put data into shape
        for i = 1:length(x(:,1,1))
            r(i,:,:) = x(i,:,:);
        end
        a = zeros(size(r)); clear x i
        
        for bands = 1:length(freqs(1,:)) 
        
        % Test
        [clust_stat_p, clust_mass_p, base_tmap] = rsa_permutation_cluster_test_2dtfr_func(r,a,bands,bands,epoch(1),epoch(2),option.perms,1,alpha,rand_perms);
        
        % work out start/stop of clusters
        start = nonzeros([1; diff(find(base_tmap>abs(tinv(alpha,(size(r,1)-1)))))>1] .* find(base_tmap>abs(tinv(alpha,(size(r,1)-1)))));
        stop = nonzeros([diff(find(base_tmap>abs(tinv(alpha,(size(r,1)-1)))))>1; 1] .* find(base_tmap>abs(tinv(alpha,(size(r,1)-1)))));
        if isempty(start)
            start = 0; stop = 0;
        else
            start = option.tfepoch(1,start+epoch(1)-1);
            stop = option.tfepoch(1,stop+epoch(1)-1);
        end
        
        % compile cluster size and permutation data
        r_clust = [];
        ii=0;
        if and(start,stop) == 0; r_clust = [];
        else
            for i = 1:length(clust_stat_p(:,1))
                if clust_stat_p(i,2) < 0.2  % only record if p < 0.2
                    ii=ii+1;
                    r_clust{ii,1} = mod_names{modl};  % model name
                    r_clust{ii,2} = masknam;       % ROI
                    r_clust{ii,3} = [num2str(freqs{bands}(1)) '-' num2str(freqs{bands}(end)) ' Hz'];       % freq band
                    r_clust{ii,4} = [num2str(start(i)) '-' num2str(stop(i)) ' ms'];       % times
                    r_clust{ii,5} = clust_stat_p(i,1);   % Cluster mass
                    r_clust{ii,6} = clust_stat_p(i,2);   % cluster p
                end
            end
        end
        clusters = [clusters; r_clust];        % cluster sizes with label index
        max_perm_hist = [max_perm_hist; clust_mass_p]; clear i
        
        end % bands        
    end % mod
end % roi

if ~isempty(clusters)
    max_perm_hist2 = max(max_perm_hist);
    for c = 1:length(clusters(:,1))
        max_p_val{c} = length(find(max_perm_hist2 > clusters{c,5}))/length(max_perm_hist2(1,:));
    end
    
    ROI_corrected_stats = [clusters, max_p_val'];
else
    ROI_corrected_stats = [clusters];
end

option.alpha = alpha;
outname = [option.rfxdir 'TFbands_RSA_timecourse_perm_stats_p' num2str(alpha) '.mat'];
save(outname,'ROI_corrected_stats','option','max_perm_hist','freqs','mod_names');
end % tfbands


if option.doTF && option.doTFbands==0
    option.tw = round((option.tw*option.srate)/option.tfstep);
    if option.tw > 1
        if mod(option.tw,2)  % must be even number
            option.tw = option.tw+1;
        end
    end        
%     if mod((option.tw*option.srate)/option.tfstep,2)  % must be even number
%         tfw = ((option.tw*option.srate)/option.tfstep)+1;
%     else
%         tfw = ((option.tw*option.srate)/option.tfstep);
%     end

for mask = 1:length(option.masknic)
    
    masknam = [option.masknic{mask}];            
    
    % Load data
    masknam = [option.masknic{mask}];
    
    if option.doPhase || option.doPhasePower % phase        
            load([option.tf_pre 'phz_' option.rsafront masknam option.midname num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'rsa_out','option');       
    else % Power                
            load([option.tf_pre option.rsafront masknam option.midname num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'rsa_out','option');        
    end
    
    % re-check some options
    options = optionsfile(option.s,1);        
    option = setfield(option,'sub_out',options.sub_out);
    option = setfield(option,'tails',options.tails);
    option = setfield(option,'alpha',options.alpha);
    option = setfield(option,'perms',options.perms);
    option = setfield(option,'stats_epoch',options.stats_epoch);
    clear options          
    
    freqs = option.fsc;
    epoch = [find(option.tfepoch == option.stats_epoch(1)) find(option.tfepoch == option.stats_epoch(2))];
    if isempty(epoch)  % use approximate
    epoch = [min(find((option.tfepoch>option.stats_epoch(1))))-1 ...
       max(find((option.tfepoch>option.stats_epoch(1))))];
    end
    times = option.tfepoch(epoch(1):epoch(2));
    
    for modl = 1:nmods
        r=[];
        % Get data
        r = rsa_out.(mod_names{modl});
        r(isnan(r)) = 0; % set nans to zero
        r(option.sub_out,:,:) = [];        
        a = zeros(size(r));                
        
        % Test
        [clust_stat_p, clust_mass_p, base_tmap] = rsa_permutation_cluster_test_2dtfr_func(r,a,1,size(r,2),epoch(1),epoch(2),option.perms,1,alpha,rand_perms);        
        clust_stat_p(find(clust_stat_p(:,1) == 0),:) = []; % get rid of clusters that are sized zero
        
        % Get time and freq bounds of clusters
%        clusters_matrix_p = bwlabel(base_tmap>abs(tinv(alpha,(size(r,1)-1))));
%        clusters_matrix_n = bwlabel(base_tmap<(tinv(alpha,(size(r,1)-1))));
%        % Get combined clusters matrix
%        clusters_matrix_n(find(clusters_matrix_n)) = clusters_matrix_n(find(clusters_matrix_n))+max(max(clusters_matrix_p));  % adjust neg labels
%        clusters_matrix = clusters_matrix_p + clusters_matrix_n;
        
        clusters_matrix = bwlabel(base_tmap>abs(tinv(alpha,(size(r,1)-1))));
        
        if sum(clusters_matrix(:)) == 0; %isempty(clusters_matrix)
            Fstart = 0;
            Fstop = 0;
            Tstart = 0;
            Tstop = 0;
        else
            for i = 1:max(max(clusters_matrix))
                [F,T] = ind2sub([size(r,2),size(r,3)],find(clusters_matrix==i));
                Fstart(i) = round(freqs(min(F)));
                Fstop(i) = round(freqs(max(F)));
                Tstart(i) = times(min(T));
                Tstop(i) = times(max(T));
            end
        end
        
        r_clust = [];
        ii=0;
        % compile cluster size and permutation data
        if min(Fstart == 0)
            r_clust = [];
        else
            for i = 1:length(clust_stat_p(:,1))
                if clust_stat_p(i,2) < 0.5 % don't both with smaller clusters
                    ii=ii+1;
                    r_clust{ii,1} = mod_names{modl};  % model name
                    r_clust{ii,2} = masknam;       % ROI                    
                    r_clust{ii,3} = [num2str(Fstart(i)) '-' num2str(Fstop(i)) ' Hz'];       % freqs
                    r_clust{ii,4} = [num2str(Tstart(i)) '-' num2str(Tstop(i)) ' ms'];       % times
                    r_clust{ii,5} = clust_stat_p(i,1);   % Cluster mass
                    r_clust{ii,6} = clust_stat_p(i,2);   % cluster p
                end
            end
        end        
        clusters = [clusters; r_clust];        % cluster sizes with label index
        max_perm_hist = [max_perm_hist; clust_mass_p];
        clear Fstart Fstop Tstart Tstop i                
    end % mod
end % roi

if ~isempty(clusters)
    max_perm_hist2 = max(max_perm_hist);
    for c = 1:length(clusters(:,1))
        max_p_val{c} = length(find(max_perm_hist2 > clusters{c,5}))/length(max_perm_hist2(1,:));
    end    
    ROI_corrected_stats = [clusters, max_p_val'];
else
    ROI_corrected_stats = [clusters];
end

outname = [option.rfxdir 'TF_RSA_timecourse_perm_stats_p' num2str(alpha) '.mat'];
save(outname,'ROI_corrected_stats','option','max_perm_hist','freqs','mod_names','times');
end % tf
