function roiRSA_permstats_3D(option)

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
        x(isnan(x)) = 0;
        x(option.sub_out,:) = [];
        % put data into shape
        for i = 1:length(x(:,1))
            r(i,modl,:) = x(i,:);
        end
    end    
    a = zeros(size(r)); clear x i
    times = [round(option.start):option.srate:round(option.stop)];
    times = times(epoch(1):epoch(2));
    
    % Test
    [clust_stat_p, clust_mass_p, base_tmap] = rsa_permutation_cluster_test_2dtfr_func(r,a,1,size(r,2),epoch(1),epoch(2),option.perms,1,alpha,rand_perms);
    clust_stat_p(find(clust_stat_p(:,1) == 0),:) = []; % get rid of clusters that are sized zero
               
    clusters_matrix = bwlabel(base_tmap>abs(tinv(alpha,(size(r,1)-1))));
    if sum(clusters_matrix(:)) == 0; %isempty(clusters_matrix)
        Fstart = 0;
        Fstop = 0;
        Tstart = 0;
        Tstop = 0;
    else
        for i = 1:max(max(clusters_matrix))
            [F,T] = ind2sub([size(r,2),size(r,3)],find(clusters_matrix==i));
            Fstart(i) = min(F);
            Fstop(i) = max(F);
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
                r_clust{ii,1} = option.models;  % model set name
                r_clust{ii,2} = masknam;       % ROI
                r_clust{ii,3} = [num2str(Fstart(i)) '-' num2str(Fstop(i)) ' ticks'];       % freqs
                r_clust{ii,4} = [num2str(Tstart(i)) '-' num2str(Tstop(i)) ' ms'];       % times
                r_clust{ii,5} = clust_stat_p(i,1);   % Cluster mass
                r_clust{ii,6} = clust_stat_p(i,2);   % cluster p
            end
        end
    end
        clusters = [clusters; r_clust];        % cluster sizes with label index
        max_perm_hist = [max_perm_hist; clust_mass_p];
        clear Fstart Fstop Tstart Tstop i              
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

outname = [option.rfxdir 'RSA_3D_timecourse_perm_stats_p' num2str(alpha) '.mat'];
save(outname,'ROI_corrected_stats','option','max_perm_hist','mod_names','times');

end % tf
