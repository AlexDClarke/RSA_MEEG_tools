function [clust_stats, clust_mass_perm, base_tmap] = rsa_permutation_cluster_test_2dtfr_func(data1,data2,fpmin,fpmax,tpmin,tpmax,perm_num,test,pval,rand_perms)

%%% Cluster based permutation t-test between conditions of 2D TF data
%
% Inputs:
%       data1 = condition 1 data (subjects x freqs x timepoint)
%       data2 = condition 2
%       fpmin = lowest frequency of interest (e.g. 5)
%       fpmax = highest frequency of interest (e.g. 60)
%       tpmin = Earliest timepoint of interest
%       tpmax = latest timepoint of interest
%       perm_num = number of permutation (e.g. 5000)
%       test = one or paired sample t-test (1 = one, 0 = paired [default])
%       pval = one-tailed sig threshold wanted (e.g. 0.025 [default])
%       rand_perms = matrix of permutations to apply to the data
%
% Output:
%       p_val = p value of largest cluster size according to the
%               permutation distribution.
%
% % Note the timepoints will relate to time-sample number, not necessarily
% % time. Min and max frequency relate to the row number, so if data begins
% % at 5Hz, the fpmin = 2 means a freq of 6Hz.
%
% Writen bu A.Ghuman (2009)
% Adapted by Alex C (07/09) to calculate p value for all clusters, and
% return permutation distribution.

if nargin == 8
    pval = 0.025 %
else if nargin == 7
        test = 0; % paired samples t-test
        pval = 0.025 %
    end
end

% Extract data specified by options
data_1_rel = data1(:,fpmin:fpmax,tpmin:tpmax);
data_2_rel = data2(:,fpmin:fpmax,tpmin:tpmax);

% Calculate t-statistic map
base_tmap = squeeze(mean(data_1_rel-data_2_rel)./(std(data_1_rel-data_2_rel)/sqrt(size(data_1_rel,1))));

% Find cluser sizes for sig clusters (positive only)
[~, clust_pos_data] = find_max_clust_mass(base_tmap,abs(tinv(pval,(size(data_1_rel,1)-1))));

% Permute data
if test
    for i = 1:perm_num
% For one-sample test (flip signs not re-label conditions)
        perm_rand = rand_perms(:,i);
        data_1_temp = [data_1_rel(find(perm_rand==1),:,:);-1*(data_1_rel(find(perm_rand==0),:,:))];
        perm_tmap = squeeze(mean(data_1_temp)./(std(data_1_temp)/sqrt(size(data_1_rel,1))));

        clust_mass_pos = find_max_clust_mass(perm_tmap,abs(tinv(pval,(size(data_1_rel,1)-1))));
        clust_mass_neg=0;
        clust_mass_perm(i) = max(clust_mass_neg,clust_mass_pos);
    end
else
    for i = 1:perm_num
% For paired-test (re-labeling)
        perm_rand = rand_perms(:,i);
%       perm_rand = round(rand(size(data_1_rel,1),1));
        data_1_temp = [data_1_rel(find(perm_rand==1),:,:);data_2_rel(find(perm_rand==0),:,:)];
        data_2_temp = [data_2_rel(find(perm_rand==1),:,:);data_1_rel(find(perm_rand==0),:,:)];
        perm_tmap=squeeze(mean(data_1_temp-data_2_temp)./(std(data_1_temp-data_2_temp)/sqrt(size(data_1_rel,1))));

        clust_mass_pos = find_max_clust_mass(perm_tmap,abs(tinv(pval,(size(data_1_rel,1)-1))));       
        clust_mass_neg=0;
        clust_mass_perm(i) = max(clust_mass_neg,clust_mass_pos);
    end
end % if test

% p-value of each cluster
for c = 1:length(clust_pos_data)
%   p_val(c) = length(find(clust_mass_perm >= clust_pos_data(c)))/perm_num;
    p_val(c) = length(find([clust_mass_perm max(clust_pos_data)] >= clust_pos_data(c)))/(perm_num+1);
end

% cluster stats output
% % Outputs cluster size as well as cluster mass although this is only for
% % the purpose of identifying which p-val refers to which cluster

clust_stats = [clust_pos_data' p_val'];
