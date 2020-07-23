% taua = rankCorr_Kendall_taua(a,b)
%
% computes the Kendall's tau a correlation coefficient between the input
% vectors (a and b). NaN entries would be removed from both.
%__________________________________________________________________________
% Copyright (C) 2012 Medical Research Council
function taua = rankCorr_Kendall_taua(a,b)

    %% preparations
    a = a(:);
    b = b(:);
    validEntryIs = ~isnan(a) & ~isnan(b); % Use bitwise & here, not &&
    a = a(validEntryIs);
    b = b(validEntryIs);
    n = numel(a);

    %% compute Kendall rank correlation coefficient tau_a
    K = 0;
    for i = 1 : n-1
        pairRelations_a = sign(a(i)-a(i+1:n));
        pairRelations_b = sign(b(i)-b(i+1:n));
        K = K + sum(pairRelations_a .* pairRelations_b);
    end
    taua = K/(n*(n-1)/2); % normalise by the total number of pairs 

end%function
