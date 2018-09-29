function nn= hist4( x, nbins )
%HIST4 Summary of this function goes here
%   Detailed explanation goes here
if nargin==1
    nbins=[10,10,10];
elseif numel(nbins)~=3
    error(message('stats:hist4:WrongNumBins'));
end

[nrows,ncols] = size(x);
if ncols ~= 3
    error(message('stats:hist4:WrongNumCols'));
end

% Special case for empty data (follows what HIST does).
if isempty(x)

    n = zeros(nbins); % Nothing to count, return nbins(1) by nbins(2) zeros
    
else
    % Bin each observation in the x-direction, and in the y-direction.
    bin = zeros(nrows,3);
    for i = 1:3
        minx = min(x(:,i));
        maxx = max(x(:,i));
        
        % If only the number of bins was given, compute edges and centers
        % for equal-sized bins spanning the data.

        if isinf(minx) || isinf(maxx)
            error(message('stats:hist3:InfData'));
        elseif minx == maxx
            minx = minx - floor(nbins(i)/2) - 0.5;
            maxx = maxx + ceil(nbins(i)/2) - 0.5;
        end
        binwidth{i} = (maxx - minx) / nbins(i);
        edges{i} = minx + binwidth{i}*(0:nbins(i));
        ctrs{i} = edges{i}(1:nbins(i)) + binwidth{i}/2;
        % Make histc mimic hist behavior:  everything < ctrs(1) gets
        % counted in first bin, everything > ctrs(end) gets counted in
        % last bin.  ctrs, edges, and binwidth do not reflect that, but
        % histcEdges does.
        histcEdges = [-Inf edges{i}(2:end-1) Inf];

        % Get the 1D bin numbers for this column of x.  Make sure +Inf
        % goes into the nth bin, not the (n+1)th.
        [dum,bin(:,i)] = histc(x(:,i),histcEdges,1);
        bin(:,i) = min(bin(:,i),nbins(i));
    end
    
    % Combine the two vectors of 1D bin counts into a grid of 2D bin
    % counts.
    %n = accumarray(bin(all(bin>0,2),:),1,nbins);  %%for test
    n=accumarray(bin,1,nbins);
end

if 0 < nargout
    nn = n;
    return
end

end

