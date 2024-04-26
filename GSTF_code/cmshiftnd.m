function y=cmshiftnd(x,shifts)

% Function to circularly shift N-D arrays
% M.A. Griswold 9/11/98
%
% FUNCTION
%   y = cmshiftnd(x,shifts)
%
% INPUT
%   x           Array                       (N-dims)
%   shifts      Vector containing shifs     (1,N)
%                   one shift for each dimension
%                   positive entry leads to a shift to negative indexes
%
% OUTPUT
%   y           Shifted array
%



if nargin < 2
   shifts=0;						% no shift
end

numDims = ndims(x);				% number of dimensions
idx = cell(1, numDims);			% creates cell array of empty matrices, 
										% one cell for each dimension
for k = 1:numDims
    m = size(x, k);
    p = ceil(shifts(k));

	if p < 0
		p=m+p;
	end

    idx{k} = [p+1:m 1:p];
end

% Use comma-separated list syntax for N-D indexing.
y = x(idx{:});
