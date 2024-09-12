function [ probingMatrix ] = createProbingMatrix( positions )
%
% Copyright (c) 2024 Hannah Scholten
%%
numVoxels = size(positions,1); % positions has size [numVoxels, 3]

% The number of voxels determines up to which order we can expand the
% measured field
if numVoxels > 15
    probingMatrix = zeros(numVoxels, 16); % expansion up to third order possible
elseif numVoxels > 8
    probingMatrix = zeros(numVoxels, 9); % expansion up to second order possible
elseif numVoxels > 3
    probingMatrix = zeros(numVoxels, 4); % expansion only to first order
else
    probingMatrix = 0;
    disp('ERROR: Not enough voxels!')
end

if numVoxels > 3
    for p=1:1:numVoxels
        x = positions(p,1);
        y = positions(p,2);
        z = positions(p,3);
        
        % The following formulas for the real-valued sperical harmonics were taken
        % from https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics
        probingMatrix(p,1) = 1;
        
        probingMatrix(p,2) = x;
        probingMatrix(p,3) = y;
        probingMatrix(p,4) = z;
        
        if numVoxels > 8
            probingMatrix(p,5) = x*y;
            probingMatrix(p,6) = y*z;
            probingMatrix(p,7) = (2*z*z-x*x-y*y);
            probingMatrix(p,8) = z*x;
            probingMatrix(p,9) = (x*x-y*y);
            
            if numVoxels > 15
                probingMatrix(p,10) = (3*x*x-y*y)*y;
                probingMatrix(p,11) = x*y*z;
                probingMatrix(p,12) = (4*z*z-x*x-y*y)*y;
                probingMatrix(p,13) = (2*z*z-3*x*x-3*y*y)*z;
                probingMatrix(p,14) = (4*z*z-x*x-y*y)*x;
                probingMatrix(p,15) = (x*x-y*y)*z;
                probingMatrix(p,16) = (x*x-3*y*y)*x;
            end
        end
    end
end

disp(['probingMatrix successfully created. size: ',num2str(size(probingMatrix))])

end










