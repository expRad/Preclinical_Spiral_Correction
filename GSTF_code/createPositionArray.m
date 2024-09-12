function [ positions ] = createPositionArray( orientation, numPE, numSlices, FOV, PosSlices )
%
% Copyright (c) 2024 Hannah Scholten
%%
positions = zeros(numPE, numPE, numSlices, 3); % [PE steps, partitions, slices, 3], last dimension is [x y z]

delta = FOV/numPE;

%% Transversal slices:
% slice = z, PE = y, RO = x
if strcmp(orientation,'dTra')
    for phase=1:1:numPE
        y = (-FOV + delta)/2 + (phase-1)*delta;
        for part=1:1:numPE
            x = (-FOV + delta)/2 + (part-1)*delta;
            for slice=1:1:numSlices
                z = PosSlices(slice);
                positions(phase,part,slice,1) = x;
                positions(phase,part,slice,2) = y;
                positions(phase,part,slice,3) = z;
            end
        end
    end
    disp('Position array for transversal slices successfully created.')
end

%% Saggital slices
% slice = x, PE = y, RO = z
if strcmp(orientation,'dSag')
    for phase=1:1:numPE
        y = (-FOV + delta)/2 + (phase-1)*delta;
        for part=1:1:numPE
            z = (-FOV + delta)/2 + (part-1)*delta;
            for slice=1:1:numSlices
                x = PosSlices(slice);
                positions(phase,part,slice,1) = x;
                positions(phase,part,slice,2) = y;
                positions(phase,part,slice,3) = z;
            end
        end
    end
    disp('Position array for sagittal slices successfully created.')
end

%% Coronal slices
% slice = y, PE = x, RO = z
if strcmp(orientation,'dCor')
    for phase=1:1:numPE
        x = (-FOV + delta)/2 + (phase-1)*delta;
        for part=1:1:numPE
            z = (-FOV + delta)/2 + (part-1)*delta;
            for slice=1:1:numSlices
                y = PosSlices(slice);
                positions(phase,part,slice,1) = x;
                positions(phase,part,slice,2) = y;
                positions(phase,part,slice,3) = z;
            end
        end
    end
    disp('Position array for coronal slices successfully created.')
end

positions = positions / 1000; % convert from mm to m
end