function [ out_signals, magnitudeTooLow ] = calcMeasOutput_Bruker_avgComplex( girf_data, numTriang, calcChannels, t_ADC)
%CALCMEASOUTPUT calculates the measured gradient output 
%%
kspace = girf_data.raw; % [ROP, slices, rep, triangles]
numROP = size(kspace,1); % number of Read Out Points
numSlices = size(kspace,2);
numRep = size(kspace,3); % number of repetitions (measurements)
numPulses = size(kspace,4);
dwelltime = 1/double(girf_data.bandwidth);
PosSlices = [-5.1; 5.1].*1e-3;
% kspace: [ROPs, slices, averages, 1+numTriangles], first measurements are reference
orientation = 'dTra'; % dTra, dSag or dCor, depending on the slice orientation
FOV = 200;
numPE = 1;

%% Average over repetitions
kspace = mean(kspace,3);
numRep = 1;
disp(['size(kspace)=',num2str(size(kspace))])

%% Get magnitude and phase data
diff_phase = zeros(size(kspace));
magnitude = abs(kspace);
for slice = 1:size(kspace,2)
    for triang = 1:size(kspace,4)
        phase_pp = interp1(t_ADC, unwrap(angle(squeeze(kspace(:,slice,:,triang))),[],1),'linear','pp');
        diff_phase_pp = fnder(phase_pp);
        diff_phase(:,slice,:,triang) = ppval(diff_phase_pp, t_ADC);
    end
end
clearvars kspace;

%% Separate data with test gradients from reference data without test gradients
diff_phase_plus = diff_phase(:,:,:,2:end); % [numROP, numSlices, numRep(1), numTriang]
diff_phase_minus = diff_phase(:,:,:,1);

%% Calculate difference between triangles and reference
delta_diff_phase = (diff_phase_plus - diff_phase_minus); % [numROP, numSlices, numRep(1), numTriang]
clearvars diff_phase_plus diff_phase_minus;

delta_diff_phase = permute(delta_diff_phase, [2,3,1,4]); % [numSlices, 1, numROP, numTriang]
delta_diff_phase = reshape(delta_diff_phase, [numSlices, numROP*numTriang]); % [numVoxels, numTimePoints]

%% Get the positions of the measured voxels
positions = createPositionArray(orientation, numPE, numSlices, FOV, PosSlices); % [numPE(PE), numPE(part), numSlices, 3]
positions = reshape(positions, [numPE*numPE*numSlices, 3])*1000; % [numVoxels, 3]

%% Sort out unusable voxels
% by thresholding the signal magnitude
mag = permute(magnitude, [2,3,1,4]); % [numSlices, 1, numROP, numTriang]
mag = reshape(mag, [numSlices, numROP*(numTriang+1)]); % [numVoxels, numTimePoints]
mag = squeeze(mean(mag(:,10:50),2));
max_mag = max(mag);

validVoxels = zeros(size(positions,1),1) + 1;

for voxel=numPE*numPE*numSlices:-1:1
    r = sqrt(positions(voxel,1)*positions(voxel,1) + positions(voxel,2)*positions(voxel,2) + positions(voxel,3)*positions(voxel,3));
    if 0
        positions(voxel,:) = [];
        delta_diff_phase_interp(voxel,:) = [];
        validVoxels(voxel) = 0;
    elseif mag(voxel) < max_mag*0.6
        positions(voxel,:) = [];
        delta_diff_phase_interp(voxel,:) = [];
        validVoxels(voxel) = 0;
    end
end

for voxel=1:1:size(positions,1)
    r = sqrt(positions(voxel,1)*positions(voxel,1) + positions(voxel,2)*positions(voxel,2) + positions(voxel,3)*positions(voxel,3));
    disp(['voxel ',num2str(voxel)])
    disp(['x=',num2str(positions(voxel,1)),', y=',num2str(positions(voxel,2)),', z=',num2str(positions(voxel,3)),'\n'])
    disp(['r = ',num2str(r)])
end

FIDs.validVoxels = validVoxels;

%% Get the probing matrix
if size(positions,1)<4
    probingMatrix = zeros(size(positions,1), calcChannels);
    for slice=1:1:size(positions,1)
        probingMatrix(slice,1) = 1;
        if calcChannels>1
            if strcmp(orientation,'dTra')
                slicePosition = positions(slice,3);
            elseif strcmp(orientation,'dSag')
                slicePosition = positions(slice,1);
            elseif strcmp(orientation,'dCor')
                slicePosition = positions(slice,2);
            end
            probingMatrix(slice,2) = slicePosition;
            if calcChannels>2
                probingMatrix(slice,3) = 2*slicePosition*slicePosition;
            end
        end
    end
else
    probingMatrix = createProbingMatrix(positions); % [numValidVoxels, 4/9/16], depending on the maximum expansion order
end

%% Calculate the output signals
gamma = 267.513*10^6; %Hz/T
out_signals = 1/gamma * (probingMatrix\delta_diff_phase); % [channels, numTimePoints]
out_signals = reshape(out_signals, [calcChannels, numROP, numTriang]); % [calcChannels, numROP, numTriang]

magnitudeTooLow = zeros(size(out_signals)) + 1;

end

