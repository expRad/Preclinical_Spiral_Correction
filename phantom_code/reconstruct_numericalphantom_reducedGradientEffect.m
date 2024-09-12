%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2024 Hannah Scholten
%
% This script demonstrates that the lowpass characteristic of the gradient
% system leads to a smaller appearance of the object if it is not
% compensated for.
%
% numerical phantom = circle
% abs(GSTF) = constant = 0.9
% angle(GSTF) = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear all;
% Change current directory to that of this .m file
mfile_name          = mfilename('fullpath');
if contains(mfile_name,'LiveEditorEvaluationHelper')
    mfile_name = matlab.desktop.editor.getActiveFilename;
end
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);
% Add necessary scripts to the MATLAB path
folder1 = "..\GSTF_code";
folder2 = "..\Spiral_Recon_NUFFT";
addpath(folder1, folder2);
% Run setup for MIRT toolbox
run ('..\MIRT_toolbox\irt\setup.m');
traj_factor = 6;
matrix_size = 192;

%% initialize data structures
meas1 = struct;
meas1.filepath = '..\phantom_data\spiral_data_p96.mat';

meas2 = struct;
meas2.filepath = '..\phantom_data\spiral_data_p16.mat';

meas3 = struct;
meas3.filepath = '..\phantom_data\spiral_data_p03.mat';

measurements = {meas1, meas2, meas3};

%% initialize arrays for results
recos_nom = zeros(matrix_size, matrix_size, length(measurements));
recos_girf = zeros(matrix_size, matrix_size, length(measurements));

%% create numerical phantom
P = zeros(matrix_size);
[xx,yy] = meshgrid((1:matrix_size),(1:matrix_size));
P = double(sqrt((xx-matrix_size/2).^2+(yy-matrix_size/2).^2)<70);

%% load GSTFs
H_x = load('..\GSTF_data\Hx_sr100_20231009.mat');
H_y = load('..\GSTF_data\Hy_sr100_20231009.mat');
H_z = load('..\GSTF_data\Hz_sr100_20231009.mat');

%% Run simulations for every spiral trajectory
for loop = 1:length(measurements)
    disp(['loop = ',num2str(loop)]);
    meas = measurements{loop};
    filepath = meas.filepath;
    
    % Initialize transfer functions
    H_1 = H_x.H_matrix;
    gstf_x = H_x.H_matrix.gstf(:,2);
    
    gstf_x = ones(size(gstf_x))*0.9;
    gstf_y = gstf_x;
    gstf_z = gstf_x;

    dwelltime_girf = H_x.H_matrix.dt;
    
    %% Load spiral data
    disp('Load data.')

    if isfile(filepath)
        load(filepath);
    else
        disp(['Error: cannot find ',filepath]);
        break
    end
    nPost = spiral_data.nPost; % data points discarded for reconstruction
    measurements{loop}.nPost = nPost;

    disp(['NUMBER of SPIRAL INTERLEAVES: ',num2str(spiral_data.numInterleaves)]);

    %% get nominal trajectory
    % calculate in-plane gradients for each spiral interleaf
    spiral_x = zeros(spiral_data.numPointsPerInterleave,spiral_data.numInterleaves);
    spiral_y = zeros(size(spiral_x));
    max_grad_strength = 1.5212; % T/m
    grad_raster_time = 8.1e-6;
    dwelltime_meas = 1/double(spiral_data.bandwidth);

    for arm = 1:spiral_data.numInterleaves
        spiral_x(:,arm) = spiral_data.spiralShape1*spiral_data.spiralInterleaveCos(arm) - spiral_data.spiralShape2*spiral_data.spiralInterleaveSin(arm);
        spiral_y(:,arm) = spiral_data.spiralShape1*spiral_data.spiralInterleaveSin(arm) + spiral_data.spiralShape2*spiral_data.spiralInterleaveCos(arm);
        % spiralshapes are normalized to 1
    end

    % interpolate to ADC time grid
    t_ADC = (0:1:(double(spiral_data.datapoints)-1))*dwelltime_meas + dwelltime_meas/2;
    measurements{loop}.t_ADC = t_ADC;
    t_GRT = ((1:1:(size(spiral_x,1)))-0.5)*grad_raster_time; % + PVM_SpiralPreSize, aber ist bei uns immer(?) 0
    t_eps = 1e-12;
    t_GRT_00 = [-t_eps, t_GRT(1)-t_eps, t_GRT, t_GRT(end)+t_eps, t_GRT(end)+grad_raster_time+t_eps];
    
    spiral_x_00 = [zeros(2,size(spiral_x,2)); spiral_x; zeros(2,size(spiral_x,2))];
    spiral_y_00 = [zeros(2,size(spiral_x,2)); spiral_y; zeros(2,size(spiral_x,2))];

    spiral_x_pp = interp1(t_GRT_00, spiral_x_00*max_grad_strength, 'linear', 'pp');
    spiral_y_pp = interp1(t_GRT_00, spiral_y_00*max_grad_strength, 'linear', 'pp');

    % integrate to get k-space trajectory
    traj_x_pp = fnint(spiral_x_pp);
    traj_y_pp = fnint(spiral_y_pp);

    %% calculate nominal trajectory
    gamma = 267.513*10^6; %Hz/T
    trajx_nom = gamma/(2*pi)*ppval(traj_x_pp, t_ADC)/1000;
    trajy_nom = gamma/(2*pi)*ppval(traj_y_pp, t_ADC)/1000;
    traj_nom = trajx_nom.' + 1i*trajy_nom.';
    measurements{loop}.traj_nom = traj_nom;

    %% interpolate to gstf time grid (so the frequency range matches)
    t_GIRF = H_1.t_axis;
    if t_GIRF(end) < t_GRT(end)
        t_GIRF_end = ceil(t_GRT(end)/dwelltime_girf);
        t_GIRF = ((0:t_GIRF_end)+0.5)*dwelltime_girf;
    end
    grad_x_tgirf = ppval(spiral_x_pp, t_GIRF);
    grad_y_tgirf = ppval(spiral_y_pp, t_GIRF);

    %% rotate in-plane gradients to XYZ coordinate system
    % and calculate concomitant field for the position x0 = y0 = z0 = 1cm
    grad_physicalX = zeros(size(grad_x_tgirf));
    grad_physicalY = zeros(size(grad_x_tgirf));
    grad_physicalZ = zeros(size(grad_x_tgirf));

    for arm = 1:spiral_data.numInterleaves
        grad_xyz = [squeeze(grad_x_tgirf(:,arm)), squeeze(grad_y_tgirf(:,arm)), squeeze(zeros(size(t_GIRF.')))].';
        gradXYZ = spiral_data.gradMatrix.' * grad_xyz;
        grad_physicalX(:,arm) = gradXYZ(1,:);
        grad_physicalY(:,arm) = gradXYZ(2,:);
        grad_physicalZ(:,arm) = gradXYZ(3,:);
    end
    
    grad_physicalX_pp = interp1(t_GIRF, grad_physicalX, 'linear','pp');
    grad_physicalY_pp = interp1(t_GIRF, grad_physicalY, 'linear','pp');
    measurements{loop}.grad_x_nom = ppval(grad_physicalX_pp, t_ADC);
    measurements{loop}.grad_y_nom = ppval(grad_physicalY_pp, t_ADC);

    %% zerofill gradients to avoid side effects by the fft calculation
    nExtra = round((1e6-size(grad_x_tgirf,1))/2);
    grad_x_ZF = cat(1, zeros(nExtra,size(spiral_x,2)), grad_physicalX, zeros(nExtra,size(spiral_x,2)));
    grad_y_ZF = cat(1, zeros(nExtra,size(spiral_y,2)), grad_physicalY, zeros(nExtra,size(spiral_x,2)));
    grad_z_ZF = cat(1, zeros(nExtra,size(spiral_x,2)), grad_physicalZ, zeros(nExtra,size(spiral_x,2)));

    gstf_x_interp = interp1(H_1.f_axis, gstf_x, linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).', 'linear');
    gstf_y_interp = interp1(H_1.f_axis, gstf_y, linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).', 'linear');
    gstf_z_interp = interp1(H_1.f_axis, gstf_z, linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).', 'linear');

    %% calculate post-correction with GSTF
    tmp1 = fft1d(grad_x_ZF,1);
    tmp2 = repmat(gstf_x_interp, [1,spiral_data.numInterleaves]);
    grad_px_girf = real(ifft1d( tmp1.*tmp2, 1));
    grad_px_girf = grad_px_girf(nExtra+1:end-nExtra,:);

    tmp1 = fft1d(grad_y_ZF,1);
    tmp2 = repmat(gstf_y_interp, [1,spiral_data.numInterleaves]);
    grad_py_girf = real(ifft1d( tmp1.*tmp2, 1));
    grad_py_girf = grad_py_girf(nExtra+1:end-nExtra,:);

    tmp1 = fft1d(grad_z_ZF,1);
    tmp2 = repmat(gstf_z_interp, [1,spiral_data.numInterleaves]);
    grad_pz_girf = real(ifft1d( tmp1.*tmp2, 1));
    grad_pz_girf = grad_pz_girf(nExtra+1:end-nExtra,:);
    
    grad_px_girf_pp = interp1(t_GIRF, grad_px_girf, 'linear','pp');
    grad_py_girf_pp = interp1(t_GIRF, grad_py_girf, 'linear','pp');
    measurements{loop}.grad_x_girf = ppval(grad_px_girf_pp, t_ADC);
    measurements{loop}.grad_y_girf = ppval(grad_py_girf_pp, t_ADC);

    % rotate back to PRS (xyz/ in-plane) coordinate system
    grad_x_girf = zeros(size(grad_x_tgirf));
    grad_y_girf = zeros(size(grad_x_tgirf));

    for arm = 1:spiral_data.numInterleaves
        gradXYZ = [squeeze(grad_px_girf(:,arm)), squeeze(grad_py_girf(:,arm)), squeeze(grad_pz_girf(:,arm))].';
        grad_xyz = spiral_data.gradMatrix * gradXYZ;
        grad_x_girf(:,arm) = grad_xyz(1,:);
        grad_y_girf(:,arm) = grad_xyz(2,:);
    end

    % interpolate back to ADC time grid for reconstruction
    grad_x_girf_pp = interp1(t_GIRF, grad_x_girf, 'linear','pp');
    grad_y_girf_pp = interp1(t_GIRF, grad_y_girf, 'linear','pp');
    traj_x_girf_pp = fnint(grad_x_girf_pp);
    traj_y_girf_pp = fnint(grad_y_girf_pp);

    clearvars grad_x_girf grad_y_girf gradXYZ grad_px_girf grad_py_girf grad_pz_girf tmp1 tmp2 

    %% calculate GSTF-"corrected" trajectory
    trajx_corr = gamma/(2*pi)*ppval(traj_x_girf_pp, t_ADC)/1000;
    trajy_corr = gamma/(2*pi)*ppval(traj_y_girf_pp, t_ADC)/1000;
    traj_corr = trajx_corr.' + 1i*trajy_corr.';
    measurements{loop}.traj_girf = traj_corr;

    %% create synthetic raw data
    disp('    Create raw data with GSTF-"corrected" trajectory')
    FT = cGrid(traj_corr/traj_factor,matrix_size);
    raw = FT*P;

    %% reconstruct with GSTF-corrected trajectory
    disp('    Reconstruction with GSTF-corrected trajectory')
    FT = cGrid(traj_corr(1:end-nPost,:)/traj_factor,matrix_size);
    reco_girf = FT'*raw(1:end-nPost,:);

    %% reconstruct with nominal trajectory
    disp('    Reconstruction with nominal trajectory')
    FT = cGrid(traj_nom(1:end-nPost,:)/traj_factor,matrix_size);
    reco_nom = FT'*raw(1:end-nPost,:);

    %% Normalize reconstructed images
    s_nom = sum(P.*abs(reco_nom), 'all') / sum(abs(reco_nom).*abs(reco_nom), 'all');
    s_girf = sum(P.*abs(reco_girf), 'all') / sum(abs(reco_girf).*abs(reco_girf), 'all');
    
    recos_nom(:,:,loop) = s_nom*reco_nom;
    recos_girf(:,:,loop) = s_girf*reco_girf;

end % loop

meas1 = measurements{1};
meas2 = measurements{2};
meas3 = measurements{3};

%% Define colors
blue = [0 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];
violet = [0.4940 0.1840 0.5560];
green = [0.4660 0.6740 0.1880];
yellow = [0.9290 0.6940 0.1250];
lblue = [0.3010 0.7450 0.9330];
dred = [0.6350 0.0780 0.1840];

%% Plot reconstructed images
figure('Units','centimeters','Position',[0 0 17.56 24],'Color','k');
colormap gray;

c_lims1 = [0 1.2];
c_lims2 = [-0.2 1];
c_lims3 = [-0.3 0.3];

% 96 interleaves, phantom
ax1 = subplot(5,3,1);
imagesc(P, c_lims1);
axis image;
xticklabels([]); yticklabels([]);
title('96 interleaves','Color','w');
ylabel('Phantom','Color','w','FontWeight','bold');
text(0,10,'(A)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');
colorbar('Color','w');

% 16 interleaves, phantom
ax2 = subplot(5,3,2);
imagesc(P, c_lims1);
axis image;
xticklabels([]); yticklabels([]);
title('16 interleaves','Color','w');
text(0,10,'(B)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');
colorbar('Color','w');

% 3 interleaves, phantom
ax3 = subplot(5,3,3);
imagesc(P, c_lims1);
axis image;
xticklabels([]); yticklabels([]);
title('3 interleaves','Color','w');
text(0,10,'(C)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');
colorbar('Color','w');

% 96 interleaves, nominal recon
ax4 = subplot(5,3,4);
imagesc(abs(recos_nom(:,:,1)), c_lims1);
axis image;
xticklabels([]); yticklabels([]);
ylabel('nominal','Color','w','FontWeight','bold');
text(0,10,'(D)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');
colorbar('Color','w');

% 16 interleaves, nominal recon
ax5 = subplot(5,3,5);
imagesc(abs(recos_nom(:,:,2)), c_lims1);
axis image;
xticklabels([]); yticklabels([]);
text(0,10,'(E)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');
colorbar('Color','w');

% 3 interleaves, nominal recon
ax6 = subplot(5,3,6);
imagesc(abs(recos_nom(:,:,3)), c_lims1);
axis image;
xticklabels([]); yticklabels([]);
text(0,10,'(F)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');
colorbar('Color','w');

% 96 interleaves, Phantom - nominal recon
ax7 = subplot(5,3,7);
imagesc(P-abs(recos_nom(:,:,1)), c_lims2);
axis image;
xticklabels([]); yticklabels([]);
ylabel('Phantom - nominal recon','Color','w','FontWeight','bold');
text(0,10,'(G)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');
colorbar('Color','w');

% 16 interleaves, Phantom - nominal recon
ax8 = subplot(5,3,8);
imagesc(P-abs(recos_nom(:,:,2)), c_lims2);
axis image;
xticklabels([]); yticklabels([]);
text(0,10,'(H)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');
colorbar('Color','w');

% 3 interleaves, Phantom - nominal recon
ax9 = subplot(5,3,9);
imagesc(P-abs(recos_nom(:,:,3)), c_lims2);
axis image;
xticklabels([]); yticklabels([]);
text(0,10,'(I)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');
colorbar('Color','w');

% 96 interleaves, GSTF recon
ax10 = subplot(5,3,10);
imagesc(abs(recos_girf(:,:,1)), c_lims1);
axis image;
xticklabels([]); yticklabels([]);
ylabel('corrected','Color','w','FontWeight','bold');
text(0,10,'(G)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');
colorbar('Color','w');

% 16 interleaves, GSTF recon
ax11 = subplot(5,3,11);
imagesc(abs(recos_girf(:,:,2)), c_lims1);
axis image;
xticklabels([]); yticklabels([]);
text(0,10,'(H)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');
colorbar('Color','w');

% 3 interleaves, GSTF recon
ax12 = subplot(5,3,12);
imagesc(abs(recos_girf(:,:,3)), c_lims1);
axis image;
xticklabels([]); yticklabels([]);
text(0,10,'(I)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');
colorbar('Color','w');

% 96 interleaves, Phantom - GSTF recon
ax13 = subplot(5,3,13);
imagesc(P-abs(recos_girf(:,:,1)), c_lims3);
axis image;
xticklabels([]); yticklabels([]);
ylabel('Phantom - corrected recon','Color','w','FontWeight','bold');
text(0,10,'(J)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');
colorbar('Color','w');

% 16 interleaves, Phantom - GSTF recon
ax14 = subplot(5,3,14);
imagesc(P-abs(recos_girf(:,:,2)), c_lims3);
axis image;
xticklabels([]); yticklabels([]);
text(0,10,'(K)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');
colorbar('Color','w');

% 3 interleaves, Phantom - GSTF recon
ax15 = subplot(5,3,15);
imagesc(P-abs(recos_girf(:,:,3)), c_lims3);
axis image;
xticklabels([]); yticklabels([]);
text(0,10,'(L)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');
colorbar('Color','w');

set(gcf, 'InvertHardcopy', 'off');







