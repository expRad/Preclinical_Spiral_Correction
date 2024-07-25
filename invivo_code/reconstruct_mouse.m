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

filepath = "..\invivo_data\spiral_data_m96.mat"; % in vivo measurement
filepath_traj = "..\invivo_data\spiral_data_traj96.mat"; % trajectory measurement in spherical phantom
filepath_cart = "..\invivo_data\cartesian_data_m.mat"; % Cartesian measurement

delay_girf = -2.9*1e-6;
delay_del = -17.8*1e-6;

%% load GSTFs
H_x = load('..\GSTF_data\Hx_sr100_20231009.mat');
H_y = load('..\GSTF_data\Hy_sr100_20231009.mat');
H_z = load('..\GSTF_data\Hz_sr100_20231009.mat');

%%
H_1 = H_x.H_matrix;
H_2 = H_y.H_matrix;
H_3 = H_z.H_matrix;
gstf_x = H_x.H_matrix.gstf(:,2);
gstf_y = H_y.H_matrix.gstf(:,2);
gstf_z = H_z.H_matrix.gstf(:,2);

gstf_xwdel = H_x.H_matrix.gstf(:,2) .* exp(1i*2*pi*H_x.H_matrix.f_axis.' *delay_girf);
gstf_ywdel = H_y.H_matrix.gstf(:,2) .* exp(1i*2*pi*H_y.H_matrix.f_axis.' *delay_girf);
gstf_zwdel = H_z.H_matrix.gstf(:,2) .* exp(1i*2*pi*H_z.H_matrix.f_axis.' *delay_girf);

dwelltime_girf = H_x.H_matrix.dt;

gamma = 267.513*10^6; %Hz/T

%%
blue = [0 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];
violet = [0.4940 0.1840 0.5560];
green = [0.4660 0.6740 0.1880];
yellow = [0.9290 0.6940 0.1250];
lblue = [0.3010 0.7450 0.9330];
dred = [0.6350 0.0780 0.1840];

%% Read data
disp('Read measurement data.')

if isfile(filepath)
    load(filepath);
else
    disp(['Error: cannot find ',filepath]);
end
if isfile(filepath_traj)
    spiral_data_traj = load(filepath_traj);
    spiral_data_traj = spiral_data_traj.spiral_data;
else
    disp(['Error: cannot find ',filepath_traj]);
end
if isfile(filepath_cart)
    load(filepath_cart);
else
    disp(['Error: cannot find ',filepath_cart]);
end
raw = spiral_data.raw;
nPost = spiral_data.nPost; % data points discarded for reconstruction

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
t_ADC = (0:1:(size(raw,1)-1))*dwelltime_meas + dwelltime_meas/2;
t_GRT = ((1:1:(size(spiral_x,1)))-0.5)*grad_raster_time;
t_eps = 1e-12;
t_GRT_00 = [-t_eps, t_GRT(1)-t_eps, t_GRT, t_GRT(end)+t_eps, t_GRT(end)+grad_raster_time+t_eps];

spiral_x_00 = [zeros(2,size(spiral_x,2)); spiral_x; zeros(2,size(spiral_x,2))];
spiral_y_00 = [zeros(2,size(spiral_x,2)); spiral_y; zeros(2,size(spiral_x,2))];

spiral_x_pp = interp1(t_GRT_00, spiral_x_00*max_grad_strength, 'linear', 'pp');
spiral_y_pp = interp1(t_GRT_00, spiral_y_00*max_grad_strength, 'linear', 'pp');

traj_x_pp = fnint(spiral_x_pp);
traj_y_pp = fnint(spiral_y_pp);

%% reconstruct with nominal trajectory
trajx_nom = gamma/(2*pi)*ppval(traj_x_pp, t_ADC)/1000;
trajy_nom = gamma/(2*pi)*ppval(traj_y_pp, t_ADC)/1000;
traj_nom = trajx_nom.' + 1i*trajy_nom.';

FT = cGrid(traj_nom(1:end-nPost,:)/traj_factor,matrix_size);
reco_nom = FT'*raw(1:end-nPost,:);

%% interpolate to gstf time grid (so the frequency range matches)
t_GIRF = H_1.t_axis;
grad_x_tgirf = ppval(spiral_x_pp, t_GIRF);
grad_y_tgirf = ppval(spiral_y_pp, t_GIRF);

%% rotate to XYZ coordinate system
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
grad_px_nom = ppval(grad_physicalX_pp, t_ADC);
grad_py_nom = ppval(grad_physicalY_pp, t_ADC);

%% zerofill gradients to avoid side effects by the fft calculation
nExtra = round((1e6-size(grad_x_tgirf,1))/2);
grad_x_ZF = cat(1, zeros(nExtra,size(spiral_x,2)), grad_physicalX, zeros(nExtra,size(spiral_x,2)));
grad_y_ZF = cat(1, zeros(nExtra,size(spiral_y,2)), grad_physicalY, zeros(nExtra,size(spiral_x,2)));
grad_z_ZF = cat(1, zeros(nExtra,size(spiral_x,2)), grad_physicalZ, zeros(nExtra,size(spiral_x,2)));

gstf_x_interp = interp1(H_1.f_axis, gstf_x, linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).', 'linear');
gstf_y_interp = interp1(H_2.f_axis, gstf_y, linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).', 'linear');
gstf_z_interp = interp1(H_3.f_axis, gstf_z, linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).', 'linear');

%% calculate post-correction with GSTF correction
tmp1 = fft1d(grad_x_ZF,1);
tmp2 = repmat(gstf_x_interp, [1,spiral_data.numInterleaves]);
grad_px_corr = real(ifft1d( tmp1.*tmp2, 1)) + H_x.H_matrix.fieldOffsets(2);
grad_px_corr = grad_px_corr(nExtra+1:end-nExtra,:);

tmp1 = fft1d(grad_y_ZF,1);
tmp2 = repmat(gstf_y_interp, [1,spiral_data.numInterleaves]);
grad_py_corr = real(ifft1d( tmp1.*tmp2, 1)) + H_y.H_matrix.fieldOffsets(2);
grad_py_corr = grad_py_corr(nExtra+1:end-nExtra,:);

tmp1 = fft1d(grad_z_ZF,1);
tmp2 = repmat(gstf_z_interp, [1,spiral_data.numInterleaves]);
grad_pz_corr = real(ifft1d( tmp1.*tmp2, 1)) + H_z.H_matrix.fieldOffsets(2);
grad_pz_corr = grad_pz_corr(nExtra+1:end-nExtra,:);

grad_px_corr_pp = interp1(t_GIRF, grad_px_corr, 'linear','pp');
grad_py_corr_pp = interp1(t_GIRF, grad_py_corr, 'linear','pp');
grad_px_girf = ppval(grad_px_corr_pp, t_ADC);
grad_py_girf = ppval(grad_py_corr_pp, t_ADC);

% rotate back to PRS (xyz/ in-plane) coordinate system
grad_x_girf = zeros(size(grad_x_tgirf));
grad_y_girf = zeros(size(grad_x_tgirf));

for arm = 1:spiral_data.numInterleaves
    gradXYZ = [squeeze(grad_px_corr(:,arm)), squeeze(grad_py_corr(:,arm)), squeeze(grad_pz_corr(:,arm))].';
    grad_xyz = spiral_data.gradMatrix * gradXYZ;
    grad_x_girf(:,arm) = grad_xyz(1,:);
    grad_y_girf(:,arm) = grad_xyz(2,:);
end

% interpolate back to ADC time grid for reconstruction
grad_x_girf_pp = interp1(t_GIRF, grad_x_girf, 'linear','pp');
grad_y_girf_pp = interp1(t_GIRF, grad_y_girf, 'linear','pp');
traj_x_girf_pp = fnint(grad_x_girf_pp);
traj_y_girf_pp = fnint(grad_y_girf_pp);

%% reconstruct with GSTF-corrected trajectory
trajx_corr = gamma/(2*pi)*ppval(traj_x_girf_pp, t_ADC)/1000;
trajy_corr = gamma/(2*pi)*ppval(traj_y_girf_pp, t_ADC)/1000;
traj_corr = trajx_corr.' + 1i*trajy_corr.';

FT = cGrid(traj_corr(1:end-nPost,:)/traj_factor,matrix_size);
reco_girf = FT'*raw(1:end-nPost,:);

%% calculate post-correction with GSTF+delay correction
gstf_xwdel_interp = interp1(H_1.f_axis, gstf_xwdel, linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).', 'linear');
gstf_ywdel_interp = interp1(H_2.f_axis, gstf_ywdel, linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).', 'linear');
gstf_zwdel_interp = interp1(H_3.f_axis, gstf_zwdel, linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).', 'linear');

tmp1 = fft1d(grad_x_ZF,1);
tmp2 = repmat(gstf_xwdel_interp, [1,spiral_data.numInterleaves]);
grad_px_corr = real(ifft1d( tmp1.*tmp2, 1)) + H_x.H_matrix.fieldOffsets(2);
grad_px_corr = grad_px_corr(nExtra+1:end-nExtra,:);

tmp1 = fft1d(grad_y_ZF,1);
tmp2 = repmat(gstf_ywdel_interp, [1,spiral_data.numInterleaves]);
grad_py_corr = real(ifft1d( tmp1.*tmp2, 1)) + H_y.H_matrix.fieldOffsets(2);
grad_py_corr = grad_py_corr(nExtra+1:end-nExtra,:);

tmp1 = fft1d(grad_z_ZF,1);
tmp2 = repmat(gstf_zwdel_interp, [1,spiral_data.numInterleaves]);
grad_pz_corr = real(ifft1d( tmp1.*tmp2, 1)) + H_z.H_matrix.fieldOffsets(2);
grad_pz_corr = grad_pz_corr(nExtra+1:end-nExtra,:);

grad_px_corr_pp = interp1(t_GIRF, grad_px_corr, 'linear','pp');
grad_py_corr_pp = interp1(t_GIRF, grad_py_corr, 'linear','pp');
grad_px_girfdel = ppval(grad_px_corr_pp, t_ADC);
grad_py_girfdel = ppval(grad_py_corr_pp, t_ADC);

% rotate back to PRS (xyz/ in-plane) coordinate system
grad_x_girf = zeros(size(grad_x_tgirf));
grad_y_girf = zeros(size(grad_x_tgirf));

for arm = 1:spiral_data.numInterleaves
    gradXYZ = [squeeze(grad_px_corr(:,arm)), squeeze(grad_py_corr(:,arm)), squeeze(grad_pz_corr(:,arm))].';
    grad_xyz = spiral_data.gradMatrix * gradXYZ;
    grad_x_girf(:,arm) = grad_xyz(1,:);
    grad_y_girf(:,arm) = grad_xyz(2,:);
end

% interpolate back to ADC time grid for reconstruction
grad_x_girfdel_pp = interp1(t_GIRF, grad_x_girf, 'linear','pp');
grad_y_girfdel_pp = interp1(t_GIRF, grad_y_girf, 'linear','pp');
traj_x_girfdel_pp = fnint(grad_x_girfdel_pp);
traj_y_girfdel_pp = fnint(grad_y_girfdel_pp);

%% reconstruct with GSTF+delay-corrected trajectory
trajx_girfdel = gamma/(2*pi)*ppval(traj_x_girfdel_pp, t_ADC)/1000;
trajy_girfdel = gamma/(2*pi)*ppval(traj_y_girfdel_pp, t_ADC)/1000;
traj_girfdel = trajx_girfdel.' + 1i*trajy_girfdel.';

FT = cGrid(traj_girfdel(1:end-nPost,:)/traj_factor,matrix_size);
reco_girf_delay = FT'*raw(1:end-nPost,:);

%% calculate post-correction with isotropic delay
gstf_x_del = exp(1i*2*pi*linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).' *delay_del);
gstf_y_del = exp(1i*2*pi*linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).' *delay_del);
gstf_z_del = exp(1i*2*pi*linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).' *delay_del);

tmp1 = fft1d(grad_x_ZF,1);
tmp2 = repmat(gstf_x_del, [1,spiral_data.numInterleaves]);
grad_px_corr = real(ifft1d( tmp1.*tmp2, 1));
grad_px_corr = grad_px_corr(nExtra+1:end-nExtra,:);

tmp1 = fft1d(grad_y_ZF,1);
tmp2 = repmat(gstf_y_del, [1,spiral_data.numInterleaves]);
grad_py_corr = real(ifft1d( tmp1.*tmp2, 1));
grad_py_corr = grad_py_corr(nExtra+1:end-nExtra,:);

tmp1 = fft1d(grad_z_ZF,1);
tmp2 = repmat(gstf_z_del, [1,spiral_data.numInterleaves]);
grad_pz_corr = real(ifft1d( tmp1.*tmp2, 1));
grad_pz_corr = grad_pz_corr(nExtra+1:end-nExtra,:);

grad_px_corr_pp = interp1(t_GIRF, grad_px_corr, 'linear','pp');
grad_py_corr_pp = interp1(t_GIRF, grad_py_corr, 'linear','pp');
grad_px_del = ppval(grad_px_corr_pp, t_ADC);
grad_py_del = ppval(grad_py_corr_pp, t_ADC);
    
% rotate back to PRS (xyz/ in-plane) coordinate system
grad_x_girf = zeros(size(grad_x_tgirf));
grad_y_girf = zeros(size(grad_x_tgirf));

for arm = 1:spiral_data.numInterleaves
    gradXYZ = [squeeze(grad_px_corr(:,arm)), squeeze(grad_py_corr(:,arm)), squeeze(grad_pz_corr(:,arm))].';
    grad_xyz = spiral_data.gradMatrix * gradXYZ;
    grad_x_girf(:,arm) = grad_xyz(1,:);
    grad_y_girf(:,arm) = grad_xyz(2,:);
end

% interpolate back to ADC time grid for reconstruction
grad_x_del_pp = interp1(t_GIRF, grad_x_girf, 'linear','pp');
grad_y_del_pp = interp1(t_GIRF, grad_y_girf, 'linear','pp');
traj_x_del_pp = fnint(grad_x_del_pp);
traj_y_del_pp = fnint(grad_y_del_pp);

%% reconstruct with delay-corrected trajectory
trajx_del = gamma/(2*pi)*ppval(traj_x_del_pp, t_ADC)/1000;
trajy_del = gamma/(2*pi)*ppval(traj_y_del_pp, t_ADC)/1000;
traj_del = trajx_del.' + 1i*trajy_del.';

FT = cGrid(traj_del(1:end-nPost,:)/traj_factor,matrix_size);
reco_delay = FT'*raw(1:end-nPost,:);

%% reconstruct images with measured trajectory (measured in situ)
trajx_meas = reshape(spiral_data.trajKx_meas, [spiral_data.datapoints, spiral_data.numInterleaves]);
trajy_meas = reshape(spiral_data.trajKy_meas, [spiral_data.datapoints, spiral_data.numInterleaves]);
trajz_meas = reshape(spiral_data.trajKz_meas, [spiral_data.datapoints, spiral_data.numInterleaves]);

traj_inPlaneX = zeros(size(trajx_meas));
traj_inPlaneY = zeros(size(trajx_meas));
traj_inPlaneZ = zeros(size(trajx_meas));
for arm = 1:spiral_data.numInterleaves
    traj_xyz_meas = [squeeze(trajx_meas(:,arm)), squeeze(trajy_meas(:,arm)), squeeze(trajz_meas(:,arm))].';
    trajXYZ_inPlane = spiral_data.gradMatrix * traj_xyz_meas;
    traj_inPlaneX(:,arm) = trajXYZ_inPlane(1,:);
    traj_inPlaneY(:,arm) = trajXYZ_inPlane(2,:);
    traj_inPlaneZ(:,arm) = trajXYZ_inPlane(3,:);
end
traj_measmouse = traj_inPlaneX + 1i*traj_inPlaneY;

FT = cGrid(traj_measmouse(1:end-nPost,:)/traj_factor,matrix_size);
reco_measmouse = FT'*raw(1:end-nPost,:);

trajx_meas_pp = interp1(t_ADC, trajx_meas, 'linear','pp');
trajy_meas_pp = interp1(t_ADC, trajy_meas, 'linear','pp');
gradx_meas_pp = fnder(trajx_meas_pp);
grady_meas_pp = fnder(trajy_meas_pp);
grad_x_measmouse = ppval(gradx_meas_pp, t_ADC).'/gamma*2*pi*1000;
grad_y_measmouse = ppval(grady_meas_pp, t_ADC).'/gamma*2*pi*1000;

%% reconstruct images with measured trajectory (measured in phantom)
trajx_meas = reshape(spiral_data_traj.trajKx_meas, [spiral_data.datapoints, spiral_data.numInterleaves]);
trajy_meas = reshape(spiral_data_traj.trajKy_meas, [spiral_data.datapoints, spiral_data.numInterleaves]);
trajz_meas = reshape(spiral_data_traj.trajKz_meas, [spiral_data.datapoints, spiral_data.numInterleaves]);

traj_inPlaneX = zeros(size(trajx_meas));
traj_inPlaneY = zeros(size(trajx_meas));
traj_inPlaneZ = zeros(size(trajx_meas));
for arm = 1:spiral_data.numInterleaves
    traj_xyz_meas = [squeeze(trajx_meas(:,arm)), squeeze(trajy_meas(:,arm)), squeeze(trajz_meas(:,arm))].';
    trajXYZ_inPlane = spiral_data.gradMatrix * traj_xyz_meas;
    traj_inPlaneX(:,arm) = trajXYZ_inPlane(1,:);
    traj_inPlaneY(:,arm) = trajXYZ_inPlane(2,:);
    traj_inPlaneZ(:,arm) = trajXYZ_inPlane(3,:);
end
traj_measphan = traj_inPlaneX + 1i*traj_inPlaneY;

FT = cGrid(traj_measphan(1:end-nPost,:)/traj_factor,matrix_size);
reco_measphan = FT'*raw(1:end-nPost,:);

trajx_meas_pp = interp1(t_ADC, trajx_meas, 'linear','pp');
trajy_meas_pp = interp1(t_ADC, trajy_meas, 'linear','pp');
gradx_meas_pp = fnder(trajx_meas_pp);
grady_meas_pp = fnder(trajy_meas_pp);
grad_x_measphan = (ppval(gradx_meas_pp, t_ADC)).'/gamma*2*pi*1000;
grad_y_measphan = (ppval(grady_meas_pp, t_ADC)).'/gamma*2*pi*1000;

%% Cartesian image
frame = 1;
reco_cart = ifft2d(squeeze(cartesian_data.raw(:,:,frame)));

%% Plot reconstructed images (Figure 8)
yellowColMap = [linspace(255, 0, 124)', linspace(255, 0, 124)', zeros(124, 1); zeros(132, 3)];
blueColMap = [zeros(132, 3); zeros(124, 1), linspace(0, 255, 124)', linspace(0, 255, 124)'];
myColorMap = uint8(blueColMap + yellowColMap);


figure('Units','centimeters','Position',[0 0 17.56 22],'Color','k');
colormap gray;
c_lims1 = [0 0.2];
c_lims2 = [-0.06 0.06];

dx = 0.29;
dy = 0.21;
x1 = 0.03;
x2 = 0.32;
x3 = 0.61;
y1 = 0.75;
y2 = 0.52;
y3 = 0.27;
y4 = 0.03;

% reco nom
ax1 = subplot('Position',[x1 y1 dx dy]);
imagesc(abs(reco_nom),c_lims1);
axis image; axis off;
title('nominal','Color','w');
set(gca,'FontName','Times','Fontsize',9);
text(-20,96, 'images','Color','w','FontWeight','bold','HorizontalAlignment','center','Rotation',90,'FontName','Times','Fontsize',10);
text(0,10,'(A)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% reco delay
ax2 = subplot('Position',[x2 y1 dx dy]);
imagesc(abs(reco_delay),c_lims1);
axis image; axis off;
title('delay','Color','w');
set(gca,'FontName','Times','Fontsize',9);
text(0,10,'(B)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% reco gstf
ax3 = subplot('Position',[x3 y1 dx dy]);
imagesc(abs(reco_girf),c_lims1);
axis image; axis off;
colorbar('Position',[0.92 y1+0.01 0.02 dy-0.02], 'Color','w');
title('GSTF','Color','w');
set(gca,'FontName','Times','Fontsize',9);
text(0,10,'(C)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% diff nom
ax4 = subplot('Position',[x1 y2 dx dy]);
imagesc(abs(reco_measphan)-abs(reco_nom),c_lims2);
colormap(ax4, myColorMap);
axis image; axis off;
text(-20,96, 'difference to (I)','Color','w','FontWeight','bold','HorizontalAlignment','center','Rotation',90,'FontName','Times','Fontsize',10);
text(0,10,'(D)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% diff delay
ax5 = subplot('Position',[x2 y2 dx dy]);
imagesc(abs(reco_measphan)-abs(reco_delay),c_lims2);
colormap(ax5, myColorMap);
axis image; axis off;
text(0,10,'(E)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% diff gstf
ax6 = subplot('Position',[x3 y2 dx dy]);
imagesc(abs(reco_measphan)-abs(reco_girf),c_lims2);
colormap(ax6, myColorMap);
axis image; axis off;
colorbar('Position',[0.92 y2+0.01 0.02 dy-0.02], 'Color','w');
set(gca,'FontName','Times','Fontsize',9);
text(0,10,'(F)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% reco gstfdel
ax7 = subplot('Position',[x1 y3 dx dy]);
imagesc(abs(reco_girf_delay),c_lims1);
% colorbar;
axis image; axis off;
title('GSTF+delay','Color','w');
set(gca,'FontName','Times','Fontsize',9);
text(-20,96, 'images','Color','w','FontWeight','bold','HorizontalAlignment','center','Rotation',90,'FontName','Times','Fontsize',10);
text(0,10,'(G)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% reco measmouse
ax8 = subplot('Position',[x2 y3 dx dy]);
imagesc(abs(reco_measmouse),c_lims1);
% colorbar;
axis image; axis off;
title('measured in vivo','Color','w');
set(gca,'FontName','Times','Fontsize',9);
text(0,10,'(H)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% reco measphan
ax9 = subplot('Position',[x3 y3 dx dy]);
imagesc(abs(reco_measphan),c_lims1);
axis image; axis off;
colorbar('Position',[0.92 y3+0.01 0.02 dy-0.02], 'Color','w');
title('measured in phantom','Color','w');
set(gca,'FontName','Times','Fontsize',9);
text(0,10,'(I)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% diff gstfdel
ax10 = subplot('Position',[x1 y4 dx dy]);
imagesc(abs(reco_measphan)-abs(reco_girf_delay),c_lims2);
colormap(ax10, myColorMap);
axis image; axis off;
text(-20,96, 'difference to (I)','Color','w','FontWeight','bold','HorizontalAlignment','center','Rotation',90,'FontName','Times','Fontsize',10);
text(0,10,'(J)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% diff measmouse
ax11 = subplot('Position',[x2 y4 dx dy]);
imagesc(abs(reco_measphan)-abs(reco_measmouse),c_lims2);
colormap(ax11, myColorMap);
axis image; axis off;
colorbar('Position',[0.62 y4+0.01 0.02 dy-0.02], 'Color','w');
set(gca,'FontName','Times','Fontsize',9);
text(0,10,'(K)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% diff measphan
ax12 = subplot('Position',[x3+0.09 y4-0.0 dx dy]);
imagesc(abs(reco_cart));
axis image; axis off;
title('Cartesian acquisition', 'Color','w');
set(gca,'FontName','Times','Fontsize',9);
text(0,10,'(L)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

set(gcf, 'InvertHardcopy', 'off');

%% plot gradients and trajectory (Figure 9)
arm = 1;

figure('Units','centimeters','Position',[0 0 17.56 22]);

ax1 = subplot('Position',[0.08 0.805 0.4 0.13]);
plot(t_ADC*1000, grad_x_measmouse(:,arm), 'Color',green,'LineWidth',1.2);
hold on;
plot(t_ADC*1000, grad_px_nom(:,arm), 'Color',blue,'LineWidth',1.2);
plot(t_ADC*1000, grad_px_del(:,arm), 'Color',orange,'LineWidth',1.2);
plot(t_ADC*1000, grad_px_girf(:,arm), 'Color',yellow,'LineWidth',1.2);
plot(t_ADC*1000, grad_px_girfdel(:,arm), 'Color',violet,'LineWidth',1.2);
xlabel('Time (ms)');
ylabel('Gradient_x (T/m)');
ylim([-0.2 0.2]);
text(1.15, 0.25,'trajectory measured in vivo', 'HorizontalAlignment','center','FontName','Times','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Times','Fontsize',9);
text(-0.4,0.25,'(A)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(-0.4,-0.35,'(C)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(-0.4,-1.5,'(E)','FontName','Arial','Fontsize',11,'FontWeight','bold');

ax2 = subplot('Position',[0.575 0.805 0.4 0.13]);
plot(t_ADC*1000, grad_px_nom(:,arm), 'Color',blue,'LineWidth',1.2,'DisplayName','nominal');
hold on;
plot(t_ADC*1000, grad_px_del(:,arm), 'Color',orange,'LineWidth',1.2,'DisplayName','delay');
plot(t_ADC*1000, grad_px_girf(:,arm), 'Color',yellow,'LineWidth',1.2,'DisplayName','GSTF');
plot(t_ADC*1000, grad_px_girfdel(:,arm), 'Color',violet,'LineWidth',1.2,'DisplayName','GSTF+delay');
plot(t_ADC*1000, grad_x_measphan(:,arm), 'Color',green,'LineWidth',1.2,'DisplayName','measured');
legend('numColumns',5,'Position',[0.08 0.97 0.895 0.025]);
xlabel('Time (ms)');
ylim([-0.2 0.2]);
text(1.15, 0.25,'trajectory measured in phantom', 'HorizontalAlignment','center','FontName','Times','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Times','Fontsize',9);
text(-0.4,0.25,'(B)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(-0.4,-0.35,'(D)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(-0.4,-1.5,'(F)','FontName','Arial','Fontsize',11,'FontWeight','bold');

ax3 = subplot('Position',[0.08 0.415 0.4 0.35]);
plot(traj_nom(:,arm),'LineWidth',1.2);
hold on;
plot(traj_del(:,arm),'LineWidth',1.2);
plot(traj_corr(:,arm),'LineWidth',1.2);
plot(traj_girfdel(:,arm),'LineWidth',1.2);
plot(traj_measmouse(:,arm),'LineWidth',1.2);
axis image;
xlim([-2 4]);
ylim([-4 2]-0.7);
xlabel('k_x (1/m)');
ylabel('k_y (1/m)');
set(gca,'FontName','Times','Fontsize',9);
rectangle('Position',[-0.05 -0.1 0.5 0.4],'LineStyle','-.','LineWidth',1.2);
rectangle('Position',[0.25 -2.6 1.55 0.6],'LineStyle','--','LineWidth',1.2);
% inset 1 measured in vivo
ax3a = axes('Position',[0.115 0.45 0.13 0.08]);
plot(traj_nom(:,arm),'LineWidth',1.2);
hold on;
plot(traj_del(:,arm),'LineWidth',1.2);
plot(traj_corr(:,arm),'LineWidth',1.2);
plot(traj_girfdel(:,arm),'LineWidth',1.2);
plot(traj_measmouse(:,arm),'LineWidth',1.2);
xlim([-0.1 0.4]+0.05);
ylim([-0.1 0.3]);
box off;
rectangle('Position',[-0.045 -0.095 0.49 0.39],'LineStyle','-.','LineWidth',1.2);
% inset 2 measured in vivo
ax3b = axes('Position',[0.3 0.45 0.17 0.08]);
plot(traj_nom(:,arm),'LineWidth',1.2);
hold on;
plot(traj_del(:,arm),'LineWidth',1.2);
plot(traj_corr(:,arm),'LineWidth',1.2);
plot(traj_girfdel(:,arm),'LineWidth',1.2);
plot(traj_measmouse(:,arm),'LineWidth',1.2);
xlim([0.25 1.8]);
ylim([-2.6 -2]);
box off;
rectangle('Position',[0.26 -2.59 1.53 0.58],'LineStyle','--','LineWidth',1.2);

ax4 = subplot('Position',[0.575 0.415 0.4 0.35]);
plot(traj_nom(:,arm),'LineWidth',1.2);
hold on;
plot(traj_del(:,arm),'LineWidth',1.2);
plot(traj_corr(:,arm),'LineWidth',1.2);
plot(traj_girfdel(:,arm),'LineWidth',1.2);
plot(traj_measphan(:,arm),'LineWidth',1.2);
axis image;
xlim([-2 4]);
ylim([-4 2]-0.7);
xlabel('k_x (1/m)');
set(gca,'FontName','Times','Fontsize',9);
rectangle('Position',[-0.05 -0.1 0.5 0.4],'LineStyle','-.','LineWidth',1.2);
rectangle('Position',[0.25 -2.6 1.55 0.6],'LineStyle','--','LineWidth',1.2);
% inset 1 measured in phantom
ax4a = axes('Position',[0.61 0.45 0.13 0.08]);
plot(traj_nom(:,arm),'LineWidth',1.2);
hold on;
plot(traj_del(:,arm),'LineWidth',1.2);
plot(traj_corr(:,arm),'LineWidth',1.2);
plot(traj_girfdel(:,arm),'LineWidth',1.2);
plot(traj_measphan(:,arm),'LineWidth',1.2);
xlim([-0.1 0.4]+0.05);
ylim([-0.1 0.3]);
box off;
rectangle('Position',[-0.045 -0.095 0.49 0.39],'LineStyle','-.','LineWidth',1.2);
% inset 2 measured in phantom
ax4b = axes('Position',[0.795 0.45 0.17 0.08]);
plot(traj_nom(:,arm),'LineWidth',1.2);
hold on;
plot(traj_del(:,arm),'LineWidth',1.2);
plot(traj_corr(:,arm),'LineWidth',1.2);
plot(traj_girfdel(:,arm),'LineWidth',1.2);
plot(traj_measphan(:,arm),'LineWidth',1.2);
xlim([0.25 1.8]);
ylim([-2.6 -2]);
box off;
rectangle('Position',[0.26 -2.59 1.53 0.58],'LineStyle','--','LineWidth',1.2);

ax5 = subplot('Position',[0.08 0.04 0.4 0.35]);
plot(traj_measmouse(:,arm)-traj_nom(:,arm),'LineWidth',1.2);
hold on;
plot(traj_measmouse(:,arm)-traj_del(:,arm),'LineWidth',1.2);
plot(traj_measmouse(:,arm)-traj_corr(:,arm),'LineWidth',1.2);
plot(traj_measmouse(:,arm)-traj_girfdel(:,arm),'LineWidth',1.2);
xlabel('\Delta k_x (1/m)');
ylabel('\Delta k_y (1/m)');
axis image;
xlim([-0.2 0.2]);
ylim([-0.2 0.2]);
set(gca,'FontName','Times','Fontsize',9);

ax6 = subplot('Position',[0.575 0.04 0.4 0.35]);
plot(traj_measphan(:,arm)-traj_nom(:,arm),'LineWidth',1.2);
hold on;
plot(traj_measphan(:,arm)-traj_del(:,arm),'LineWidth',1.2);
plot(traj_measphan(:,arm)-traj_corr(:,arm),'LineWidth',1.2);
plot(traj_measphan(:,arm)-traj_girfdel(:,arm),'LineWidth',1.2);
xlabel('\Delta k_x (1/m)');
axis image;
xlim([-0.2 0.2]);
ylim([-0.2 0.2]);
set(gca,'FontName','Times','Fontsize',9);

linkaxes([ax1 ax2],'x');
xlim(ax1, [0 2.3]);
%% Calculate RMSEs of the trajectories
rmse_nom = sqrt(mean(abs(traj_measphan(:,arm)-traj_nom(:,arm)).^2));
disp(['rmse_nom = ',num2str(rmse_nom)]);
rmse_del = sqrt(mean(abs(traj_measphan(:,arm)-traj_del(:,arm)).^2));
disp(['rmse_del = ',num2str(rmse_del)]);
rmse_girf = sqrt(mean(abs(traj_measphan(:,arm)-traj_corr(:,arm)).^2));
disp(['rmse_girf = ',num2str(rmse_girf)]);
rmse_girfdel = sqrt(mean(abs(traj_measphan(:,arm)-traj_girfdel(:,arm)).^2));
disp(['rmse_girfdel = ',num2str(rmse_girfdel)]);














