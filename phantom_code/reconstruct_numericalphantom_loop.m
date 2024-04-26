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
meas1.delay_girf = -4.0e-6;
meas1.delay = -19.3e-6;

meas2 = struct;
meas2.filepath = '..\phantom_data\spiral_data_p16.mat';
meas2.delay_girf = -4.8e-6;
meas2.delay = -21.1e-6;

meas3 = struct;
meas3.filepath = '..\phantom_data\spiral_data_p03.mat';
meas3.delay_girf = -1.0e-6;
meas3.delay = -17.1e-6;

measurements = {meas1, meas2, meas3};

%% initialize arrays for results
recos_nom = zeros(matrix_size, matrix_size, length(measurements));
recos_girf = zeros(matrix_size, matrix_size, length(measurements));
recos_girfdel = zeros(matrix_size, matrix_size, length(measurements));
recos_del = zeros(matrix_size, matrix_size, length(measurements));
recos_meas = zeros(matrix_size, matrix_size, length(measurements));

%% load numerical phantom
P = phantom('Modified Shepp-Logan',matrix_size);
P(find(P==1)) = 0.8;
P(find(abs(P-0.2)<0.001)) = 0.8;
P(find(abs(P-0.3)<0.001)) = 0.1;

%% load GSTFs
H_x = load('..\GSTF_data\Hx_sr100_20231009.mat');
H_y = load('..\GSTF_data\Hy_sr100_20231009.mat');
H_z = load('..\GSTF_data\Hz_sr100_20231009.mat');

%% Run simulations for every spiral trajectory
for loop = 1:length(measurements)
    disp(['loop = ',num2str(loop)]);
    meas = measurements{loop};
    filepath = meas.filepath;
    
    % Initialize delays and transfer functions
    delay_girf = meas.delay_girf;
    delay = meas.delay;

    delay_y = delay_girf;
    delay_z = delay_girf;
    H_1 = H_x.H_matrix;
    H_2 = H_y.H_matrix;
    H_3 = H_z.H_matrix;
    gstf_x = H_x.H_matrix.gstf(:,2);
    gstf_y = H_y.H_matrix.gstf(:,2);
    gstf_z = H_z.H_matrix.gstf(:,2);
    
    % GSTFs for GSTF+delay correction
    gstf_xwdel = H_x.H_matrix.gstf(:,2) .* exp(1i*2*pi*H_x.H_matrix.f_axis.' *delay_girf);
    gstf_ywdel = H_y.H_matrix.gstf(:,2) .* exp(1i*2*pi*H_y.H_matrix.f_axis.' *delay_y);
    gstf_zwdel = H_z.H_matrix.gstf(:,2) .* exp(1i*2*pi*H_z.H_matrix.f_axis.' *delay_z);

    delay_y2 = delay;
    delay_z2 = delay;

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

    %% reconstruct images with measured trajectory
    disp(['NUMBER of SPIRAL INTERLEAVES: ',num2str(spiral_data.numInterleaves)]);
    disp('    Reconstruction with measured trajectory');
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
    traj_meas = traj_inPlaneX + 1i*traj_inPlaneY;
    measurements{loop}.traj_meas = traj_meas;

    % create synthetic raw data
    FT = cGrid(traj_meas/traj_factor,matrix_size);
    raw = FT*P;
    % reconstruct image
    FT = cGrid(traj_meas(1:end-nPost,:)/traj_factor,matrix_size);
    reco_meas = FT'*raw(1:end-nPost,:);

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

    %% reconstruct with nominal trajectory
    disp('    Reconstruction with nominal trajectory')
    gamma = 267.513*10^6; %Hz/T
    trajx_nom = gamma/(2*pi)*ppval(traj_x_pp, t_ADC)/1000;
    trajy_nom = gamma/(2*pi)*ppval(traj_y_pp, t_ADC)/1000;
    traj_nom = trajx_nom.' + 1i*trajy_nom.';
    measurements{loop}.traj_nom = traj_nom;

    FT = cGrid(traj_nom(1:end-nPost,:)/traj_factor,matrix_size);
    reco_nom = FT'*raw(1:end-nPost,:);

    %% interpolate to gstf time grid (so the frequency range matches)
    t_GIRF = H_1.t_axis;
    if t_GIRF(end) < t_GRT(end)
        t_GIRF_end = ceil(t_GRT(end)/dwelltime_girf);
        t_GIRF = ((0:t_GIRF_end)+0.5)*dwelltime_girf;
    end
    grad_x_tgirf = ppval(spiral_x_pp, t_GIRF);
    grad_y_tgirf = ppval(spiral_y_pp, t_GIRF);

    %% rotate in-plane gradients to XYZ coordinate system
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
    gstf_y_interp = interp1(H_2.f_axis, gstf_y, linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).', 'linear');
    gstf_z_interp = interp1(H_3.f_axis, gstf_z, linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).', 'linear');

    %% calculate post-correction with GSTF
    tmp1 = fft1d(grad_x_ZF,1);
    tmp2 = repmat(gstf_x_interp, [1,spiral_data.numInterleaves]);
    grad_px_girf = real(ifft1d( tmp1.*tmp2, 1)) + H_x.H_matrix.fieldOffsets(2);
    grad_px_girf = grad_px_girf(nExtra+1:end-nExtra,:);

    tmp1 = fft1d(grad_y_ZF,1);
    tmp2 = repmat(gstf_y_interp, [1,spiral_data.numInterleaves]);
    grad_py_girf = real(ifft1d( tmp1.*tmp2, 1)) + H_y.H_matrix.fieldOffsets(2);
    grad_py_girf = grad_py_girf(nExtra+1:end-nExtra,:);

    tmp1 = fft1d(grad_z_ZF,1);
    tmp2 = repmat(gstf_z_interp, [1,spiral_data.numInterleaves]);
    grad_pz_girf = real(ifft1d( tmp1.*tmp2, 1)) + H_z.H_matrix.fieldOffsets(2);
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

    %% reconstruct with GSTF-corrected trajectory
    disp('    Reconstruction with GSTF-corrected trajectory')
    trajx_corr = gamma/(2*pi)*ppval(traj_x_girf_pp, t_ADC)/1000;
    trajy_corr = gamma/(2*pi)*ppval(traj_y_girf_pp, t_ADC)/1000;
    traj_corr = trajx_corr.' + 1i*trajy_corr.';
    measurements{loop}.traj_girf = traj_corr;

    FT = cGrid(traj_corr(1:end-nPost,:)/traj_factor,matrix_size);
    reco_girf = FT'*raw(1:end-nPost,:);

    %% calculate post-correction with GSTF+delay
    gstf_xwdel_interp = interp1(H_1.f_axis, gstf_xwdel, linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).', 'linear');
    gstf_ywdel_interp = interp1(H_2.f_axis, gstf_ywdel, linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).', 'linear');
    gstf_zwdel_interp = interp1(H_3.f_axis, gstf_zwdel, linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).', 'linear');

    tmp1 = fft1d(grad_x_ZF,1);
    tmp2 = repmat(gstf_xwdel_interp, [1,spiral_data.numInterleaves]);
    grad_px_girf = real(ifft1d( tmp1.*tmp2, 1)) + H_x.H_matrix.fieldOffsets(2);
    grad_px_girf = grad_px_girf(nExtra+1:end-nExtra,:);

    tmp1 = fft1d(grad_y_ZF,1);
    tmp2 = repmat(gstf_ywdel_interp, [1,spiral_data.numInterleaves]);
    grad_py_girf = real(ifft1d( tmp1.*tmp2, 1)) + H_y.H_matrix.fieldOffsets(2);
    grad_py_girf = grad_py_girf(nExtra+1:end-nExtra,:);

    tmp1 = fft1d(grad_z_ZF,1);
    tmp2 = repmat(gstf_zwdel_interp, [1,spiral_data.numInterleaves]);
    grad_pz_girf = real(ifft1d( tmp1.*tmp2, 1)) + H_z.H_matrix.fieldOffsets(2);
    grad_pz_girf = grad_pz_girf(nExtra+1:end-nExtra,:);
    
    grad_px_girf_pp = interp1(t_GIRF, grad_px_girf, 'linear','pp');
    grad_py_girf_pp = interp1(t_GIRF, grad_py_girf, 'linear','pp');
    measurements{loop}.grad_x_girfdel = ppval(grad_px_girf_pp, t_ADC);
    measurements{loop}.grad_y_girfdel = ppval(grad_py_girf_pp, t_ADC);

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
    grad_x_girfdel_pp = interp1(t_GIRF, grad_x_girf, 'linear','pp');
    grad_y_girfdel_pp = interp1(t_GIRF, grad_y_girf, 'linear','pp');
    traj_x_girfdel_pp = fnint(grad_x_girfdel_pp);
    traj_y_girfdel_pp = fnint(grad_y_girfdel_pp);

    %% reconstruct with GSTF+delay-corrected trajectory
    disp('    Reconsturction with GSTF+delay-corrected trajectory')
    trajx_girfdel = gamma/(2*pi)*ppval(traj_x_girfdel_pp, t_ADC)/1000;
    trajy_girfdel = gamma/(2*pi)*ppval(traj_y_girfdel_pp, t_ADC)/1000;
    traj_girfdel = trajx_girfdel.' + 1i*trajy_girfdel.';
    measurements{loop}.traj_girfdel = traj_girfdel;

    FT = cGrid(traj_girfdel(1:end-nPost,:)/traj_factor,matrix_size);
    reco_girf_delay = FT'*raw(1:end-nPost,:);

    %% calculate post-correction with isotropic delay
    gstf_x_del = exp(1i*2*pi*linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).' *delay);
    gstf_y_del = exp(1i*2*pi*linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).' *delay_y2);
    gstf_z_del = exp(1i*2*pi*linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).' *delay_z2);

    tmp1 = fft1d(grad_x_ZF,1);
    tmp2 = repmat(gstf_x_del, [1,spiral_data.numInterleaves]);
    grad_px_girf = real(ifft1d( tmp1.*tmp2, 1));
    grad_px_girf = grad_px_girf(nExtra+1:end-nExtra,:);

    tmp1 = fft1d(grad_y_ZF,1);
    tmp2 = repmat(gstf_y_del, [1,spiral_data.numInterleaves]);
    grad_py_girf = real(ifft1d( tmp1.*tmp2, 1));
    grad_py_girf = grad_py_girf(nExtra+1:end-nExtra,:);

    tmp1 = fft1d(grad_z_ZF,1);
    tmp2 = repmat(gstf_z_del, [1,spiral_data.numInterleaves]);
    grad_pz_girf = real(ifft1d( tmp1.*tmp2, 1));
    grad_pz_girf = grad_pz_girf(nExtra+1:end-nExtra,:);
    
    grad_px_girf_pp = interp1(t_GIRF, grad_px_girf, 'linear','pp');
    grad_py_girf_pp = interp1(t_GIRF, grad_py_girf, 'linear','pp');
    measurements{loop}.grad_x_del = ppval(grad_px_girf_pp, t_ADC);
    measurements{loop}.grad_y_del = ppval(grad_py_girf_pp, t_ADC);

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
    grad_x_del_pp = interp1(t_GIRF, grad_x_girf, 'linear','pp');
    grad_y_del_pp = interp1(t_GIRF, grad_y_girf, 'linear','pp');
    traj_x_del_pp = fnint(grad_x_del_pp);
    traj_y_del_pp = fnint(grad_y_del_pp);

    %% reconstruct with delay-corrected trajectory
    disp('    Reconstruction with delay-corrected trajectory')
    trajx_del = gamma/(2*pi)*ppval(traj_x_del_pp, t_ADC)/1000;
    trajy_del = gamma/(2*pi)*ppval(traj_y_del_pp, t_ADC)/1000;
    traj_del = trajx_del.' + 1i*trajy_del.';
    measurements{loop}.traj_del = traj_del;

    FT = cGrid(traj_del(1:end-nPost,:)/traj_factor,matrix_size);
    reco_delay = FT'*raw(1:end-nPost,:);

    %% Derive measured gradient waveforms from measured trajectories for plotting
    trajx_meas_pp = interp1(t_ADC, trajx_meas, 'linear','pp');
    trajy_meas_pp = interp1(t_ADC, trajy_meas, 'linear','pp');
    gradx_meas_pp = fnder(trajx_meas_pp);
    grady_meas_pp = fnder(trajy_meas_pp);
    grad_x_meas = ppval(gradx_meas_pp, t_ADC).'/gamma*2*pi*1000;
    grad_y_meas = ppval(grady_meas_pp, t_ADC).'/gamma*2*pi*1000;
    measurements{loop}.grad_x_meas = grad_x_meas;
    measurements{loop}.grad_y_meas = grad_y_meas;
    
    %% Normalize reconstructed images
    s_nom = sum(abs(reco_meas).*abs(reco_nom), 'all') / sum(abs(reco_nom).*abs(reco_nom), 'all');
    s_girf = sum(abs(reco_meas).*abs(reco_girf), 'all') / sum(abs(reco_girf).*abs(reco_girf), 'all');
    s_girfdel = sum(abs(reco_meas).*abs(reco_girf_delay), 'all') / sum(abs(reco_girf_delay).*abs(reco_girf_delay), 'all');
    s_delay = sum(abs(reco_meas).*abs(reco_delay), 'all') / sum(abs(reco_delay).*abs(reco_delay), 'all');
    
    recos_nom(:,:,loop) = s_nom*reco_nom;
    recos_girf(:,:,loop) = s_girf*reco_girf;
    recos_girfdel(:,:,loop) = s_girfdel*reco_girf_delay;
    recos_del(:,:,loop) = s_delay*reco_delay;
    recos_meas(:,:,loop) = reco_meas;

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

%% Plot reconstructed images (Figure 7)
figure('Units','centimeters','Position',[0 0 17.56 24],'Color','k');
colormap gray;
c_lims1 = [0 3.5e-5];
c_lims2 = [0 3.5e-5];
c_lims3 = [0 3.5e-5];

dx = 0.31;
dy = 0.195;
x1 = 0.06;
x2 = 0.37;
x3 = 0.68;
y1 = 0.78;
y2 = 0.585;
y3 = 0.39;
y4 = 0.195;
y5 = 0;

% 96 interleaves, nominal
ax1 = subplot('Position',[x1 y1 dx dy]);
imagesc(abs(recos_nom(:,:,1)), c_lims1);
axis image; axis off;
xticklabels([]); yticklabels([]);
title('96 interleaves','Color','w');
ylabel('nominal','Color','w','FontWeight','bold');
text(-30,96, 'nominal','Color','w','FontWeight','bold','HorizontalAlignment','center','Rotation',90,'FontName','Times','Fontsize',10);
set(gca,'FontName','Times','Fontsize',9);
text(0,10,'(A)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% 16 interleaves, nominal
ax2 = subplot('Position',[x2 y1 dx dy]);
imagesc(abs(recos_nom(:,:,2)), c_lims2);
axis image; axis off;
title('16 interleaves','Color','w');
set(gca,'FontName','Times','Fontsize',9);
text(0,10,'(B)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% 3 interleaves, nominal
ax3 = subplot('Position',[x3 y1 dx dy]);
imagesc(abs(recos_nom(:,:,3)), c_lims3);
axis image; axis off;
title('3 interleaves','Color','w');
set(gca,'FontName','Times','Fontsize',9);
text(0,10,'(C)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% 96 interleaves, delay
ax4 = subplot('Position',[x1 y2 dx dy]);
imagesc(abs(recos_del(:,:,1)), c_lims1);
axis image; axis off;
xticklabels([]); yticklabels([]);
ylabel('delay','Color','w','FontWeight','bold');
text(-30,96, 'delay','Color','w','FontWeight','bold','HorizontalAlignment','center','Rotation',90,'FontName','Times','Fontsize',10);
set(gca,'FontName','Times','Fontsize',9);
text(0,10,'(D)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% 16 interleaves, delay
ax5 = subplot('Position',[x2 y2 dx dy]);
imagesc(abs(recos_del(:,:,2)), c_lims2);
axis image; axis off;
text(0,10,'(E)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% 3 interleaves, delay
ax6 = subplot('Position',[x3 y2 dx dy]);
imagesc(abs(recos_del(:,:,3)), c_lims3);
axis image; axis off;
text(0,10,'(F)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% 96 interleaves, gstf
ax7 = subplot('Position',[x1 y3 dx dy]);
imagesc(abs(recos_girf(:,:,1)), c_lims1);
axis image; axis off;
xticklabels([]); yticklabels([]);
ylabel('GSTF','Color','w','FontWeight','bold');
text(-30,96, 'GSTF','Color','w','FontWeight','bold','HorizontalAlignment','center','Rotation',90,'FontName','Times','Fontsize',10);
set(gca,'FontName','Times','Fontsize',9);
text(0,10,'(G)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% 16 interleaves, gstf
ax8 = subplot('Position',[x2 y3 dx dy]);
imagesc(abs(recos_girf(:,:,2)), c_lims2);
axis image; axis off;
text(0,10,'(H)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% 3 interleaves, gstf
ax9 = subplot('Position',[x3 y3 dx dy]);
imagesc(abs(recos_girf(:,:,3)), c_lims3);
axis image; axis off;
text(0,10,'(I)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% 96 interleaves, gstf + delay
ax10 = subplot('Position',[x1 y4 dx dy]);
imagesc(abs(recos_girfdel(:,:,1)), c_lims1);
axis image; axis off;
xticklabels([]); yticklabels([]);
ylabel('GSTF+delay','Color','w','FontWeight','bold');
text(-30,96, 'GSTF+delay','Color','w','FontWeight','bold','HorizontalAlignment','center','Rotation',90,'FontName','Times','Fontsize',10);
set(gca,'FontName','Times','Fontsize',9);
text(0,10,'(J)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% 16 interleaves, gstf + delay
ax11 = subplot('Position',[x2 y4 dx dy]);
imagesc(abs(recos_girfdel(:,:,2)), c_lims2);
axis image; axis off;
text(0,10,'(K)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% 3 interleaves, gstf + delay
ax12 = subplot('Position',[x3 y4 dx dy]);
imagesc(abs(recos_girfdel(:,:,3)), c_lims3);
axis image; axis off;
text(0,10,'(L)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% 96 interleaves, measured
ax13 = subplot('Position',[x1 y5 dx dy]);
imagesc(abs(recos_meas(:,:,1)), c_lims1);
axis image; axis off;
xticklabels([]); yticklabels([]);
ylabel('measured','Color','w','FontWeight','bold');
text(-30,96, 'measured','Color','w','FontWeight','bold','HorizontalAlignment','center','Rotation',90,'FontName','Times','Fontsize',10);
set(gca,'FontName','Times','Fontsize',9);
text(0,10,'(M)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% 16 interleaves, measured
ax14 = subplot('Position',[x2 y5 dx dy]);
imagesc(abs(recos_meas(:,:,2)), c_lims2);
axis image; axis off;
text(0,10,'(N)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

% 3 interleaves, measured
ax15 = subplot('Position',[x3 y5 dx dy]);
imagesc(abs(recos_meas(:,:,3)), c_lims3);
axis image; axis off;
text(0,10,'(O)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

set(gcf, 'InvertHardcopy', 'off');

%% Plot the line profiles (Figure 4, but with numerical phantom)
figure('Units','centimeters','Position',[0 0 17.56 11]);
colormap gray;
x_line = 46;
y_line = 86;
i = 1;

ax1 = subplot('Position',[0.055 0.68 0.144 0.23]);
imagesc(abs(recos_meas(:,:,1)), c_lims1);
axis image; axis off;
yline(x_line,'LineWidth',1.5,'Color',blue);
xline(y_line,'LineWidth',1.5,'Color',blue);
text(-20,96, '96 interleaves','FontWeight','bold','HorizontalAlignment','center','Rotation',90,'FontName','Times','Fontsize',10);
i = i+1;
set(gca,'FontName','Times','Fontsize',9);
text(-70,10,'(A)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(200,10,'(B)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(740,10,'(C)','FontName','Arial','Fontsize',11,'FontWeight','bold');

ax2 = subplot('Position',[0.055 0.39 0.144 0.23]);
imagesc(abs(recos_meas(:,:,2)), c_lims2);
axis image; axis off;
yline(x_line,'LineWidth',1.5,'Color',blue);
xline(y_line,'LineWidth',1.5,'Color',blue);
text(-20,96, '16 interleaves','FontWeight','bold','HorizontalAlignment','center','Rotation',90,'FontName','Times','Fontsize',10);
i = i+1;
set(gca,'FontName','Times','Fontsize',9);
text(-70,10,'(D)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(200,10,'(E)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(740,10,'(F)','FontName','Arial','Fontsize',11,'FontWeight','bold');

ax3 = subplot('Position',[0.055 0.1 0.144 0.23]);
imagesc(abs(recos_meas(:,:,3)), c_lims3);
axis image; axis off;
yline(x_line,'LineWidth',1.5,'Color',blue);
xline(y_line,'LineWidth',1.5,'Color',blue);
text(-20,96, '3 interleaves','FontWeight','bold','HorizontalAlignment','center','Rotation',90,'FontName','Times','Fontsize',10);
i = i+1;
set(gca,'FontName','Times','Fontsize',9);
text(-70,10,'(G)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(200,10,'(H)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(740,10,'(I)','FontName','Arial','Fontsize',11,'FontWeight','bold');

ax4 = subplot('Position',[0.3 0.68 0.3 0.23]);
plot(squeeze(abs(recos_meas(x_line,:,1))), 'DisplayName','measured', 'Color',orange,'LineWidth',1.2);
hold on;
plot(squeeze(abs(recos_del(x_line,:,1))), 'DisplayName','delay', 'Color',yellow,'LineWidth',1.2);
plot(squeeze(abs(recos_girf(x_line,:,1))), 'DisplayName','GSTF', 'Color',lblue,'LineWidth',1.2);
plot(squeeze(abs(recos_girfdel(x_line,:,1))), 'DisplayName','GSTF+delay', 'Color',violet,'LineWidth',1.2);
ylabel('image intensity (a.u.)');
text(92,0.29, 'horizontal profiles','FontWeight','bold','HorizontalAlignment','center','Rotation',0,'FontName','Times','Fontsize',10);
i = i+1;
set(gca,'FontName','Times','Fontsize',9);

ax5 = subplot('Position',[0.3 0.39 0.3 0.23]);
plot(squeeze(abs(recos_meas(x_line,:,2))), 'DisplayName','measured', 'Color',orange,'LineWidth',1.2);
hold on;
plot(squeeze(abs(recos_del(x_line,:,2))), 'DisplayName','delay', 'Color',yellow,'LineWidth',1.2);
plot(squeeze(abs(recos_girf(x_line,:,2))), 'DisplayName','GSTF', 'Color',lblue,'LineWidth',1.2);
plot(squeeze(abs(recos_girfdel(x_line,:,2))), 'DisplayName','GSTF+delay', 'Color',violet,'LineWidth',1.2);
leg = legend('numColumns',4,'Position',[0.3 0.95 0.68 0.04]);
ylabel('image intensity (a.u.)');
i = i+1;
set(gca,'FontName','Times','Fontsize',9);

ax6 = subplot('Position',[0.3 0.1 0.3 0.23]);
plot(squeeze(abs(recos_meas(x_line,:,3))), 'DisplayName','measured', 'Color',orange,'LineWidth',1.2);
hold on;
plot(squeeze(abs(recos_del(x_line,:,3))), 'DisplayName','delay', 'Color',yellow,'LineWidth',1.2);
plot(squeeze(abs(recos_girf(x_line,:,3))), 'DisplayName','GSTF', 'Color',lblue,'LineWidth',1.2);
plot(squeeze(abs(recos_girfdel(x_line,:,3))), 'DisplayName','GSTF and delay', 'Color',violet,'LineWidth',1.2);
ylabel('image intensity (a.u.)');
xlabel('pixel number');
i = i+1;
set(gca,'FontName','Times','Fontsize',9);

ax7 = subplot('Position',[0.68 0.68 0.3 0.23]);
plot(squeeze(abs(recos_meas(:,y_line,1))), 'DisplayName','measured', 'Color',orange,'LineWidth',1.2);
hold on;
plot(squeeze(abs(recos_del(:,y_line,1))), 'DisplayName','delay', 'Color',yellow,'LineWidth',1.2);
plot(squeeze(abs(recos_girf(:,y_line,1))), 'DisplayName','GSTF', 'Color',lblue,'LineWidth',1.2);
plot(squeeze(abs(recos_girfdel(:,y_line,1))), 'DisplayName','GSTF and delay', 'Color',violet,'LineWidth',1.2);
text(96,0.29, 'vertical profiles','FontWeight','bold','HorizontalAlignment','center','Rotation',0,'FontName','Times','Fontsize',10);
i = i+1;
set(gca,'FontName','Times','Fontsize',9);

ax8 = subplot('Position',[0.68 0.39 0.3 0.23]);
plot(squeeze(abs(recos_meas(:,y_line,2))), 'DisplayName','measured', 'Color',orange,'LineWidth',1.2);
hold on;
plot(squeeze(abs(recos_del(:,y_line,2))), 'DisplayName','delay', 'Color',yellow,'LineWidth',1.2);
plot(squeeze(abs(recos_girf(:,y_line,2))), 'DisplayName','GSTF', 'Color',lblue,'LineWidth',1.2);
plot(squeeze(abs(recos_girfdel(:,y_line,2))), 'DisplayName','GSTF and delay', 'Color',violet,'LineWidth',1.2);
i = i+1;
set(gca,'FontName','Times','Fontsize',9);

ax9 = subplot('Position',[0.68 0.1 0.3 0.23]);
plot(squeeze(abs(recos_meas(:,y_line,3))), 'DisplayName','measured', 'Color',orange,'LineWidth',1.2);
hold on;
plot(squeeze(abs(recos_del(:,y_line,3))), 'DisplayName','delay', 'Color',yellow,'LineWidth',1.2);
plot(squeeze(abs(recos_girf(:,y_line,3))), 'DisplayName','GSTF', 'Color',lblue,'LineWidth',1.2);
plot(squeeze(abs(recos_girfdel(:,y_line,3))), 'DisplayName','GSTF and delay', 'Color',violet,'LineWidth',1.2);
xlabel('pixel number');
i = i+1;
set(gca,'FontName','Times','Fontsize',9);

linkaxes([ax4 ax5 ax6 ax7 ax8 ax9],'x');
xlim(ax4, [1,192]);

%% Plot gradients (Figure 5)
meas = meas2;
arm = 16;

figure('Units','centimeters','Position',[0 0 17.56 24]);

ax1 = subplot('Position',[0.11 0.86 0.86 0.1]);
plot(meas.t_ADC*1000, meas.grad_x_nom(:,arm),'LineWidth',1.2,'DisplayName','nominal');
hold on;
plot(meas.t_ADC*1000, meas.grad_x_del(:,arm),'LineWidth',1.2,'DisplayName','delay');
plot(meas.t_ADC*1000, meas.grad_x_girf(:,arm),'LineWidth',1.2,'DisplayName','GSTF');
plot(meas.t_ADC*1000, meas.grad_x_girfdel(:,arm),'LineWidth',1.2,'DisplayName','GSTF+delay');
plot(meas.t_ADC*1000, meas.grad_x_meas(:,arm),'LineWidth',1.2,'DisplayName','measured');
ylabel('Gradient_x (T/m)');
legend('numColumns',5,'Position',[0.11 0.97 0.86 0.02]);
rectangle('Position',[0.02 -0.145 0.96 0.29],'LineStyle','--','LineWidth',1.2,'EdgeColor','k');
set(gca,'FontName','Times','Fontsize',9);
text(-1.25,0.2,'(A)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(-1.25,-0.3,'(B)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(-1.25,-0.9,'(C)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(-1.25,-1.6,'(D)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(-1.25,-2.25,'(E)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(-1.25,-2.95,'(F)','FontName','Arial','Fontsize',11,'FontWeight','bold');

ax2 = subplot('Position',[0.11 0.73 0.86 0.1]);
plot(meas.t_ADC*1000, meas.grad_x_meas(:,arm)-meas.grad_x_nom(:,arm),'LineWidth',1.2);
hold on;
plot(meas.t_ADC*1000, meas.grad_x_meas(:,arm)-meas.grad_x_del(:,arm),'LineWidth',1.2);
plot(meas.t_ADC*1000, meas.grad_x_meas(:,arm)-meas.grad_x_girf(:,arm),'LineWidth',1.2);
plot(meas.t_ADC*1000, meas.grad_x_meas(:,arm)-meas.grad_x_girfdel(:,arm),'LineWidth',1.2);
rectangle('Position',[0.02 -0.03 0.96 0.06],'LineStyle','-.','LineWidth',1.2,'EdgeColor','k');
ylabel('\Delta Gradient_x (T/m)');
ylim([-0.08 0.08]);
set(gca,'FontName','Times','Fontsize',9);

ax3 = subplot('Position',[0.11 0.54 0.86 0.14]);
plot(meas.t_ADC*1000, meas.grad_x_nom(:,arm),'LineWidth',1.2);
hold on;
plot(meas.t_ADC*1000, meas.grad_x_del(:,arm),'LineWidth',1.2);
plot(meas.t_ADC*1000, meas.grad_x_girf(:,arm),'LineWidth',1.2);
plot(meas.t_ADC*1000, meas.grad_x_girfdel(:,arm),'LineWidth',1.2);
plot(meas.t_ADC*1000, meas.grad_x_meas(:,arm),'LineWidth',1.2);
yline(0);
xline(0.38,'--');
xline(0.7,'--');
ylabel('Gradient_x (T/m)');
ylim([-0.15 0.15]);
rectangle('Position',[0.002 -0.145 0.996 0.29],'LineStyle','--','LineWidth',1.2,'EdgeColor','k');
annotation('arrow',0.775+[0 -0.02],0.605+[0 0.02],'Color',green,'LineWidth',1.5);
annotation('arrow',0.505+[0 -0.02],0.62+[0 -0.02],'Color',green,'LineWidth',1.5);
set(gca,'FontName','Times','Fontsize',9);
box off;

ax3_inset1 = axes('Position',[0.55 0.605 0.1 0.07]);
plot(meas.t_ADC*1000, meas.grad_x_nom(:,arm),'LineWidth',1.2);
hold on;
plot(meas.t_ADC*1000, meas.grad_x_del(:,arm),'LineWidth',1.2);
plot(meas.t_ADC*1000, meas.grad_x_girf(:,arm),'LineWidth',1.2);
plot(meas.t_ADC*1000, meas.grad_x_girfdel(:,arm),'LineWidth',1.2);
plot(meas.t_ADC*1000, meas.grad_x_meas(:,arm),'LineWidth',1.4);
xlim([0.36 0.46]);
ylim([-0.06 0.01]);

ax3_inset2 = axes('Position',[0.82 0.565 0.1 0.07]);
plot(meas.t_ADC*1000, meas.grad_x_nom(:,arm),'LineWidth',1.2);
hold on;
plot(meas.t_ADC*1000, meas.grad_x_del(:,arm),'LineWidth',1.2);
plot(meas.t_ADC*1000, meas.grad_x_girf(:,arm),'LineWidth',1.2);
plot(meas.t_ADC*1000, meas.grad_x_girfdel(:,arm),'LineWidth',1.2);
plot(meas.t_ADC*1000, meas.grad_x_meas(:,arm),'LineWidth',1.4);
xlim([0.69 0.79]);
ylim([-0.02 0.07]);

% ax4 = subplot(4,1,4);
ax4 = subplot('Position',[0.11 0.39 0.86 0.12]);
plot(meas.t_ADC*1000, meas.grad_x_meas(:,arm)-meas.grad_x_nom(:,arm),'LineWidth',1.2);
hold on;
plot(meas.t_ADC*1000, meas.grad_x_meas(:,arm)-meas.grad_x_del(:,arm),'LineWidth',1.2);
plot(meas.t_ADC*1000, meas.grad_x_meas(:,arm)-meas.grad_x_girf(:,arm),'LineWidth',1.2);
plot(meas.t_ADC*1000, meas.grad_x_meas(:,arm)-meas.grad_x_girfdel(:,arm),'LineWidth',1.2);
yline(0);
xline(0.38,'--');
xline(0.7,'--');
ylabel('\Delta Gradient_x (T/m)');
ylim([-0.032 0.032]);
rectangle('Position',[0.002 -0.031 0.996 0.062],'LineStyle','-.','LineWidth',1.2,'EdgeColor','k');
annotation('arrow',0.545+[0 -0.02],0.475+[0 -0.02],'Color',yellow,'LineWidth',1.5);
annotation('arrow',0.638+[0 0.02],0.423+[0 0.02],'Color',yellow,'LineWidth',1.5);
annotation('arrow',0.555+[0 -0.02],0.417+[0 0.02],'Color',violet,'LineWidth',1.5);
annotation('arrow',0.625+[0 0.02],0.485+[0 -0.02],'Color',violet,'LineWidth',1.5);
set(gca,'FontName','Times','Fontsize',9);
box off;

ax5 = subplot('Position',[0.11 0.2 0.86 0.14]);
plot(meas.t_ADC*1000, meas.grad_y_nom(:,arm),'LineWidth',1.2,'DisplayName','nominal','Color',blue);
hold on;
plot(meas.t_ADC*1000, meas.grad_y_del(:,arm),'LineWidth',1.2,'DisplayName','delay','Color',orange);
plot(meas.t_ADC*1000, meas.grad_y_girf(:,arm),'LineWidth',1.2,'DisplayName','GSTF','Color',yellow);
plot(meas.t_ADC*1000, meas.grad_y_girfdel(:,arm),'LineWidth',1.2,'DisplayName','GSTF with delay','Color',violet);
plot(meas.t_ADC*1000, meas.grad_y_meas(:,arm),'LineWidth',1.5,'DisplayName','measured','Color',green);
ylabel('Gradient_y (T/m)');
yline(0);
xline(0.26,'--');
xline(0.54,'--');
xline(0.9,'--');
ylim([-0.15 0.15]);
annotation('arrow',0.41+[0 -0.02],0.278+[0 -0.02],'Color',green,'LineWidth',1.5);
annotation('arrow',0.652+[0 -0.02],0.265+[0 0.02],'Color',green,'LineWidth',1.5);
annotation('arrow',0.955+[0 -0.02],0.278+[0 -0.02],'Color',green,'LineWidth',1.5);
set(gca,'FontName','Times','Fontsize',9);

ax6 = subplot('Position',[0.11 0.05 0.86 0.12]);
plot(meas.t_ADC*1000, meas.grad_y_meas(:,arm)-meas.grad_y_nom(:,arm),'LineWidth',1.2,'DisplayName','nominal','Color',blue);
hold on;
plot(meas.t_ADC*1000, meas.grad_y_meas(:,arm)-meas.grad_y_del(:,arm),'LineWidth',1.2,'DisplayName','delay','Color',orange);
plot(meas.t_ADC*1000, meas.grad_y_meas(:,arm)-meas.grad_y_girf(:,arm),'LineWidth',1.2,'DisplayName','GSTF','Color',yellow);
plot(meas.t_ADC*1000, meas.grad_y_meas(:,arm)-meas.grad_y_girfdel(:,arm),'LineWidth',1.2,'DisplayName','GSTF with delay','Color',violet);
yline(0);
xline(0.26,'--');
xline(0.54,'--');
xline(0.9,'--');
ylabel('\Delta Gradient_y (T/m)');
xlabel('Time (ms)');
ylim([-0.032 0.032]);
set(gca,'FontName','Times','Fontsize',9);

linkaxes([ax1 ax2],'x');
xlim(ax1, [0 10.2]);
linkaxes([ax3 ax4 ax5 ax6],'x');
xlim(ax3,[0 1]);

%% Plot trajectories (Figure 6)
figure('Units','centimeters','Position',[0 0 17.56 13]);

ax1 = subplot('Position',[0.07 0.31 0.4 0.7]);
plot(meas.traj_nom(:,arm),'LineWidth',1.2,'DisplayName','nominal');
hold on;
plot(meas.traj_del(:,arm),'LineWidth',1.2,'DisplayName','delay');
plot(meas.traj_girf(:,arm),'LineWidth',1.2,'DisplayName','GSTF');
plot(meas.traj_girfdel(:,arm),'LineWidth',1.2,'DisplayName','GSTF+delay');
plot(meas.traj_meas(:,arm),'LineWidth',1.2,'DisplayName','measured');
axis image;
xlabel('k_x (1/m)');
ylabel('k_y (1/m)');
xlim([-3.4 3.4]);
ylim([-3.4 3.4]);
legend('numColumns',5,'Position',[0.07 0.95 0.91 0.04]);
rectangle('Position',[-1 1.5 1 1],'LineStyle','--');
rectangle('Position',[0 -2.7 1 1],'LineStyle','-.');
set(gca,'FontName','Times','Fontsize',9);
text(-4.4,3.3,'(A)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(3.9,3.3,'(B)','FontName','Arial','Fontsize',11,'FontWeight','bold');

ax1a = axes('Position',[0.04 0.05 0.21 0.25]);
plot(meas.traj_nom(:,arm),'LineWidth',1.2);
hold on;
plot(meas.traj_del(:,arm),'LineWidth',1.2);
plot(meas.traj_girf(:,arm),'LineWidth',1.2);
plot(meas.traj_girfdel(:,arm),'LineWidth',1.2);
plot(meas.traj_meas(:,arm),'LineWidth',1.2);
axis image; box off;
xlim([-1 0]);
ylim([1.5 2.5]);
rectangle('Position',[-1 1.5 1 1],'LineStyle','--','LineWidth',1.2);
set(gca,'FontName','Times','Fontsize',9);
text(-0.97,2.4,'(C)','FontName','Arial','Fontsize',11,'FontWeight','bold');

ax1b = axes('Position',[0.28 0.05 0.21 0.25]);
plot(meas.traj_nom(:,arm),'LineWidth',1.2);
hold on;
plot(meas.traj_del(:,arm),'LineWidth',1.2);
plot(meas.traj_girf(:,arm),'LineWidth',1.2);
plot(meas.traj_girfdel(:,arm),'LineWidth',1.2);
plot(meas.traj_meas(:,arm),'LineWidth',1.2);
axis image; box off;
xlim([0 1]);
ylim([-2.7 -1.7]);
rectangle('Position',[0 -2.7 1 1],'LineStyle','-.','LineWidth',1.2);
set(gca,'FontName','Times','Fontsize',9);
text(0.03,-1.8,'(D)','FontName','Arial','Fontsize',11,'FontWeight','bold');


ax2 = subplot('Position',[0.58 0.31 0.4 0.7]);
plot(meas.traj_meas(:,arm)-meas.traj_nom(:,arm),'LineWidth',1.2);
hold on;
plot(meas.traj_meas(:,arm)-meas.traj_del(:,arm),'LineWidth',1.2);
plot(meas.traj_meas(:,arm)-meas.traj_girf(:,arm),'LineWidth',1.2);
plot(meas.traj_meas(:,arm)-meas.traj_girfdel(:,arm),'LineWidth',1.2);
axis image;
xlabel('\Delta k_x (1/m)');
ylabel('\Delta k_y (1/m)');
xlim([-0.18 0.18]);
ylim([-0.18 0.18]);
xline(0,'--');
yline(0,'--');
set(gca,'FontName','Times','Fontsize',9);

ax2a = axes('Position',[0.55 0.05 0.21 0.25]);
idx_end = 20;
plot(meas.traj_meas(1:idx_end,arm)-meas.traj_nom(1:idx_end,arm),'LineWidth',1.2);
hold on;
plot(meas.traj_meas(1:idx_end,arm)-meas.traj_del(1:idx_end,arm),'LineWidth',1.2);
plot(meas.traj_meas(1:idx_end,arm)-meas.traj_girf(1:idx_end,arm),'LineWidth',1.2);
plot(meas.traj_meas(1:idx_end,arm)-meas.traj_girfdel(1:idx_end,arm),'LineWidth',1.2);
axis image;
xlim([-0.008 0.008]);
ylim([-0.008 0.008]);
xline(0,'--');
yline(0,'--');
ax2a.XAxis.Exponent = 0;
ax2a.YAxis.Exponent = 0;
text(0,-0.0069,'  0 - 0.1 ms','FontName','Times','Fontsize',9);
set(gca,'FontName','Times','Fontsize',9);
text(-0.0072,0.0064,'(E)','FontName','Arial','Fontsize',11,'FontWeight','bold');

ax2b = axes('Position',[0.79 0.05 0.21 0.25]);
idx_end = 200;
plot(meas.traj_meas(1:idx_end,arm)-meas.traj_nom(1:idx_end,arm),'LineWidth',1.2);
hold on;
plot(meas.traj_meas(1:idx_end,arm)-meas.traj_del(1:idx_end,arm),'LineWidth',1.2);
plot(meas.traj_meas(1:idx_end,arm)-meas.traj_girf(1:idx_end,arm),'LineWidth',1.2);
plot(meas.traj_meas(1:idx_end,arm)-meas.traj_girfdel(1:idx_end,arm),'LineWidth',1.2);
axis image;
xlim([-0.06 0.06]);
ylim([-0.06 0.06]);
xline(0,'--');
yline(0,'--');
text(0,-0.052,'  0 - 1 ms','FontName','Times','Fontsize',9);
set(gca,'FontName','Times','Fontsize',9);
text(-0.055,0.048,'(F)','FontName','Arial','Fontsize',11,'FontWeight','bold');

rmse_nom = sqrt(mean(abs(meas.traj_meas(:,arm)-meas.traj_nom(:,arm)).^2));
disp(['rmse_nom = ',num2str(rmse_nom)]);
rmse_del = sqrt(mean(abs(meas.traj_meas(:,arm)-meas.traj_del(:,arm)).^2));
disp(['rmse_del = ',num2str(rmse_del)]);
rmse_girf = sqrt(mean(abs(meas.traj_meas(:,arm)-meas.traj_girf(:,arm)).^2));
disp(['rmse_girf = ',num2str(rmse_girf)]);
rmse_girfdel = sqrt(mean(abs(meas.traj_meas(:,arm)-meas.traj_girfdel(:,arm)).^2));
disp(['rmse_girfdel = ',num2str(rmse_girfdel)]);
disp(' ');
rmse_nom = sqrt(mean(abs(meas.traj_meas(1:20,arm)-meas.traj_nom(1:20,arm)).^2));
disp(['rmse_nom 0.1ms = ',num2str(rmse_nom)]);
rmse_del = sqrt(mean(abs(meas.traj_meas(1:20,arm)-meas.traj_del(1:20,arm)).^2));
disp(['rmse_del 0.1ms = ',num2str(rmse_del)]);
rmse_girf = sqrt(mean(abs(meas.traj_meas(1:20,arm)-meas.traj_girf(1:20,arm)).^2));
disp(['rmse_girf 0.1ms = ',num2str(rmse_girf)]);
rmse_girfdel = sqrt(mean(abs(meas.traj_meas(1:20,arm)-meas.traj_girfdel(1:20,arm)).^2));
disp(['rmse_girfdel 0.1ms = ',num2str(rmse_girfdel)]);
disp(' ');
rmse_nom = sqrt(mean(abs(meas.traj_meas(1:200,arm)-meas.traj_nom(1:200,arm)).^2));
disp(['rmse_nom 1ms = ',num2str(rmse_nom)]);
rmse_del = sqrt(mean(abs(meas.traj_meas(1:200,arm)-meas.traj_del(1:200,arm)).^2));
disp(['rmse_del 1ms = ',num2str(rmse_del)]);
rmse_girf = sqrt(mean(abs(meas.traj_meas(1:200,arm)-meas.traj_girf(1:200,arm)).^2));
disp(['rmse_girf 1ms = ',num2str(rmse_girf)]);
rmse_girfdel = sqrt(mean(abs(meas.traj_meas(1:200,arm)-meas.traj_girfdel(1:200,arm)).^2));
disp(['rmse_girfdel 1ms = ',num2str(rmse_girfdel)]);










