% Copyright (c) 2024 Hannah Scholten

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
path2save = "..\invivo_data\RMSEs_and_delays_m96.mat";

% initialize arrays for results
delays_girf = ((-10:0))*1e-6;
delays_del = ((-10:0)-15)*1e-6;
sum_abs_err_girf = zeros(size(delays_girf));
sum_abs_err_del = zeros(size(delays_del));
sum_gradx_err_girf = zeros(size(delays_girf));
sum_grady_err_girf = zeros(size(delays_girf));
sum_gradz_err_girf = zeros(size(delays_girf));
sum_gradx_err_del = zeros(size(delays_del));
sum_grady_err_del = zeros(size(delays_del));
sum_gradz_err_del = zeros(size(delays_del));
sum_traj_err_girf = zeros(size(delays_girf));
sum_traj_err_del = zeros(size(delays_del));
recos_girf = zeros(matrix_size, matrix_size, length(delays_girf));
recos_del = zeros(matrix_size, matrix_size, length(delays_del));

%% load GSTFs
H_x = load('..\GSTF_data\Hx_sr100_20231009.mat');
H_y = load('..\GSTF_data\Hy_sr100_20231009.mat');
H_z = load('..\GSTF_data\Hz_sr100_20231009.mat');

%% Reconstruct images with different delays
for loop = 1:length(delays_girf)
    disp(['loop = ',num2str(loop)]);

    delay_x = delays_girf(loop); % for GSTF+delay correction
    delay_x2 = delays_del(loop); % for isotropic delay correction
    delay_y = delay_x;
    delay_z = delay_x;
    H_1 = H_x.H_matrix;
    H_2 = H_y.H_matrix;
    H_3 = H_z.H_matrix;
    % Add delays to GSTFs
    gstf_x = H_x.H_matrix.gstf(:,2) .* exp(1i*2*pi*H_x.H_matrix.f_axis.' *delay_x);
    gstf_y = H_y.H_matrix.gstf(:,2) .* exp(1i*2*pi*H_y.H_matrix.f_axis.' *delay_y);
    gstf_z = H_z.H_matrix.gstf(:,2) .* exp(1i*2*pi*H_z.H_matrix.f_axis.' *delay_z);

    delay_y2 = delay_x2;
    delay_z2 = delay_x2;

    dwelltime_girf = H_x.H_matrix.dt;

    %% Load data
    disp('Read measurement data.')

    if isfile(filepath)
        load(filepath);
        spiral_data_traj = load(filepath_traj);
        spiral_data_traj = spiral_data_traj.spiral_data;
    else
        disp(['Error: cannot find ',filepath]);
        break
    end
    raw = spiral_data.raw;
    nPost = spiral_data.nPost; % data points discarded for reconstruction

    %% reconstruct images with measured trajectory
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
    traj_meas = traj_inPlaneX + 1i*traj_inPlaneY;

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
    t_GRT = ((1:1:(size(spiral_x,1)))-0.5)*grad_raster_time; % + PVM_SpiralPreSize, aber ist bei uns immer(?) 0
    t_eps = 1e-12;
    t_GRT_00 = [-t_eps, t_GRT(1)-t_eps, t_GRT, t_GRT(end)+t_eps, t_GRT(end)+grad_raster_time+t_eps];
    
    spiral_x_00 = [zeros(2,size(spiral_x,2)); spiral_x; zeros(2,size(spiral_x,2))];
    spiral_y_00 = [zeros(2,size(spiral_x,2)); spiral_y; zeros(2,size(spiral_x,2))];

    spiral_x_pp = interp1(t_GRT_00, spiral_x_00*max_grad_strength, 'linear', 'pp');
    spiral_y_pp = interp1(t_GRT_00, spiral_y_00*max_grad_strength, 'linear', 'pp');

    traj_x_pp = fnint(spiral_x_pp);
    traj_y_pp = fnint(spiral_y_pp);

    %% reconstruct with nominal trajectory
    gamma = 267.513*10^6; %Hz/T
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
        gradXYZ = spiral_data.gradMatrix * grad_xyz;
        grad_physicalX(:,arm) = gradXYZ(1,:);
        grad_physicalY(:,arm) = gradXYZ(2,:);
        grad_physicalZ(:,arm) = gradXYZ(3,:);
    end

    %% zerofill gradients to avoid side effects by the fft calculation
    nExtra = round((1e6-size(grad_x_tgirf,1))/2);
    grad_x_ZF = cat(1, zeros(nExtra,size(spiral_x,2)), grad_physicalX, zeros(nExtra,size(spiral_x,2)));
    grad_y_ZF = cat(1, zeros(nExtra,size(spiral_y,2)), grad_physicalY, zeros(nExtra,size(spiral_x,2)));
    grad_z_ZF = cat(1, zeros(nExtra,size(spiral_x,2)), grad_physicalZ, zeros(nExtra,size(spiral_x,2)));

    gstf_x_interp = interp1(H_1.f_axis, gstf_x, linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).', 'linear');
    gstf_y_interp = interp1(H_2.f_axis, gstf_y, linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).', 'linear');
    gstf_z_interp = interp1(H_3.f_axis, gstf_z, linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).', 'linear');

    %% calculate post-correction with GSTF+delay
    tmp1 = fft1d(grad_x_ZF,1);
    tmp2 = repmat(gstf_x_interp, [1,spiral_data.numInterleaves]);
    grad_px_girf = real(ifft1d( tmp1.*tmp2, 1));% + H_x.H_matrix.fieldOffsets(2);
    grad_px_girf = grad_px_girf(nExtra+1:end-nExtra,:);

    tmp1 = fft1d(grad_y_ZF,1);
    tmp2 = repmat(gstf_y_interp, [1,spiral_data.numInterleaves]);
    grad_py_girf = real(ifft1d( tmp1.*tmp2, 1));% + H_y.H_matrix.fieldOffsets(2);
    grad_py_girf = grad_py_girf(nExtra+1:end-nExtra,:);

    tmp1 = fft1d(grad_z_ZF,1);
    tmp2 = repmat(gstf_z_interp, [1,spiral_data.numInterleaves]);
    grad_pz_girf = real(ifft1d( tmp1.*tmp2, 1));% + H_z.H_matrix.fieldOffsets(2);
    grad_pz_girf = grad_pz_girf(nExtra+1:end-nExtra,:);

    % rotate back to PRS (xyz/ in-plane) coordinate system
    grad_x_girf = zeros(size(grad_x_tgirf));
    grad_y_girf = zeros(size(grad_x_tgirf));

    for arm = 1:spiral_data.numInterleaves
        gradXYZ = [squeeze(grad_px_girf(:,arm)), squeeze(grad_py_girf(:,arm)), squeeze(grad_pz_girf(:,arm))].';
        grad_xyz = spiral_data.gradMatrix.' * gradXYZ;
        grad_x_girf(:,arm) = grad_xyz(1,:);
        grad_y_girf(:,arm) = grad_xyz(2,:);
    end

    % interpolate back to ADC time grid for reconstruction
    grad_x_girf_pp = interp1(t_GIRF, grad_x_girf, 'linear','pp');
    grad_y_girf_pp = interp1(t_GIRF, grad_y_girf, 'linear','pp');
    traj_x_girf_pp = fnint(grad_x_girf_pp);
    traj_y_girf_pp = fnint(grad_y_girf_pp);

    clearvars grad_x_girf grad_y_girf gradXYZ tmp1 tmp2 

    %% reconstruct with GSTF+delay-corrected trajectory
    trajx_corr = gamma/(2*pi)*ppval(traj_x_girf_pp, t_ADC)/1000;
    trajy_corr = gamma/(2*pi)*ppval(traj_y_girf_pp, t_ADC)/1000;
    traj_corr = trajx_corr.' + 1i*trajy_corr.';

    delta_traj_girf = (traj_corr(1:end-nPost,:) - traj_meas(1:end-nPost,:))/max(abs(traj_nom),[],'all');
    sum_traj_err_girf(loop) = sqrt(sum(abs(delta_traj_girf).^2,'all')/size(delta_traj_girf,1)/size(delta_traj_girf,2));

    FT = cGrid(traj_corr(1:end-nPost,:)/traj_factor,matrix_size);
    reco_girf = FT'*raw(1:end-nPost,:);

    %% calculate post-correction with isotropic delay
    gstf_x_del = exp(1i*2*pi*linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).' *delay_x2);
    gstf_y_del = exp(1i*2*pi*linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).' *delay_y2);
    gstf_z_del = exp(1i*2*pi*linspace(H_1.f_axis(1),H_1.f_axis(end),size(grad_x_ZF,1)).' *delay_z2);

    tmp1 = fft1d(grad_x_ZF,1);
    tmp2 = repmat(gstf_x_del, [1,spiral_data.numInterleaves]);
    grad_px_del = real(ifft1d( tmp1.*tmp2, 1));
    grad_px_del = grad_px_del(nExtra+1:end-nExtra,:);

    tmp1 = fft1d(grad_y_ZF,1);
    tmp2 = repmat(gstf_y_del, [1,spiral_data.numInterleaves]);
    grad_py_del = real(ifft1d( tmp1.*tmp2, 1));
    grad_py_del = grad_py_del(nExtra+1:end-nExtra,:);

    tmp1 = fft1d(grad_z_ZF,1);
    tmp2 = repmat(gstf_z_del, [1,spiral_data.numInterleaves]);
    grad_pz_del = real(ifft1d( tmp1.*tmp2, 1));
    grad_pz_del = grad_pz_del(nExtra+1:end-nExtra,:);

    % rotate back to PRS (xyz/ in-plane) coordinate system
    grad_x_del = zeros(size(grad_x_tgirf));
    grad_y_del = zeros(size(grad_x_tgirf));

    for arm = 1:spiral_data.numInterleaves
        gradXYZ = [squeeze(grad_px_del(:,arm)), squeeze(grad_py_del(:,arm)), squeeze(grad_pz_del(:,arm))].';
        grad_xyz = spiral_data.gradMatrix.' * gradXYZ;
        grad_x_del(:,arm) = grad_xyz(1,:);
        grad_y_del(:,arm) = grad_xyz(2,:);
    end

    % interpolate back to ADC time grid for reconstruction
    grad_x_del_pp = interp1(t_GIRF, grad_x_del, 'linear','pp');
    grad_y_del_pp = interp1(t_GIRF, grad_y_del, 'linear','pp');
    traj_x_del_pp = fnint(grad_x_del_pp);
    traj_y_del_pp = fnint(grad_y_del_pp);

    %% reconstruct with delay-corrected trajectory
    trajx_del = gamma/(2*pi)*ppval(traj_x_del_pp, t_ADC)/1000;
    trajy_del = gamma/(2*pi)*ppval(traj_y_del_pp, t_ADC)/1000;
    traj_del = trajx_del.' + 1i*trajy_del.';

    delta_traj_del = (traj_del(1:end-nPost,:) - traj_meas(1:end-nPost,:))/max(abs(traj_nom),[],'all');
    sum_traj_err_del(loop) = sqrt(sum(abs(delta_traj_del).^2,'all')/size(delta_traj_del,1)/size(delta_traj_del,2));

    FT = cGrid(traj_del(1:end-nPost,:)/traj_factor,matrix_size);
    reco_delay = FT'*raw(1:end-nPost,:);

    %% Compare measured and predicted gradient waveforms
    % Calculate measured gradient waveforms as derivative of measured trajectory
    trajx_meas_pp = interp1(t_ADC, trajx_meas, 'linear','pp');
    trajy_meas_pp = interp1(t_ADC, trajy_meas, 'linear','pp');
    trajz_meas_pp = interp1(t_ADC, trajz_meas, 'linear','pp');
    gradx_meas_pp = fnder(trajx_meas_pp);
    grady_meas_pp = fnder(trajy_meas_pp);
    gradz_meas_pp = fnder(trajz_meas_pp);
    % ... on ADC time grid
    grad_x_meas = ppval(gradx_meas_pp, t_ADC).'/gamma*2*pi*1000;
    grad_y_meas = ppval(grady_meas_pp, t_ADC).'/gamma*2*pi*1000;
    grad_z_meas = ppval(gradz_meas_pp, t_ADC).'/gamma*2*pi*1000;
    % ... on GIRF time grid
    grad_x_meas_tGIRF = ppval(gradx_meas_pp, t_GIRF).'/gamma*2*pi*1000;
    grad_y_meas_tGIRF = ppval(grady_meas_pp, t_GIRF).'/gamma*2*pi*1000;
    grad_z_meas_tGIRF = ppval(gradz_meas_pp, t_GIRF).'/gamma*2*pi*1000;

    % Calculate RMSEs between the gradient waveforms on each axis
    delta_Gx_girf = (grad_px_del - grad_x_meas_tGIRF)./max(grad_physicalX,[],1);
    sum_gradx_err_girf(loop) = sqrt(sum(abs(delta_Gx_girf).^2,'all')/size(delta_Gx_girf,1)/size(delta_Gx_girf,2));
    delta_Gy_girf = (grad_py_del - grad_y_meas_tGIRF)./max(grad_physicalY,[],1);
    sum_grady_err_girf(loop) = sqrt(sum(abs(delta_Gy_girf).^2,'all')/size(delta_Gy_girf,1)/size(delta_Gy_girf,2));
    delta_Gz_girf = (grad_pz_del - grad_z_meas_tGIRF)./max(grad_physicalZ,[],1);
    sum_gradz_err_girf(loop) = sqrt(sum(abs(delta_Gz_girf).^2,'all')/size(delta_Gz_girf,1)/size(delta_Gz_girf,2));

    delta_Gx_del = (grad_px_del - grad_x_meas_tGIRF)./max(grad_physicalX,[],1);
    sum_gradx_err_del(loop) = sqrt(sum(abs(delta_Gx_del).^2,'all')/size(delta_Gx_del,1)/size(delta_Gx_del,2));
    delta_Gy_del = (grad_py_del - grad_y_meas_tGIRF)./max(grad_physicalY,[],1);
    sum_grady_err_del(loop) = sqrt(sum(abs(delta_Gy_del).^2,'all')/size(delta_Gy_del,1)/size(delta_Gy_del,2));
    delta_Gz_del = (grad_pz_del - grad_z_meas_tGIRF)./max(grad_physicalZ,[],1);
    sum_gradz_err_del(loop) = sqrt(sum(abs(delta_Gz_del).^2,'all')/size(delta_Gz_del,1)/size(delta_Gz_del,2));
    
    %% Calculate nRMSEs between corrected and ground truth images
    recos_girf(:,:,loop) = reco_girf;
    recos_del(:,:,loop) = reco_delay;
    
    s = sum(abs(reco_meas).*abs(reco_delay), 'all') / sum(abs(reco_delay).*abs(reco_delay), 'all');
    sum_abs_err_del(loop) = sqrt(sum((abs(reco_meas)-s*abs(reco_delay)).^2, 'all')/matrix_size^2);
    s = sum(abs(reco_meas).*abs(reco_girf), 'all') / sum(abs(reco_girf).*abs(reco_girf), 'all');
    sum_abs_err_girf(loop) = sqrt(sum((abs(reco_meas)-s*abs(reco_girf)).^2, 'all')/matrix_size^2);

end % loop
save(path2save, 'delays_girf','delays_del','sum_abs_err_del','sum_abs_err_girf', ...
    'sum_gradx_err_girf','sum_grady_err_girf','sum_gradz_err_girf','sum_gradx_err_del','sum_grady_err_del','sum_gradz_err_del', ...
    'sum_traj_err_girf','sum_traj_err_del');

%% Fit polynomial to error-vs-delay
% Initial guess for the parameters
initial_guess = [1, 1, 1, 1, 1];
% polynomial of 4th order
quad_function = @(params, x) params(1)*x.^4 + params(2)*x.^3 + params(3)*x.^2 + params(4)*x + params(5);
% Fit the data using lsqcurvefit
params_del = lsqcurvefit(quad_function, initial_guess, delays_del(1:11)*1e5, sum_abs_err_del(1:11)*100);
params_girf = lsqcurvefit(quad_function, initial_guess, delays_girf(1:11)*1e5, sum_abs_err_girf(1:11)*100);

% Find extrema
extrema_del = roots([4*params_del(1) 3*params_del(2) 2*params_del(3) params_del(4)]);
for i=1:3
    if imag(extrema_del(i))==0
        if 12*params_del(1)*extrema_del(i)^2 + 6*params_del(2)*extrema_del(i) + 2*params_del(3) > 0
            disp(['Min del = ',num2str(extrema_del(i)*1e-5)]);
            min_del = extrema_del(i)*1e-5;
        end
    end
end
extrema_girf = roots([4*params_girf(1) 3*params_girf(2) 2*params_girf(3) params_girf(4)]);
for i=1:3
    if imag(extrema_girf(i))==0
        if 12*params_girf(1)*extrema_girf(i)^2 + 6*params_girf(2)*extrema_girf(i) + 2*params_girf(3) > 0
            disp(['Min girf = ',num2str(extrema_girf(i)*1e-5)]);
            min_girf = extrema_girf(i)*1e-5;
        end
    end
end

delays_girf_fine = ((-100:0))*1e-7*1e5;
delays_del_fine = ((-100:0)-150)*1e-7*1e5;

fit_del = params_del(1)*delays_del_fine.^4 + params_del(2)*delays_del_fine.^3 + params_del(3)*delays_del_fine.^2 + params_del(4)*delays_del_fine + params_del(5);
fit_girf = params_girf(1)*delays_girf_fine.^4 + params_girf(2)*delays_girf_fine.^3 + params_girf(3)*delays_girf_fine.^2 + params_girf(4)*delays_girf_fine + params_girf(5);

fit_del = fit_del/100;
fit_girf = fit_girf/100;

%% Plot different RMSEs
blue = [0 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];
violet = [0.4940 0.1840 0.5560];
green = [0.4660 0.6740 0.1880];
yellow = [0.9290 0.6940 0.1250];
lblue = [0.3010 0.7450 0.9330];
dred = [0.6350 0.0780 0.1840];

figure;
subplot(2,1,1);
plot(delays_del, sum_abs_err_del, '-o', 'DisplayName','nRMSE(magnitude images)');
hold on;
plot(delays_del_fine*1e-5, fit_del, '--', 'DisplayName','Fit');
plot(delays_del, sum_gradx_err_del, '*', 'DisplayName','RMSE(G_x)');
plot(delays_del, sum_grady_err_del, '*', 'DisplayName','RMSE(G_y)');
plot(delays_del, sum_gradz_err_del, '*', 'DisplayName','RMSE(G_z)');
plot(delays_del, sum_traj_err_del, '+', 'DisplayName','RMSE(k)');
xline(min_del, 'DisplayName','min(nRMSE)');
title('(n)RMSEs with delay-corrected trajectory');
xlabel('Delay (s)');
ylabel('(n)RMSE (a.u.)');
legend('Location','northoutside','NumColumns',7);

subplot(2,1,2);
plot(delays_girf, sum_abs_err_girf, '-o');
hold on;
plot(delays_girf_fine*1e-5, fit_girf, '--');
plot(delays_girf, sum_gradx_err_girf, '*');
plot(delays_girf, sum_grady_err_girf, '*');
plot(delays_girf, sum_gradz_err_girf, '*');
plot(delays_girf, sum_traj_err_girf, '+');
xline(min_girf);
title('(n)RMSEs with GSTF+delay correction');
xlabel('Delay (s)');
ylabel('(n)RMSE (a.u.)');













