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
addpath(folder1);

%% Load GSTFs of x-, y-, and z-axis
H_x_sr100_data = load('..\GSTF_data\Hx_sr100_20231009.mat');
H_y_sr100_data = load('..\GSTF_data\Hy_sr100_20231009.mat');
H_z_sr100_data = load('..\GSTF_data\Hz_sr100_20231009.mat');

Hx_sr100 = H_x_sr100_data.H_matrix;
Hy_sr100 = H_y_sr100_data.H_matrix;
Hz_sr100 = H_z_sr100_data.H_matrix;

%% Define paths to spiral trajectory data
meas1 = struct;
meas1.filepath = '..\phantom_data\spiral_data_p96.mat';

meas2 = struct;
meas2.filepath = '..\phantom_data\spiral_data_p16.mat';

meas3 = struct;
meas3.filepath = '..\phantom_data\spiral_data_p03.mat';

measurements = {meas1, meas2, meas3};

%% Forward calculation of triangles

% Load prescribed gradient shapes
triang_x_sr100_data = load('..\GSTF_data\girf_data_x.mat').girf_data;
triang_x_sr100 = triang_x_sr100_data.gradShapeGIRF.';

% Define time grids
grad_raster_time = 8.1e-6;
dwelltime_meas = 1/double(triang_x_sr100_data.bandwidth);
t_ADC = (0:1:triang_x_sr100_data.datapoints-1)*dwelltime_meas + dwelltime_meas/2;
t_GRT = ((1:1:(size(triang_x_sr100,1)))-0.5)*grad_raster_time +triang_x_sr100_data.nPre*1e-6;
t_eps = 1e-12;
t_GRT_00 = [-t_eps, t_GRT(1)-t_eps, t_GRT, t_GRT(end)+t_eps, t_GRT(end)+grad_raster_time+t_eps];

triang_x_sr100_00 = [zeros(2,size(triang_x_sr100,2)); triang_x_sr100; zeros(2,size(triang_x_sr100,2))];
triang_x_sr100_pp = interp1(t_GRT_00, triang_x_sr100_00, 'linear', 'pp');

t_GIRF = Hx_sr100.t_axis;
dwelltime_girf = Hx_sr100.dt;
if t_GIRF(end) < t_GRT(end)
    t_GIRF_end = ceil(t_GRT(end)/dwelltime_girf);
    t_GIRF = ((0:t_GIRF_end)+0.5)*dwelltime_girf;
end
% Interpolate triangles to time grid of the GIRF
triang_x_sr100_tgirf = ppval(triang_x_sr100_pp, t_GIRF);

% Zero-filling to avoid Gibbs-artifacts
nExtra = round((1e6-size(triang_x_sr100_tgirf,1))/2);
triang_x_sr100_ZF = cat(1, zeros(nExtra,size(triang_x_sr100,2)), triang_x_sr100_tgirf, zeros(nExtra,size(triang_x_sr100,2)));

% Interpolate GSTF to zero-filled frequency grid
gstf_x_sr100_interp = interp1(Hx_sr100.f_axis, Hx_sr100.gstf(:,2), linspace(Hx_sr100.f_axis(1),Hx_sr100.f_axis(end),size(triang_x_sr100_ZF,1)).', 'linear');

% Calculate GSTF-predicted triangle waveforms
tmp1 = fft1d(triang_x_sr100_ZF,1);
tmp2 = repmat(gstf_x_sr100_interp, [1,size(triang_x_sr100,2)]);
triang_x_sr100_girf = real(ifft1d( tmp1.*tmp2, 1)) + Hx_sr100.fieldOffsets(2);
triang_x_sr100_girf = triang_x_sr100_girf(nExtra+1:end-nExtra,:);

% Calculate measured triangle evolutions from FID data
[triang_x_sr100_meas, ~] = calcMeasOutput_Bruker_avgComplex(triang_x_sr100_data, 13, 2, t_ADC);

%% Load spiral trajectories
for loop = 1:length(measurements)
    meas = measurements{loop};
    filepath = meas.filepath;

    if isfile(filepath)
        load(filepath);
    else
        disp(['Error: cannot find ',filepath]);
        break;
    end

    spiral_x = zeros(spiral_data.numPointsPerInterleave,spiral_data.numInterleaves);
    spiral_y = zeros(size(spiral_x));
    max_grad_strength = 1.5212; % T/m
    grad_raster_time = 8.1e-6;
    dwelltime_meas = 1/double(spiral_data.bandwidth);
    gamma = 267.513*10^6; %Hz/T

    % Calculate gradient waveform for each spiral interleaf
    for arm = 1:spiral_data.numInterleaves
        spiral_x(:,arm) = spiral_data.spiralShape1*spiral_data.spiralInterleaveCos(arm) - spiral_data.spiralShape2*spiral_data.spiralInterleaveSin(arm);
        spiral_y(:,arm) = spiral_data.spiralShape1*spiral_data.spiralInterleaveSin(arm) + spiral_data.spiralShape2*spiral_data.spiralInterleaveCos(arm);
        % spiralshapes are normalized to 1
    end

    % Put together complex gradient waveform and calculate frequency spectrum
    grad_cplx = spiral_data.spiralShape1(1:end-spiral_data.nPost) + 1i*spiral_data.spiralShape2(1:end-spiral_data.nPost);
    measurements{loop}.grad_cplx = grad_cplx;
    spectrum_grad = fft1d(grad_cplx,1);
    F_GRT = 1/grad_raster_time;
    f_GRT = linspace(-F_GRT/2, F_GRT/2, size(spectrum_grad,1));
    measurements{loop}.spectrum = spectrum_grad;

    measurements{loop}.spiralShape = spiral_data.spiralShape1;
    measurements{loop}.spiral_x = spiral_x;
    measurements{loop}.spiral_y = spiral_y;
    measurements{loop}.f_axis = f_GRT;

    % Integrate spiral waveforms to obtain k-space trajectory
    t_GRT_spir = ((1:1:(size(spiral_x,1)))-0.5)*grad_raster_time;
    t_eps = 1e-12;
    t_GRT_00 = [-t_eps, t_GRT_spir(1)-t_eps, t_GRT_spir, t_GRT_spir(end)+t_eps, t_GRT_spir(end)+grad_raster_time+t_eps];
    
    spiral_x_00 = [zeros(2,size(spiral_x,2)); spiral_x; zeros(2,size(spiral_x,2))];
    spiral_y_00 = [zeros(2,size(spiral_x,2)); spiral_y; zeros(2,size(spiral_x,2))];

    spiral_x_pp = interp1(t_GRT_00, spiral_x_00*max_grad_strength, 'linear', 'pp');
    spiral_y_pp = interp1(t_GRT_00, spiral_y_00*max_grad_strength, 'linear', 'pp');

    traj_x_pp = fnint(spiral_x_pp);
    traj_y_pp = fnint(spiral_y_pp);
    traj_x = gamma/(2*pi)*ppval(traj_x_pp, t_GRT_spir)/1000;
    traj_y = gamma/(2*pi)*ppval(traj_y_pp, t_GRT_spir)/1000;
    traj = (traj_x + 1i*traj_y);
    traj = traj./max(abs(traj(:,1:end-spiral_data.nPost)),[],'all')/2;
    measurements{loop}.traj = traj;
    measurements{loop}.nPost = spiral_data.nPost;
end

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

%% Plot Figure 1
figure('Units','centimeters','Position',[1 1 17.56 24]);%19

ax1 = subplot('Position',[0.08 0.81 0.42 0.18]);
plot(Hx_sr100.f_axis/1000, abs(Hx_sr100.gstf(:,2)), 'DisplayName','|GSTF_x|','LineWidth',1.2);
hold on;
plot(Hy_sr100.f_axis/1000, abs(Hy_sr100.gstf(:,2)), 'DisplayName','|GSTF_y|','LineWidth',1.2,'Color',green);
plot(Hz_sr100.f_axis/1000, abs(Hz_sr100.gstf(:,2)), 'DisplayName','|GSTF_z|','LineWidth',1.2,'Color',violet);
xlabel('Frequency (kHz)');
ylabel('Magnitude (1)');
leg = legend('numColumns',3,'Location','north');
leg.ItemTokenSize = [15 5];
ylim([0 1.3]);
set(gca,'FontName','Times','Fontsize',9);

text(-33,1.25,'(A)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(26,1.25,'(B)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(-33,-0.37,'(C)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(6,-0.37,'(D)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(44,-0.37,'(E)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(-33,-2.6,'(F)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(6.5,-2.6,'(G)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(47.5,-2.6,'(H)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(-33,-4.1,'(I)','FontName','Arial','Fontsize',11,'FontWeight','bold');
text(31.5,-4.1,'(J)','FontName','Arial','Fontsize',11,'FontWeight','bold');

ax2 = subplot('Position',[0.56 0.81 0.42 0.18]);
plot(Hx_sr100.f_axis/1000, angle(Hx_sr100.gstf(:,2)), 'DisplayName','\angle GSTF_x','LineWidth',1.2);
hold on;
plot(Hy_sr100.f_axis/1000, angle(Hy_sr100.gstf(:,2)), 'DisplayName','\angle GSTF_y','LineWidth',1.2,'Color',green);
plot(Hz_sr100.f_axis/1000, angle(Hz_sr100.gstf(:,2)), 'DisplayName','\angle GSTF_z','LineWidth',1.2,'Color',violet);
xlabel('Frequency (kHz)');
ylabel('Phase (rad)');
leg = legend('Location','northeast');
leg.ItemTokenSize = [13 5];
ylim([-2.5 2.5]);
set(gca,'FontName','Times','Fontsize',9);

ax3 = subplot('Position',[0.08 0.51 0.26 0.24]);
plot(t_GRT*1e3, triang_x_sr100(:,end:-1:1), 'LineWidth',1.1,'Color',[0.5 0.5 0.5]+0.1);
hold on;
plot([-0.001 0],[0 0]);plot([-0.001 0],[0 0]);
plot(t_ADC*1e3, squeeze(triang_x_sr100_meas(2,:,end:-1:1)), 'LineWidth',1.2);
xline(0.10265,'--');
xline(0.15935,'--');
ylabel('Gradient (T/m)');
xlabel('Time (ms)');
ylim([0 2]);
set(gca,'FontName','Times','Fontsize',9);

ax4 = subplot('Position',[0.39 0.51 0.26 0.24]);
plot(t_ADC*1e3, squeeze(triang_x_sr100_meas(2,:,end:-1:1)), 'LineWidth',1.1,'Color',[0.5 0.5 0.5]+0.1);
hold on;
plot([-0.001 0],[0 0]);plot([-0.001 0],[0 0]);
plot(t_GIRF*1e3, triang_x_sr100_girf(:,end:-1:1),'--', 'LineWidth',1.2);
xline(0.10265,'--');
xline(0.15935,'--');
xlabel('Time (ms)');
ylim([0 2]);
set(gca,'FontName','Times','Fontsize',9);

ax5 = subplot('Position',[0.72 0.51 0.26 0.24]);
plot([-0.001 0],[0 0]);hold on;
plot(t_ADC*1e3, squeeze(triang_x_sr100_meas(2,:,end:-1:1))-triang_x_sr100_girf(1:length(t_ADC),end:-1:1), 'LineWidth',1.2);
xline(0.10265,'--');
xline(0.15935,'--');
xlabel('Time (ms)');
ylim([-0.043 0.043]);
set(gca,'FontName','Times','Fontsize',9);

linkaxes([ax1 ax2], 'x');
xlim(ax1,[-25 25]);

linkaxes([ax3 ax4 ax5],'x');
xlim(ax3,[0 0.35]);

ax6 = subplot('Position',[0.08 0.3 0.22 0.16]);
plot(meas1.traj(1,1:end-meas1.nPost).', '--', 'LineWidth',1.8,'Color',orange);
axis image;
xlabel('k_x');
ylabel('k_y');

ax7 = subplot('Position',[0.42 0.3 0.22 0.16]);
plot(meas2.traj(1,1:end-meas2.nPost).', '-.', 'LineWidth',1.8,'Color',yellow);
axis image;
xlabel('k_x');

ax8 = subplot('Position',[0.76 0.3 0.22 0.16]);
plot(meas3.traj(1,1:end-meas3.nPost).', 'LineWidth',1.1,'Color',lblue);
axis image;
xlabel('k_x');

linkaxes([ax6 ax7 ax8],'xy');
xlim(ax6,[-1 1]*0.55);
ylim(ax6,[-1 1]*0.55);

ax9 = subplot('Position',[0.08 0.048 0.46 0.2]);
plot(Hx_sr100.f_axis/1000, abs(Hx_sr100.gstf(:,2)), 'DisplayName','|GSTF_x|','LineWidth',1.2);
hold on;
plot(meas1.f_axis/1000, abs(meas1.spectrum)/max(abs(meas1.spectrum)),'--', 'DisplayName','|FFT(spiral_9_6)|','LineWidth',1.2,'Color',orange);
plot(Hy_sr100.f_axis/1000, abs(Hy_sr100.gstf(:,2)), 'DisplayName','|GSTF_y|','LineWidth',1.2,'Color',green);
plot(meas2.f_axis/1000, abs(meas2.spectrum)/max(abs(meas2.spectrum)),'-.', 'DisplayName','|FFT(spiral_1_6)|','LineWidth',1.2,'Color',yellow);
plot(Hz_sr100.f_axis/1000, abs(Hz_sr100.gstf(:,2)), 'DisplayName','|GSTF_z|','LineWidth',1.2,'Color',violet);
plot(meas3.f_axis/1000, abs(meas3.spectrum)/max(abs(meas3.spectrum)),':', 'DisplayName','|FFT(spiral_3)|','LineWidth',1.5,'Color',lblue);
xlabel('Frequency (kHz)');
ylabel('Magnitude (1)');
leg = legend('numColumns',3,'Location','north');
leg.ItemTokenSize = [15 5];
ylim([0 1.55]);
set(gca,'FontName','Times','Fontsize',9);
xlim([-25 25]);
rectangle(ax9, 'Position',[0 0 4 1.05], 'EdgeColor',orange);
rectangle(ax9, 'Position',[0 0.8 4 0.22], 'EdgeColor',blue);
text(4.5,1,'(J)','FontName','Arial','Fontsize',8);

ax10 = subplot('Position',[0.62 0.048 0.3 0.2]);
yyaxis right;
plot(meas1.f_axis/1000, abs(meas1.spectrum)/max(abs(meas1.spectrum)),'--', 'DisplayName','|FFT(spiral_9_6)|','LineWidth',1.2,'Color',orange);
hold on;
plot(meas2.f_axis/1000, abs(meas2.spectrum)/max(abs(meas2.spectrum)),'-.', 'DisplayName','|FFT(spiral_1_6)|','LineWidth',1.2,'Color',yellow);
plot(meas3.f_axis/1000, abs(meas3.spectrum)/max(abs(meas3.spectrum)),':', 'DisplayName','|FFT(spiral_3)|','LineWidth',1.5,'Color',lblue);
ylim([0 1.05]);
ylabel('spiral spectral intensity (a.u.)');
yyaxis left;
plot(Hx_sr100.f_axis/1000, abs(Hx_sr100.gstf(:,2)),'-', 'DisplayName','|GSTF_x|','LineWidth',1.2, 'Color',blue);
plot(Hy_sr100.f_axis/1000, abs(Hy_sr100.gstf(:,2)),'-', 'DisplayName','|GSTF_y|','LineWidth',1.2,'Color',green);
plot(Hz_sr100.f_axis/1000, abs(Hz_sr100.gstf(:,2)),'-', 'DisplayName','|GSTF_z|','LineWidth',1.2,'Color',violet);
xlabel('Frequency (kHz)');
ylabel('GSTF magnitude (1)');
ylim([0.8 1.02]);
set(gca,'FontName','Times','Fontsize',9);
xlim([0 4]);









