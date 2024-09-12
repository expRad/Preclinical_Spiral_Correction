% Copyright (c) 2024 Hannah Scholten

clear all;
% Change current directory to that of this .m file
mfile_name          = mfilename('fullpath');
if contains(mfile_name,'LiveEditorEvaluationHelper')
    mfile_name = matlab.desktop.editor.getActiveFilename;
end
clear all;
% Change current directory to that of this .m file
mfile_name          = mfilename('fullpath');
if contains(mfile_name,'LiveEditorEvaluationHelper')
    mfile_name = matlab.desktop.editor.getActiveFilename;
end
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);
% Add necessary scripts to the MATLAB path
folder1 = "..\Spiral_Recon_NUFFT";
addpath(folder1);
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
recos_meas = zeros(matrix_size, matrix_size, length(measurements));
recos_meas_wB0 = zeros(matrix_size, matrix_size, length(measurements));

%% Reconstruct images for each spiral trajectory
for loop = 1:length(measurements)
    disp(['loop = ',num2str(loop)]);
    meas = measurements{loop};
    filepath = meas.filepath;
    
    %% Load spiral data
    disp('Read measurement data.')

    if isfile(filepath)
        load(filepath);
    else
        disp(['Error: cannot find ',filepath]);
        break
    end
    raw = spiral_data.raw;
    nPost = spiral_data.nPost; % data points discarded for reconstruction
    measurements{loop}.nPost = nPost;
    disp(['NUMBER of SPIRAL INTERLEAVES: ',num2str(spiral_data.numInterleaves)]);

    %% reconstruct images with measured trajectory
    trajx_meas = reshape(spiral_data.trajKx_meas, [spiral_data.datapoints, spiral_data.numInterleaves]);
    trajy_meas = reshape(spiral_data.trajKy_meas, [spiral_data.datapoints, spiral_data.numInterleaves]);
    trajz_meas = reshape(spiral_data.trajKz_meas, [spiral_data.datapoints, spiral_data.numInterleaves]);

    % measured phases due to B0 eddy currents
    trajBx_meas = reshape(spiral_data.trajBx_meas, [spiral_data.datapoints, spiral_data.numInterleaves]);
    trajBy_meas = reshape(spiral_data.trajBy_meas, [spiral_data.datapoints, spiral_data.numInterleaves]);
    trajBz_meas = reshape(spiral_data.trajBz_meas, [spiral_data.datapoints, spiral_data.numInterleaves]);
    trajB0_meas = trajBx_meas + trajBy_meas + trajBz_meas;

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

    FT = cGrid(traj_meas(1:end-nPost,:)/traj_factor,matrix_size); % NUFFT operator
    % reconstruction without B0 eddy current correction
    reco_meas = FT'*raw(1:end-nPost,:);
    % reconstruction with B0 eddy current correction
    reco_meas_wB0 = FT'*(raw(1:end-nPost,:).*exp(-1i*trajB0_meas(1:end-nPost,:)));

    %% Normalize reconstructed images
    s_meas_wB0 = sum(abs(reco_meas).*abs(reco_meas_wB0), 'all') / sum(abs(reco_meas_wB0).*abs(reco_meas_wB0), 'all');
    
    recos_meas(:,:,loop) = reco_meas;
    recos_meas_wB0(:,:,loop) = s_meas_wB0*reco_meas_wB0;

end % loop

meas1 = measurements{1};
meas2 = measurements{2};
meas3 = measurements{3};

%% Plot difference between reconstructions with and without B0-correction (measured trajectory)
figure('Units','centimeters','Position',[0 0 17.56 24],'Color','k');
colormap gray;
c_lims1 = [0 0.2];
c_lims2 = [0 0.1];
c_lims3 = [0 0.05];
c_lims1a = [-0.2 0.2];
c_lims2a = [-0.1 0.1];
c_lims3a = [-0.05 0.05];

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

ax1 = subplot('Position',[x1 y1 dx dy]);
imagesc(abs(recos_meas(:,:,1)), c_lims1);
axis image; axis off;
xticklabels([]); yticklabels([]);
title('96 interleaves','Color','w');
ylabel('without B_0','Color','w','FontWeight','bold');
text(-30,96, 'without B_0','Color','w','FontWeight','bold','HorizontalAlignment','center','Rotation',90,'FontName','Times','Fontsize',10);
set(gca,'FontName','Times','Fontsize',9);
text(0,10,'(A)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

ax2 = subplot('Position',[x2 y1 dx dy]);
imagesc(abs(recos_meas(:,:,2)), c_lims2);
axis image; axis off;
title('16 interleaves','Color','w');
set(gca,'FontName','Times','Fontsize',9);
text(0,10,'(B)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

ax3 = subplot('Position',[x3 y1 dx dy]);
imagesc(abs(recos_meas(:,:,3)), c_lims3);
axis image; axis off;
title('3 interleaves','Color','w');
set(gca,'FontName','Times','Fontsize',9);
text(0,10,'(C)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

ax4 = subplot('Position',[x1 y2 dx dy]);
imagesc(abs(recos_meas_wB0(:,:,1)), c_lims1);
axis image; axis off;
xticklabels([]); yticklabels([]);
ylabel('delay','Color','w','FontWeight','bold');
text(-30,96, 'with B_0','Color','w','FontWeight','bold','HorizontalAlignment','center','Rotation',90,'FontName','Times','Fontsize',10);
set(gca,'FontName','Times','Fontsize',9);
text(0,10,'(D)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

ax5 = subplot('Position',[x2 y2 dx dy]);
imagesc(abs(recos_meas_wB0(:,:,2)), c_lims2);
axis image; axis off;
text(0,10,'(E)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

ax6 = subplot('Position',[x3 y2 dx dy]);
imagesc(abs(recos_meas_wB0(:,:,3)), c_lims3);
axis image; axis off;
text(0,10,'(F)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

ax7 = subplot('Position',[x1 y3 dx dy]);
imagesc(abs(recos_meas_wB0(:,:,1))-abs(recos_meas(:,:,1)), c_lims1a);
axis image; axis off;
xticklabels([]); yticklabels([]);
ylabel('GSTF','Color','w','FontWeight','bold');
text(-30,96, 'difference','Color','w','FontWeight','bold','HorizontalAlignment','center','Rotation',90,'FontName','Times','Fontsize',10);
set(gca,'FontName','Times','Fontsize',9);
text(0,10,'(G)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

ax8 = subplot('Position',[x2 y3 dx dy]);
imagesc(abs(recos_meas_wB0(:,:,2))-abs(recos_meas(:,:,2)), c_lims2a);
axis image; axis off;
text(0,10,'(H)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

ax9 = subplot('Position',[x3 y3 dx dy]);
imagesc(abs(recos_meas_wB0(:,:,3))-abs(recos_meas(:,:,3)), c_lims3a);
axis image; axis off;
text(0,10,'(I)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

ax10 = subplot('Position',[x1 y4 dx dy]);
imagesc(abs(recos_meas_wB0(:,:,1))-abs(recos_meas(:,:,1)), c_lims1a/10);
axis image; axis off;
xticklabels([]); yticklabels([]);
ylabel('GSTF+delay','Color','w','FontWeight','bold');
text(-30,96, 'colorscale/10','Color','w','FontWeight','bold','HorizontalAlignment','center','Rotation',90,'FontName','Times','Fontsize',10);
set(gca,'FontName','Times','Fontsize',9);
text(0,10,'(J)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

ax11 = subplot('Position',[x2 y4 dx dy]);
imagesc(abs(recos_meas_wB0(:,:,2))-abs(recos_meas(:,:,2)), c_lims2a/10);
axis image; axis off;
text(0,10,'(K)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

ax12 = subplot('Position',[x3 y4 dx dy]);
imagesc(abs(recos_meas_wB0(:,:,3))-abs(recos_meas(:,:,3)), c_lims3a/10);
axis image; axis off;
text(0,10,'(L)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

ax13 = subplot('Position',[x1 y5 dx dy]);
imagesc(abs(recos_meas_wB0(:,:,1))-abs(recos_meas(:,:,1)), c_lims1a/100);
axis image; axis off;
xticklabels([]); yticklabels([]);
ylabel('measured','Color','w','FontWeight','bold');
text(-30,96, 'colorscale/100','Color','w','FontWeight','bold','HorizontalAlignment','center','Rotation',90,'FontName','Times','Fontsize',10);
set(gca,'FontName','Times','Fontsize',9);
text(0,10,'(M)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

ax14 = subplot('Position',[x2 y5 dx dy]);
imagesc(abs(recos_meas_wB0(:,:,2))-abs(recos_meas(:,:,2)), c_lims2a/100);
axis image; axis off;
text(0,10,'(N)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

ax15 = subplot('Position',[x3 y5 dx dy]);
imagesc(abs(recos_meas_wB0(:,:,3))-abs(recos_meas(:,:,3)), c_lims3a/100);
axis image; axis off;
text(0,10,'(O)','FontName','Arial','Fontsize',11,'FontWeight','bold','Color','w');

set(gcf, 'InvertHardcopy', 'off');

%% Plot difference between reconstructions with measured trajectories

reco_ref = squeeze(recos_meas(:,:,3));
reco_96 = squeeze(recos_meas(:,:,1));
reco_16 = squeeze(recos_meas(:,:,2));
reco_3 = squeeze(recos_meas(:,:,3));

s_96 = sum(abs(reco_ref).*abs(reco_96), 'all') / sum(abs(reco_96).*abs(reco_96), 'all');
s_16 = sum(abs(reco_ref).*abs(reco_16), 'all') / sum(abs(reco_16).*abs(reco_16), 'all');
s_3 = sum(abs(reco_ref).*abs(reco_3), 'all') / sum(abs(reco_3).*abs(reco_3), 'all');

c_lims = [0 0.08];

figure('Color','k');
colormap gray;

subplot(2,3,1);
imagesc(s_96*abs(reco_96),c_lims);
axis image; axis off;
title('96 interleaves','Color','w');
colorbar('Color','w');

subplot(2,3,2);
imagesc(s_16*abs(reco_16),c_lims);
axis image; axis off;
title('16 interleaves','Color','w');
colorbar('Color','w');

subplot(2,3,3);
imagesc(s_3*abs(reco_3),c_lims);
axis image; axis off;
title('3 interleaves','Color','w');
colorbar('Color','w');

subplot(2,3,4);
imagesc(abs(s_96*abs(reco_96)-abs(reco_ref)),c_lims);
axis image; axis off;
title('difference to 3 interleaves','Color','w');
colorbar('Color','w');

subplot(2,3,5);
imagesc(abs(s_16*abs(reco_16)-abs(reco_ref)),c_lims);
axis image; axis off;
title('difference to 3 interleaves','Color','w');
colorbar('Color','w');

subplot(2,3,6);
imagesc(abs(s_3*abs(reco_3)-abs(reco_ref)),c_lims);
axis image; axis off;
title('difference to 3 interleaves','Color','w');
colorbar('Color','w');









