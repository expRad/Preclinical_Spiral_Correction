clear all;

% Load data
data_p96 = load('..\phantom_data\RMSEs_and_delays_p96.mat');
data_p16 = load('..\phantom_data\RMSEs_and_delays_p16.mat');
data_p3 = load('..\phantom_data\RMSEs_and_delays_p03.mat');
data_mouse = load('..\invivo_data\RMSEs_and_delays_m96.mat');

data = {data_p96, data_p16, data_p3, data_mouse};
vals = zeros(2,4);

%% Fit polynomial to error-vs-delay for each acquisition
for i=1:length(data)
    delays_del = data{i}.delays_del;
    delays_girf = data{i}.delays_girf;
    sum_abs_err_del = data{i}.sum_abs_err_del;
    sum_abs_err_girf = data{i}.sum_abs_err_girf;
    % sum_gradx_err_girf = data{i}.sum_gradx_err_girf;
    % sum_grady_err_girf = data{i}.sum_grady_err_girf;
    % sum_gradz_err_girf = data{i}.sum_gradz_err_girf;
    % sum_gradx_err_del = data{i}.sum_gradx_err_del;
    % sum_grady_err_del = data{i}.sum_grady_err_del;
    % sum_gradz_err_del = data{i}.sum_gradz_err_del;
    % sum_traj_err_girf = data{i}.sum_traj_err_girf;
    % sum_traj_err_del = data{i}.sum_traj_err_del;
    
    % Initial guess for the parameters
    initial_guess = [1, 1, 1, 1, 1];
    % polynomial of 4th order
    quad_function = @(params, x) params(1)*x.^4 + params(2)*x.^3 + params(3)*x.^2 + params(4)*x + params(5);
    
    % Fit the data using lsqcurvefit
    params_del_abs = lsqcurvefit(quad_function, initial_guess, delays_del(1:11)*1e5, sum_abs_err_del(1:11)*100);
    params_girf_abs = lsqcurvefit(quad_function, initial_guess, delays_girf(1:11)*1e5, sum_abs_err_girf(1:11)*100);

    % % for gradient waveform RMSEs
    % params_del_gradx = lsqcurvefit(quad_function, initial_guess, delays_del(1:11)*1e5, sum_gradx_err_del(1:11)*100);
    % params_del_grady = lsqcurvefit(quad_function, initial_guess, delays_del(1:11)*1e5, sum_grady_err_del(1:11)*100);
    % params_del_gradz = lsqcurvefit(quad_function, initial_guess, delays_del(1:11)*1e5, sum_gradz_err_del(1:11)*100);
    % params_girf_gradx = lsqcurvefit(quad_function, initial_guess, delays_girf(1:11)*1e5, sum_gradx_err_girf(1:11)*100);
    % params_girf_grady = lsqcurvefit(quad_function, initial_guess, delays_girf(1:11)*1e5, sum_grady_err_girf(1:11)*100);
    % params_girf_gradz = lsqcurvefit(quad_function, initial_guess, delays_girf(1:11)*1e5, sum_gradz_err_girf(1:11)*100);
    % % for trajectory RMSE
    % params_del_traj = lsqcurvefit(quad_function, initial_guess, delays_del(1:11)*1e5, sum_traj_err_del(1:11)*100);
    % params_girf_traj = lsqcurvefit(quad_function, initial_guess, delays_girf(1:11)*1e5, sum_traj_err_girf(1:11)*100);

    % Find extrema for delay correction
    extrema_del_abs = roots([4*params_del_abs(1) 3*params_del_abs(2) 2*params_del_abs(3) params_del_abs(4)]);
    data{i}.extrema_del_abs = extrema_del_abs*1e-5;
    for j=1:3
        if imag(extrema_del_abs(j))==0
            if 12*params_del_abs(1)*extrema_del_abs(j)^2 + 6*params_del_abs(2)*extrema_del_abs(j) + 2*params_del_abs(3) > 0
                disp(['Min del = ',num2str(extrema_del_abs(j)*1e-5)]);
                data{i}.min_del = extrema_del_abs(j)*1e-5;
            end
        end
    end
    vals(1,i) = data{i}.min_del;

    % % for gradient waveform RMSEs
    % extrema_del_gradx = roots([4*params_del_gradx(1) 3*params_del_gradx(2) 2*params_del_gradx(3) params_del_gradx(4)]);
    % disp(['Min del gradx = ',num2str(extrema_del_gradx(2)*1e-5)]);
    % data{i}.extrema_del_gradx = extrema_del_gradx*1e-5;
    % extrema_del_grady = roots([4*params_del_grady(1) 3*params_del_grady(2) 2*params_del_grady(3) params_del_grady(4)]);
    % disp(['Min del grady = ',num2str(extrema_del_grady(2)*1e-5)]);
    % data{i}.extrema_del_grady = extrema_del_grady*1e-5;
    % extrema_del_gradz = roots([4*params_del_gradz(1) 3*params_del_gradz(2) 2*params_del_gradz(3) params_del_gradz(4)]);
    % disp(['Min del gradz = ',num2str(extrema_del_gradz(2)*1e-5)]);
    % data{i}.extrema_del_gradz = extrema_del_gradz*1e-5;
    % % for trajectory RMSE
    % extrema_del_traj = roots([4*params_del_traj(1) 3*params_del_traj(2) 2*params_del_traj(3) params_del_traj(4)]);
    % disp(['Min del traj = ',num2str(extrema_del_traj(2)*1e-5)]);
    % data{i}.extrema_del_traj = extrema_del_traj*1e-5;

    extrema_girf_abs = roots([4*params_girf_abs(1) 3*params_girf_abs(2) 2*params_girf_abs(3) params_girf_abs(4)]);
    disp(['Min girf abs = ',num2str(extrema_girf_abs(2)*1e-5)]);
    data{i}.extrema_girf_abs = extrema_girf_abs*1e-5;
    for j=1:3
        if imag(extrema_girf_abs(j))==0
            if 12*params_girf_abs(1)*extrema_girf_abs(j)^2 + 6*params_girf_abs(2)*extrema_girf_abs(j) + 2*params_girf_abs(3) > 0
                disp(['Min girf = ',num2str(extrema_girf_abs(j)*1e-5)]);
                data{i}.min_girf = extrema_girf_abs(j)*1e-5;
            end
        end
    end
    vals(2,i) = data{i}.min_girf;

    % % for gradient waveform RMSEs
    % extrema_girf_gradx = roots([4*params_girf_gradx(1) 3*params_girf_gradx(2) 2*params_girf_gradx(3) params_girf_gradx(4)]);
    % disp(['Min girf gradx = ',num2str(extrema_girf_gradx(2)*1e-5)]);
    % data{i}.extrema_girf_gradx = extrema_girf_gradx*1e-5;
    % extrema_girf_grady = roots([4*params_girf_grady(1) 3*params_girf_grady(2) 2*params_girf_grady(3) params_girf_grady(4)]);
    % disp(['Min girf grady = ',num2str(extrema_girf_grady(2)*1e-5)]);
    % data{i}.extrema_girf_grady = extrema_girf_grady*1e-5;
    % extrema_girf_gradz = roots([4*params_girf_gradz(1) 3*params_girf_gradz(2) 2*params_girf_gradz(3) params_girf_gradz(4)]);
    % disp(['Min girf gradz = ',num2str(extrema_girf_gradz(2)*1e-5)]);
    % data{i}.extrema_girf_gradz = extrema_girf_gradz*1e-5;
    % % for trajectory RMSE
    % extrema_girf_traj = roots([4*params_girf_traj(1) 3*params_girf_traj(2) 2*params_girf_traj(3) params_girf_traj(4)]);
    % disp(['Min girf traj = ',num2str(extrema_girf_traj(2)*1e-5)]);
    % data{i}.extrema_girf_traj = extrema_girf_traj*1e-5;

    delays_girf_fine = ((-100:0))*1e-7*1e5;
    delays_del_fine = ((-100:0)-150)*1e-7*1e5;

    data{i}.delays_girf_fine = delays_girf_fine*1e-5;
    data{i}.delays_del_fine = delays_del_fine*1e-5;

    fit_del_abs = params_del_abs(1)*delays_del_fine.^4 + params_del_abs(2)*delays_del_fine.^3 + params_del_abs(3)*delays_del_fine.^2 + params_del_abs(4)*delays_del_fine + params_del_abs(5);
    % fit_del_gradx = params_del_gradx(1)*delays_del_fine.^4 + params_del_gradx(2)*delays_del_fine.^3 + params_del_gradx(3)*delays_del_fine.^2 + params_del_gradx(4)*delays_del_fine + params_del_gradx(5);
    % fit_del_grady = params_del_grady(1)*delays_del_fine.^4 + params_del_grady(2)*delays_del_fine.^3 + params_del_grady(3)*delays_del_fine.^2 + params_del_grady(4)*delays_del_fine + params_del_grady(5);
    % fit_del_gradz = params_del_gradz(1)*delays_del_fine.^4 + params_del_gradz(2)*delays_del_fine.^3 + params_del_gradz(3)*delays_del_fine.^2 + params_del_gradz(4)*delays_del_fine + params_del_gradz(5);
    % fit_del_traj = params_del_traj(1)*delays_del_fine.^4 + params_del_traj(2)*delays_del_fine.^3 + params_del_traj(3)*delays_del_fine.^2 + params_del_traj(4)*delays_del_fine + params_del_traj(5);

    fit_girf_abs = params_girf_abs(1)*delays_girf_fine.^4 + params_girf_abs(2)*delays_girf_fine.^3 + params_girf_abs(3)*delays_girf_fine.^2 + params_girf_abs(4)*delays_girf_fine + params_girf_abs(5);
    % fit_girf_gradx = params_girf_gradx(1)*delays_girf_fine.^4 + params_girf_gradx(2)*delays_girf_fine.^3 + params_girf_gradx(3)*delays_girf_fine.^2 + params_girf_gradx(4)*delays_girf_fine + params_girf_gradx(5);
    % fit_girf_grady = params_girf_grady(1)*delays_girf_fine.^4 + params_girf_grady(2)*delays_girf_fine.^3 + params_girf_grady(3)*delays_girf_fine.^2 + params_girf_grady(4)*delays_girf_fine + params_girf_grady(5);
    % fit_girf_gradz = params_girf_gradz(1)*delays_girf_fine.^4 + params_girf_gradz(2)*delays_girf_fine.^3 + params_girf_gradz(3)*delays_girf_fine.^2 + params_girf_gradz(4)*delays_girf_fine + params_girf_gradz(5);
    % fit_girf_traj = params_girf_traj(1)*delays_girf_fine.^4 + params_girf_traj(2)*delays_girf_fine.^3 + params_girf_traj(3)*delays_girf_fine.^2 + params_girf_traj(4)*delays_girf_fine + params_girf_traj(5);

    fit_del_abs = fit_del_abs/100;
    % fit_del_gradx = fit_del_gradx/100;
    % fit_del_grady = fit_del_grady/100;
    % fit_del_gradz = fit_del_gradz/100;
    % fit_del_traj = fit_del_traj/100;

    fit_girf_abs = fit_girf_abs/100;
    % fit_girf_gradx = fit_girf_gradx/100;
    % fit_girf_grady = fit_girf_grady/100;
    % fit_girf_gradz = fit_girf_gradz/100;
    % fit_girf_traj = fit_girf_traj/100;

    data{i}.fit_del_abs = fit_del_abs;
    % data{i}.fit_del_gradx = fit_del_gradx;
    % data{i}.fit_del_grady = fit_del_grady;
    % data{i}.fit_del_gradz = fit_del_gradz;
    % data{i}.fit_del_traj = fit_del_traj;

    data{i}.fit_girf_abs = fit_girf_abs;
    % data{i}.fit_girf_gradx = fit_girf_gradx;
    % data{i}.fit_girf_grady = fit_girf_grady;
    % data{i}.fit_girf_gradz = fit_girf_gradz;
    % data{i}.fit_girf_traj = fit_girf_traj;
end
%% Define colors
blue = [0 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];
violet = [0.4940 0.1840 0.5560];
green = [0.4660 0.6740 0.1880];
yellow = [0.9290 0.6940 0.1250];
lblue = [0.3010 0.7450 0.9330];
dred = [0.6350 0.0780 0.1840];

%% Plot Figure 2
if 1
    dx = 0.86;
    dy = 0.13;
        
    fig = figure('Units','centimeters','Position',[5 5 17.56 17]);
    i = 1;
    ax1 = subplot('Position',[0.11 0.782 dx dy]);
    plot(data{i}.delays_del*1e6, data{i}.sum_abs_err_del, 'o', 'DisplayName','RMSE (delay corr.)');
    hold on;
    plot(data{i}.delays_girf*1e6, data{i}.sum_abs_err_girf, 'o', 'DisplayName','RMSE (GSTF+delay corr.)','Color',yellow);
    plot(data{i}.delays_del_fine*1e6, data{i}.fit_del_abs, '-', 'DisplayName','fit (delay corr.)','LineWidth',1.2,'Color',orange);
    plot(data{i}.delays_girf_fine*1e6, data{i}.fit_girf_abs, '-', 'DisplayName','fit (GSTF+delay corr.)','LineWidth',1.2);
    xline(data{i}.min_del*1e6,'--', 'DisplayName','optimum delay for isotropic delay correction','LineWidth',1.2);
    xline(data{i}.min_girf*1e6,'-.', 'DisplayName','optimum delay for GSTF+delay correction','LineWidth',1.2);
    ylabel('nRMSE (a.u.)');
    leg = legend('numColumns',3,'Position',[0.11 0.936 0.86 0.0235]);
    leg.ItemTokenSize = [17 5];
    ylim([0.005 0.04]);
    set(gca,'FontName','Times','Fontsize',9);
    text(-29.2,0.043,'(A)','FontName','Arial','Fontsize',11,'FontWeight','bold');
    text(-29.2,-0.005,'(B)','FontName','Arial','Fontsize',11,'FontWeight','bold');
    text(-29.2,-0.058,'(C)','FontName','Arial','Fontsize',11,'FontWeight','bold');
    text(-29.2,-0.104,'(D)','FontName','Arial','Fontsize',11,'FontWeight','bold');
    text(-29.2,-0.163,'(E)','FontName','Arial','Fontsize',11,'FontWeight','bold');
    text(-20,0.06, 'delay correction','FontName','Times','Fontsize',9,'FontWeight','bold','HorizontalAlignment','center');
    text(-5,0.06, 'GSTF+delay correction','FontName','Times','Fontsize',9,'FontWeight','bold','HorizontalAlignment','center');
    i = i+1;
    
    ax2 = subplot('Position',[0.11 0.613 dx dy]);
    plot(data{i}.delays_del*1e6, data{i}.sum_abs_err_del, 'o', 'DisplayName','RMSE delay');
    hold on;
    plot(data{i}.delays_del_fine*1e6, data{i}.fit_del_abs, '-', 'DisplayName','fit delay','LineWidth',1.2);
    xline(data{i}.min_del*1e6,'--', 'DisplayName','minimum delay','LineWidth',1.2);
    plot(data{i}.delays_girf*1e6, data{i}.sum_abs_err_girf, 'o', 'DisplayName','RMSE GSTF');
    plot(data{i}.delays_girf_fine*1e6, data{i}.fit_girf_abs, '-', 'DisplayName','fit GSTF','LineWidth',1.2);
    xline(data{i}.min_girf*1e6,'-.', 'DisplayName','minimum GSTF','LineWidth',1.2);
    ylabel('nRMSE (a.u.)');
    leg.ItemTokenSize = [17 5];
    ylim([0.005 0.02]);
    set(gca,'FontName','Times','Fontsize',9);
    i = i+1;
    
    ax3 = subplot('Position',[0.11 0.418 dx dy]);
    plot(data{i}.delays_del*1e6, data{i}.sum_abs_err_del, 'o', 'DisplayName','RMSE delay');
    hold on;
    plot(data{i}.delays_del_fine*1e6, data{i}.fit_del_abs, '-', 'DisplayName','fit delay','LineWidth',1.2);
    xline(data{i}.min_del*1e6,'--', 'DisplayName','minimum delay','LineWidth',1.2);
    plot(data{i}.delays_girf*1e6, data{i}.sum_abs_err_girf, 'o', 'DisplayName','RMSE GSTF');
    plot(data{i}.delays_girf_fine*1e6, data{i}.fit_girf_abs, '-', 'DisplayName','fit GSTF','LineWidth',1.2);
    xline(data{i}.min_girf*1e6,'-.', 'DisplayName','minimum GSTF','LineWidth',1.2);
    ylabel('nRMSE (a.u.)');
    leg.ItemTokenSize = [17 5];
    ylim([0.006 0.013]);
    set(gca,'FontName','Times','Fontsize',9);
    i = i+1;
    
    ax4 = subplot('Position',[0.11 0.245 dx dy]);
    plot(data{i}.delays_del*1e6, data{i}.sum_abs_err_del, 'o', 'DisplayName','RMSE delay');
    hold on;
    plot(data{i}.delays_del_fine*1e6, data{i}.fit_del_abs, '-', 'DisplayName','fit delay','LineWidth',1.2);
    xline(data{i}.min_del*1e6,'--', 'DisplayName','minimum delay','LineWidth',1.2);
    plot(data{i}.delays_girf*1e6, data{i}.sum_abs_err_girf, 'o', 'DisplayName','RMSE GSTF');
    plot(data{i}.delays_girf_fine*1e6, data{i}.fit_girf_abs, '-', 'DisplayName','fit GSTF','LineWidth',1.2);
    xline(data{i}.min_girf*1e6,'-.', 'DisplayName','minimum GSTF','LineWidth',1.2);
    xlabel('Delay \tau (\mus)');
    ylabel('nRMSE (a.u.)');
    leg.ItemTokenSize = [17 5];
    set(gca,'FontName','Times','Fontsize',9);
    i = i+1;
    
    linkaxes([ax1 ax2 ax3 ax4],'x');
    xlim(ax1,[-26 1]);

    vals = vals.';
    Mu = char(956);
    Astr = reshape(cellstr(num2str(vals(:)*1e6,'%.1f')),size(vals));
    for i=1:4
        for j=1:2
            Astr{i,j} = ['              ',Astr{i,j},' ',Mu,'s'];
        end
    end
    correction = {'delay correction','GSTF+delay correction'};
    vars = {'96 interleaves, phantom (A)','16 interleaves, phantom (B)','3 interleaves, phantom (C)','96 interleaves, in vivo (D)'};
    uit = uitable(fig,'Data',Astr,'ColumnName',correction,'RowName',vars,'Units','normalized','Position',[0.05 0.015 0.94 0.153],'ColumnWidth',{158,158},'FontName','Times','Fontsize',11);

end






