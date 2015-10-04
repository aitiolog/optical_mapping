%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optical mapping APD analysis script
%
% by Klemen Ziberna and Alexandra S. Mighiu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data

file_name = 'test_data_1.roi';
file_folder = 'roi';

file_name2 = file_name(1:end-4);
%raw_data = importdata([file_folder,file_name]);

raw_data = importdata(file_name);

n_rois = size(raw_data, 2)-1;
t = transpose(raw_data(:,1));

for i=1:n_rois
   raw_signal{i} = transpose(raw_data(:,i+1));   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIGNAL PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%

t_max = max(t);

% Invert signals (i.e. multiply by -1)

for i=1:n_rois
    F_raw{i} = (-1)*raw_signal{i};
    F_raw_max{i} =  max(F_raw{i});
    F_raw_min{i} =  min(F_raw{i});
    F_raw_range{i} = abs(F_raw_max{i}-F_raw_min{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signal smoothing

% Baseline drift correction (using the polynomial fit)

opol = 6; %order of the fitted polynomial (customize!)

for i=1:n_rois
    [p,s,mu] = polyfit(t,F_raw{i},opol);
    f_y{i} = polyval(p,t,[],mu);
    F_BC{i} = F_raw{i} - f_y{i};
end

% Moving average filter

windowSize = 15; %customize!
b = (1/windowSize)*ones(1,windowSize);
a = 1;

for i=1:n_rois
    F_avg{i} = filter(b, a, F_BC{i});
    %F_avg_max{i} =  max(F_avg{i});
    %F_avg_min{i} =  min(F_avg{i});
    F_avg_max{i} =  max(F_avg{i}(2*windowSize:(end-2*windowSize)));
    F_avg_min{i} =  min(F_avg{i}(2*windowSize:(end-2*windowSize))); 
    F_avg_range{i} = abs(F_avg_max{i}-F_avg_min{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%

% First derivative (dF/dt)

for i=1:n_rois
    dF_raw_dt{i} = diff(F_raw{i})./diff(t);
    dF_raw_dt_max{i} = max(dF_raw_dt{i});
    dF_raw_dt_min{i} = min(dF_raw_dt{i});
    dF_raw_dt_range{i} = abs(dF_raw_dt_max{i}-dF_raw_dt_min{i});
end

for i=1:n_rois
    dF_avg_dt{i} = diff(F_avg{i})./diff(t);
    %dF_avg_dt_max{i} = max(dF_avg_dt{i});
    %dF_avg_dt_min{i} = min(dF_avg_dt{i});
    dF_avg_dt_max{i} =  max(dF_avg_dt{i}(2*windowSize:(end-2*windowSize)));
    dF_avg_dt_min{i} =  min(dF_avg_dt{i}(2*windowSize:(end-2*windowSize))); 
    dF_avg_dt_range{i} = abs(dF_avg_dt_max{i}-dF_avg_dt_min{i});
end

tt = t(1:end-1)+diff(t)./2; % Associated time values (central difference)


% Find peaks

F_peak_MinHeight_set = 0.5; % min peak fold increase above baseline
F_peak_MinDistance = 50; % min peak separation (customizable)

for i=1:n_rois
    [F_avg_peaks{i}, F_avg_peaks_ind{i}] = findpeaks(F_avg{i}, ...
        'MinPeakHeight', (F_avg_min{1}+F_avg_range{1}*F_peak_MinHeight_set), ...
        'MinPeakDistance', F_peak_MinDistance);
    N_F_avg_peaks{i} = length(F_avg_peaks{i});
    mean_F_interval{i} = mean(diff(t(F_avg_peaks_ind{i})));
end

% Determine the baseline before the peak

baseline_det_rel_offset = 0.3; %relative (to mean PP) time before the peak
baseline_det_abs_offset = 25; %absolute offset before the peak (ms)

for i=1:n_rois
    for j=1:N_F_avg_peaks{i}   
        
        baseline_start_ind{i}{j}= F_avg_peaks_ind{i}(j)-...
            round(baseline_det_rel_offset*mean_F_interval{i})-baseline_det_abs_offset;
        
        if baseline_start_ind{i}{j} <1
            baseline_start_ind{i}{j} = 1; %in case start border is out of range
        end
        
        baseline_end_ind{i}{j}= F_avg_peaks_ind{i}(j)-baseline_det_abs_offset;
        baseline_mean{i}{j} = mean(F_avg{i}(baseline_start_ind{i}{j}:baseline_end_ind{i}{j}));
        
    end    
end



% Calculate max_peak
% Calculate mean (or RMS)
% Calculate t0 (PP interval)

% Calculate APD10, APD50, APD90
% Calculate max_upstroke (max dF/dt)




%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAPHS
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signal smoothing plots

% signal_fig = figure(1);
% 
% subplot(5,4,1:3);
% %set(gcf,'Visible','off'); %prevent figures to pop up on screen
% plot(t,F_raw{1}, 'b.-'), grid;
% title([file_name2, ' - Optical Mapping - Raw Signal'],'Interpreter', 'none');
% xlim([0,t_max]);
% ylim([F_raw_min{1}-0.1*F_raw_range{1}, F_raw_max{1}+0.1*F_raw_range{1}]);
% xlabel 'Time (ms)', ylabel 'F (AU)';
% 
% subplot(5,4,5:7);
% %set(gcf,'Visible','off'); %prevent figures to pop up on screen
% plot(t,f_y{1}, 'g.-'), grid;
% title([file_name2, ' - Optical Mapping - Baseline Correction'],'Interpreter', 'none');
% xlim([0,t_max]);
% ylim([F_raw_min{1}-0.1*F_raw_range{1}, F_raw_max{1}+0.1*F_raw_range{1}]);
% xlabel 'Time (ms)', ylabel 'F (AU)';
% 
% subplot(5,4,9:11);
% %set(gcf,'Visible','off'); %prevent figures to pop up on screen
% %plot(t,F_avg{1}, 'b.-'), grid;
% plot(t,F_avg{1}, 'b.-', t(F_avg_peaks_ind{1}), F_avg_peaks{1}, 'oc'), grid;
% 
% for f=1:N_F_avg_peaks{1}
% line([t(baseline_start_ind{1}{f}), t(baseline_end_ind{1}{f})],...
%     [baseline_mean{1}{f}, baseline_mean{1}{f}],...
%     'LineWidth', 2, 'Color', [0 0.6 0]) %baseline detection plot
% end
% 
% 
% title([file_name2, ' - Optical Mapping - Baseline Correction + Smoothing'],'Interpreter', 'none');
% xlim([0,t_max]);
% ylim([F_avg_min{1}-0.1*F_avg_range{1}, F_avg_max{1}+0.1*F_avg_range{1}]);
% xlabel 'Time (ms)', ylabel 'F (AU)';
% 
% subplot(5,4,13:15);
% %set(gcf,'Visible','off'); %prevent figures to pop up on screen
% plot(tt,dF_raw_dt{1}, 'r.-'), grid;
% title([file_name2, ' - Derivative - Raw Signal'],'Interpreter', 'none');
% xlim([0,t_max]);
% ylim([dF_raw_dt_min{1}-0.1*dF_raw_dt_range{1}, dF_raw_dt_max{1}+0.1*dF_raw_dt_range{1}]);
% xlabel 'Time (ms)', ylabel 'dF/dt (AU/ms)';
% 
% subplot(5,4,17:19);
% %set(gcf,'Visible','off'); %prevent figures to pop up on screen
% plot(tt,dF_avg_dt{1}, 'r.-'), grid;
% title([file_name2, ' - Derivative - Baseline Correction + Smoothing'],'Interpreter', 'none');
% xlim([0,t_max]);
% ylim([dF_avg_dt_min{1}-0.1*dF_avg_dt_range{1}, dF_avg_dt_max{1}+0.1*dF_avg_dt_range{1}]);
% xlabel 'Time (ms)', ylabel 'dF/dt (AU/ms)';
% 
% % Zoomed in plots
% 
% subplot(5,4,4);
% %set(gcf,'Visible','off'); %prevent figures to pop up on screen
% plot(t,F_raw{1}, 'b.-'), grid;
% %title([file_name2, ' - Optical Mapping - Raw Signal'],'Interpreter', 'none');
% xlim([t(F_avg_peaks_ind{1}(round(N_F_avg_peaks{1}/2)))-100,...
%     t(F_avg_peaks_ind{1}(round(N_F_avg_peaks{1}/2)))+150]);
% ylim([F_raw_min{1}-0.1*F_raw_range{1}, F_raw_max{1}+0.1*F_raw_range{1}]);
% xlabel 'Time (ms)', ylabel 'F (AU)';
% 
% subplot(5,4,8);
% %set(gcf,'Visible','off'); %prevent figures to pop up on screen
% plot(t,f_y{1}, 'g.-'), grid;
% %title([file_name2, ' - Optical Mapping - Baseline Correction'],'Interpreter', 'none');
% xlim([t(F_avg_peaks_ind{1}(round(N_F_avg_peaks{1}/2)))-100,...
%     t(F_avg_peaks_ind{1}(round(N_F_avg_peaks{1}/2)))+150]);
% ylim([F_raw_min{1}-0.1*F_raw_range{1}, F_raw_max{1}+0.1*F_raw_range{1}]);
% xlabel 'Time (ms)', ylabel 'F (AU)';
% 
% subplot(5,4,12);
% %set(gcf,'Visible','off'); %prevent figures to pop up on screen
% %plot(t,F_avg{1}, 'b.-'), grid;
% plot(t,F_avg{1}, 'b.-', t(F_avg_peaks_ind{1}), F_avg_peaks{1}, 'oc'), grid;
% %title([file_name2, ' - Optical Mapping - Baseline Correction + Smoothing'],'Interpreter', 'none');
% xlim([t(F_avg_peaks_ind{1}(round(N_F_avg_peaks{1}/2)))-100,...
%     t(F_avg_peaks_ind{1}(round(N_F_avg_peaks{1}/2)))+150]);
% ylim([F_avg_min{1}-0.1*F_avg_range{1}, F_avg_max{1}+0.1*F_avg_range{1}]);
% xlabel 'Time (ms)', ylabel 'F (AU)';
% 
% subplot(5,4,16);
% %set(gcf,'Visible','off'); %prevent figures to pop up on screen
% plot(tt,dF_raw_dt{1}, 'r.-'), grid;
% %title([file_name2, ' - Derivative - Raw Signal'],'Interpreter', 'none');
% xlim([t(F_avg_peaks_ind{1}(round(N_F_avg_peaks{1}/2)))-100,...
%     t(F_avg_peaks_ind{1}(round(N_F_avg_peaks{1}/2)))+150]);
% ylim([dF_raw_dt_min{1}-0.1*dF_raw_dt_range{1}, dF_raw_dt_max{1}+0.1*dF_raw_dt_range{1}]);
% xlabel 'Time (ms)', ylabel 'dF/dt (AU/ms)';
% 
% subplot(5,4,20);
% %set(gcf,'Visible','off'); %prevent figures to pop up on screen
% plot(tt,dF_avg_dt{1}, 'r.-'), grid;
% %title([file_name2, ' - Derivative - Baseline Correction + Smoothing'],'Interpreter', 'none');
% xlim([t(F_avg_peaks_ind{1}(round(N_F_avg_peaks{1}/2)))-100,...
%     t(F_avg_peaks_ind{1}(round(N_F_avg_peaks{1}/2)))+150]);
% ylim([dF_avg_dt_min{1}-0.1*dF_avg_dt_range{1}, dF_avg_dt_max{1}+0.1*dF_avg_dt_range{1}]);
% xlabel 'Time (ms)', ylabel 'dF/dt (AU/ms)';
% 
% % Save figure
% 
% %saveas(signal_fig, [file_name2, '_signal'], 'fig'); 
% 
% %%print('-djpeg', '-r300', [file_name2,'_signal.jpg']); %simple jpg export
% %r = 300; % pixels per inch
% %set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 3000 6000]/r); % ->resolution: 1500x3000
% %print(gcf, '-djpeg', [file_name2,'_signal.jpg']); %jpg export

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Electrophysiological variables

% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop for generating figures

for i=1:n_rois

signal_fig(i) = figure(i);

subplot(5,4,1:3);
%set(gcf,'Visible','off'); %prevent figures to pop up on screen
plot(t,F_raw{i}, 'b.-'), grid;
title([file_name2, ' - Optical Mapping - Raw Signal'],'Interpreter', 'none');
xlim([0,t_max]);
ylim([F_raw_min{i}-0.1*F_raw_range{i}, F_raw_max{i}+0.1*F_raw_range{i}]);
xlabel 'Time (ms)', ylabel 'F (AU)';

subplot(5,4,5:7);
%set(gcf,'Visible','off'); %prevent figures to pop up on screen
plot(t,f_y{i}, 'g.-'), grid;
title([file_name2, ' - Optical Mapping - Baseline Correction'],'Interpreter', 'none');
xlim([0,t_max]);
ylim([F_raw_min{i}-0.1*F_raw_range{i}, F_raw_max{i}+0.1*F_raw_range{i}]);
xlabel 'Time (ms)', ylabel 'F (AU)';

subplot(5,4,9:11);
%set(gcf,'Visible','off'); %prevent figures to pop up on screen
%plot(t,F_avg{i}, 'b.-'), grid;
plot(t,F_avg{i}, 'b.-', t(F_avg_peaks_ind{i}), F_avg_peaks{i}, 'oc'), grid;

for f=1:N_F_avg_peaks{i}
line([t(baseline_start_ind{i}{f}), t(baseline_end_ind{i}{f})],...
    [baseline_mean{i}{f}, baseline_mean{i}{f}],...
    'LineWidth', 2, 'Color', [0 0.6 0]) %baseline detection plot
end


title([file_name2, ' - Optical Mapping - Baseline Correction + Smoothing'],'Interpreter', 'none');
xlim([0,t_max]);
ylim([F_avg_min{i}-0.1*F_avg_range{i}, F_avg_max{i}+0.1*F_avg_range{i}]);
xlabel 'Time (ms)', ylabel 'F (AU)';

subplot(5,4,13:15);
%set(gcf,'Visible','off'); %prevent figures to pop up on screen
plot(tt,dF_raw_dt{i}, 'r.-'), grid;
title([file_name2, ' - Derivative - Raw Signal'],'Interpreter', 'none');
xlim([0,t_max]);
ylim([dF_raw_dt_min{i}-0.1*dF_raw_dt_range{i}, dF_raw_dt_max{i}+0.1*dF_raw_dt_range{i}]);
xlabel 'Time (ms)', ylabel 'dF/dt (AU/ms)';

subplot(5,4,17:19);
%set(gcf,'Visible','off'); %prevent figures to pop up on screen
plot(tt,dF_avg_dt{i}, 'r.-'), grid;
title([file_name2, ' - Derivative - Baseline Correction + Smoothing'],'Interpreter', 'none');
xlim([0,t_max]);
ylim([dF_avg_dt_min{i}-0.1*dF_avg_dt_range{i}, dF_avg_dt_max{i}+0.1*dF_avg_dt_range{i}]);
xlabel 'Time (ms)', ylabel 'dF/dt (AU/ms)';

% Zoomed in plots

subplot(5,4,4);
%set(gcf,'Visible','off'); %prevent figures to pop up on screen
plot(t,F_raw{i}, 'b.-'), grid;
%title([file_name2, ' - Optical Mapping - Raw Signal'],'Interpreter', 'none');
xlim([t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))-100,...
    t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))+150]);
ylim([F_raw_min{i}-0.1*F_raw_range{i}, F_raw_max{i}+0.1*F_raw_range{i}]);
xlabel 'Time (ms)', ylabel 'F (AU)';

subplot(5,4,8);
%set(gcf,'Visible','off'); %prevent figures to pop up on screen
plot(t,f_y{i}, 'g.-'), grid;
%title([file_name2, ' - Optical Mapping - Baseline Correction'],'Interpreter', 'none');
xlim([t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))-100,...
    t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))+150]);
ylim([F_raw_min{i}-0.1*F_raw_range{i}, F_raw_max{i}+0.1*F_raw_range{i}]);
xlabel 'Time (ms)', ylabel 'F (AU)';

subplot(5,4,12);
%set(gcf,'Visible','off'); %prevent figures to pop up on screen
%plot(t,F_avg{i}, 'b.-'), grid;
plot(t,F_avg{i}, 'b.-', t(F_avg_peaks_ind{i}), F_avg_peaks{i}, 'oc'), grid;
%title([file_name2, ' - Optical Mapping - Baseline Correction + Smoothing'],'Interpreter', 'none');
xlim([t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))-100,...
    t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))+150]);
ylim([F_avg_min{i}-0.1*F_avg_range{i}, F_avg_max{i}+0.1*F_avg_range{i}]);
xlabel 'Time (ms)', ylabel 'F (AU)';

subplot(5,4,16);
%set(gcf,'Visible','off'); %prevent figures to pop up on screen
plot(tt,dF_raw_dt{i}, 'r.-'), grid;
%title([file_name2, ' - Derivative - Raw Signal'],'Interpreter', 'none');
xlim([t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))-100,...
    t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))+150]);
ylim([dF_raw_dt_min{i}-0.1*dF_raw_dt_range{i}, dF_raw_dt_max{i}+0.1*dF_raw_dt_range{i}]);
xlabel 'Time (ms)', ylabel 'dF/dt (AU/ms)';

subplot(5,4,20);
%set(gcf,'Visible','off'); %prevent figures to pop up on screen
plot(tt,dF_avg_dt{i}, 'r.-'), grid;
%title([file_name2, ' - Derivative - Baseline Correction + Smoothing'],'Interpreter', 'none');
xlim([t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))-100,...
    t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))+150]);
ylim([dF_avg_dt_min{i}-0.1*dF_avg_dt_range{i}, dF_avg_dt_max{i}+0.1*dF_avg_dt_range{i}]);
xlabel 'Time (ms)', ylabel 'dF/dt (AU/ms)';

% Save figure

saveas(signal_fig(i), [file_name2, '_', num2str(i), '_signal_roi'], 'fig'); 

%print('-djpeg', '-r300', [file_name2,'_signal.jpg']); %simple jpg export
r = 300; % pixels per inch
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 3000 6000]/r); % ->resolution: 1500x3000
print(gcf, '-djpeg', [file_name2, '_', num2str(i), '_signal.jpg']); %jpg export


% Save the figure in the Matlab format

%saveas(signal_fig, [file_name2], 'fig'); 

% Export plot to jpg

%print('-djpeg', '-r300', [file_name2,'.jpg']);

end