%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optical mapping APD analysis script
%
% by Klemen Ziberna and Alexandra S. Mighiu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% v2: added feature to save results in csv file
%
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
t = transpose(raw_data(:,1)); %time

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

opol = 3; %order of the fitted polynomial (customize!)

for i=1:n_rois
    [p,s,mu] = polyfit(t,F_raw{i},opol);
    f_y{i} = polyval(p,t,[],mu);
    F_BC{i} = F_raw{i} - f_y{i};
end

% Filter

% Option A: Moving average (simple, but bad option)

% Note: some of the functions bellow are dependent upon windowSize

% windowSize = 15; %customize!
% b = (1/windowSize)*ones(1,windowSize);
% a = 1;
% 
% for i=1:n_rois
%     F_avg{i} = filter(b, a, F_BC{i});
%     %F_avg_max{i} =  max(F_avg{i});
%     %F_avg_min{i} =  min(F_avg{i});
%
%     %To avoid dips at the edges when using moving average
%     F_avg_max{i} =  max(F_avg{i}(2*windowSize:(end-2*windowSize))); 
%     F_avg_min{i} =  min(F_avg{i}(2*windowSize:(end-2*windowSize))); 
%     F_avg_range{i} = abs(F_avg_max{i}-F_avg_min{i});
% end


% Option B: IIR filter averaging (Butteworth)
% better approach
% Recommended values: 5th order, 0.6
% which corresponds to 100Hz filter using 1kHz sampling rate

% Formula to calculate the value 'a'
% a = (2*pi*Frequency_cutoff/sampling_frequency)


[b,a] = butter(6,0.4);
freqz(b,a);

for i=1:n_rois
    F_avg{i} = filtfilt(b, a, F_BC{i});
    F_avg_max{i} =  max(F_avg{i});
    F_avg_min{i} =  min(F_avg{i});
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
    dF_avg_dt_max{i} = max(dF_avg_dt{i});
    dF_avg_dt_min{i} = min(dF_avg_dt{i});
    
    % Use it with moving average
    %dF_avg_dt_max{i} =  max(dF_avg_dt{i}(2*windowSize:(end-2*windowSize)));
    %dF_avg_dt_min{i} =  min(dF_avg_dt{i}(2*windowSize:(end-2*windowSize))); 
    dF_avg_dt_range{i} = abs(dF_avg_dt_max{i}-dF_avg_dt_min{i});
end

tt = t(1:end-1)+diff(t)./2; % Associated time values (central difference)


% Find peaks for analysis

F_peak_MinHeight_set = 0.7; % min peak fold increase above baseline
F_peak_MinDistance = 70; % min peak separation (customizable)

T_BEGIN_CUTOFF = 50; % time cutoff for the first AP to be analysed (default 50)
T_END_CUTOFF = 100; % time cutoff for the last AP to be analysed (default 100)

for i=1:n_rois
    [F_avg_peaks{i}, F_avg_peaks_ind{i}] = findpeaks(F_avg{i}, ...
        'MinPeakHeight', (F_avg_min{1}+F_avg_range{1}*F_peak_MinHeight_set), ...
        'MinPeakDistance', F_peak_MinDistance);
    mean_F_interval{i} = mean(diff(t(F_avg_peaks_ind{i})));
    N_F_avg_peaks{i} = length(F_avg_peaks{i});
    
    if t(F_avg_peaks_ind{i}(N_F_avg_peaks{i}))-t_max < T_END_CUTOFF; 
       N_F_avg_peaks{i} = N_F_avg_peaks{i} -1; %omit the last AP if too close to the end
    end
    
    %%%%% TO-DO
    % - add another border condition if the AP is too close to the beginning
%     
%     % while loop (t < T_BEGIN_CUTOFF)
%       % A = A(2:end) and N = N-1
%
%     if t(F_avg_peaks_ind{i}(1)) < T_BEGIN_CUTOFF;
%        % code that will do something
%        % for loop ind2 -> ind1
%        % clever way to delete the first value from the vector
%        
%     end
    
end

% Intervals between the signals

% F_interval{i} = diff(t(F_avg_peaks_ind{i}));



% Determine the baseline before the peak

baseline_det_rel_offset = 0.15; %relative (to mean PP) time before the peak
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


% Max dF/dt and activation time

AP_border_left = 50; % n of points before every peak for max dF/dt detection!

for i=1:n_rois
    for j=1:N_F_avg_peaks{i}
        [dF_dt_max{i}(j), dF_dt_max_rel_ind{i}(j)] = ...
            max(dF_avg_dt{i}(F_avg_peaks_ind{i}(j)-AP_border_left:F_avg_peaks_ind{i}(j)));
        t_act_ind{i}(j) = (F_avg_peaks_ind{i}(j) - AP_border_left...
            + dF_dt_max_rel_ind{i}(j)); %absolute time index
        t_act{i}(j) = t(t_act_ind{i}(j)); %absolute activation times
    end          
end


% Rise time

for i=1:n_rois
    for j=1:N_F_avg_peaks{i}
        AP_rise_time{i}(j) = t(F_avg_peaks_ind{i}(j)) - t_act{i}(j);
    end          
end

% APD calculations

for i=1:n_rois
    for j=1:N_F_avg_peaks{i}
        
        F_max{i}(j) = F_avg_peaks{i}(j);
        AP_amplitude{i}(j) = F_max{i}(j) - baseline_mean{i}{j};
        
        F_AP_30{i}(j) = F_max{i}(j) - 0.3*AP_amplitude{i}(j);
        APD30_rel_ind{i}(j) = ...
            find(F_avg{i}(F_avg_peaks_ind{i}(j):end)<F_AP_30{i}(j),1); %relative index to peak
        APD30_ind{i}(j) = F_avg_peaks_ind{i}(j) + APD30_rel_ind{i}(j); %absolute index
        APD30{i}(j) = t(APD30_ind{i}(j)) - t_act{i}(j); %absolute APD time
        
        F_AP_50{i}(j) = F_max{i}(j) - 0.5*AP_amplitude{i}(j);
        APD50_rel_ind{i}(j) = ...
            find(F_avg{i}(F_avg_peaks_ind{i}(j):end)<F_AP_50{i}(j),1); %relative index to peak
        APD50_ind{i}(j) = F_avg_peaks_ind{i}(j) + APD50_rel_ind{i}(j); %absolute index
        APD50{i}(j) = t(APD50_ind{i}(j)) - t_act{i}(j); %absolute APD time
        
        F_AP_80{i}(j) = F_max{i}(j) - 0.8*AP_amplitude{i}(j);
        APD80_rel_ind{i}(j) = ...
            find(F_avg{i}(F_avg_peaks_ind{i}(j):end)<F_AP_80{i}(j),1); %relative index to peak
        APD80_ind{i}(j) = F_avg_peaks_ind{i}(j) + APD80_rel_ind{i}(j); %absolute index
        APD80{i}(j) = t(APD80_ind{i}(j)) - t_act{i}(j); %absolute APD time
        
    end          
end



%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE RESULTS IN CSV FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%


% TO-DO
% variables:
% - APD30, APD50, APD80
% - AP amplitude
% - interval between APs
%
% Format:
% - one results file per loaded ROI file
% - one line for one variable (all results horizontaly)
% - naming: file name, ROI number, variable, N measurements, measurements



%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAPHS
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signal smoothing plots

% Loop for generating figures

for i=1:n_rois

signal_fig(i) = figure(i);

subplot(5,4,1:3);
%set(gcf,'Visible','off'); %prevent figures to pop up on screen
plot(t,F_raw{i}, 'b.-'), grid;
title([file_name2, '_ROI', num2str(i), ' - Optical Mapping - Raw Signal'],...
    'Interpreter', 'none');
xlim([0,t_max]);
ylim([F_raw_min{i}-0.1*F_raw_range{i}, F_raw_max{i}+0.1*F_raw_range{i}]);
xlabel 'Time (ms)', ylabel 'F (AU)';

subplot(5,4,5:7);
%set(gcf,'Visible','off'); %prevent figures to pop up on screen
plot(t,f_y{i}, 'g.-'), grid;
title([file_name2, '_ROI', num2str(i), ' - Optical Mapping - Baseline Correction'],...
    'Interpreter', 'none');
xlim([0,t_max]);
ylim([F_raw_min{i}-0.1*F_raw_range{i}, F_raw_max{i}+0.1*F_raw_range{i}]);
xlabel 'Time (ms)', ylabel 'F (AU)';

subplot(5,4,9:11);
%set(gcf,'Visible','off'); %prevent figures to pop up on screen
%plot(t,F_avg{i}, 'b.-'), grid;
plot(t,F_avg{i}, 'b.-', t(F_avg_peaks_ind{i}(1:N_F_avg_peaks{i})), ...
    F_avg_peaks{i}(1:N_F_avg_peaks{i}), 'oc'), grid;

for f=1:N_F_avg_peaks{i}
line([t(baseline_start_ind{i}{f}), t(baseline_end_ind{i}{f})],...
    [baseline_mean{i}{f}, baseline_mean{i}{f}],...
    'LineWidth', 2, 'Color', [0 0.6 0]) %baseline detection plot
end

title([file_name2, '_ROI', num2str(i), ' - Optical Mapping - Baseline Correction + Smoothing'],...
    'Interpreter', 'none');
xlim([0,t_max]);
ylim([F_avg_min{i}-0.1*F_avg_range{i}, F_avg_max{i}+0.1*F_avg_range{i}]);
xlabel 'Time (ms)', ylabel 'F (AU)';

subplot(5,4,13:15);
%set(gcf,'Visible','off'); %prevent figures to pop up on screen
plot(tt,dF_raw_dt{i}, 'r.-'), grid;
title([file_name2, '_ROI', num2str(i), ' - Derivative - Raw Signal'],'Interpreter', 'none');
xlim([0,t_max]);
ylim([dF_raw_dt_min{i}-0.1*dF_raw_dt_range{i}, dF_raw_dt_max{i}+0.1*dF_raw_dt_range{i}]);
xlabel 'Time (ms)', ylabel 'dF/dt (AU/ms)';

subplot(5,4,17:19);
%set(gcf,'Visible','off'); %prevent figures to pop up on screen
%plot(tt,dF_avg_dt{i}, 'r.-'), grid;
plot(tt,dF_avg_dt{i}, 'r.-', t(t_act_ind{i}(1:N_F_avg_peaks{i})), ...
    dF_dt_max{i}(1:N_F_avg_peaks{i}), 'oc'), grid;
title([file_name2, '_ROI', num2str(i), ' - Derivative - Baseline Correction + Smoothing'],...
    'Interpreter', 'none');
xlim([0,t_max]);
ylim([dF_avg_dt_min{i}-0.1*dF_avg_dt_range{i}, dF_avg_dt_max{i}+0.1*dF_avg_dt_range{i}]);
xlabel 'Time (ms)', ylabel 'dF/dt (AU/ms)';

% Zoomed in plots

subplot(5,4,4);
%set(gcf,'Visible','off'); %prevent figures to pop up on screen
plot(t,F_raw{i}, 'b.-'), grid;
title(['Peak: ', num2str(round(N_F_avg_peaks{i}/2))],'Interpreter', 'none');
xlim([t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))-100,...
    t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))+150]);
ylim([F_raw_min{i}-0.1*F_raw_range{i}, F_raw_max{i}+0.1*F_raw_range{i}]);
xlabel 'Time (ms)', ylabel 'F (AU)';

subplot(5,4,8);
%set(gcf,'Visible','off'); %prevent figures to pop up on screen
plot(t,f_y{i}, 'g.-'), grid;
title(['Peak: ', num2str(round(N_F_avg_peaks{i}/2))],'Interpreter', 'none');
xlim([t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))-100,...
    t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))+150]);
ylim([F_raw_min{i}-0.1*F_raw_range{i}, F_raw_max{i}+0.1*F_raw_range{i}]);
xlabel 'Time (ms)', ylabel 'F (AU)';

subplot(5,4,12);
%set(gcf,'Visible','off'); %prevent figures to pop up on screen
%plot(t,F_avg{i}, 'b.-'), grid;
plot(t,F_avg{i}, 'b.-', t(F_avg_peaks_ind{i}), F_avg_peaks{i}, 'oc'), grid;
title(['Peak: ', num2str(round(N_F_avg_peaks{i}/2))],'Interpreter', 'none');
xlim([t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))-100,...
    t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))+150]);
ylim([F_avg_min{i}-0.1*F_avg_range{i}, F_avg_max{i}+0.1*F_avg_range{i}]);
xlabel 'Time (ms)', ylabel 'F (AU)';

subplot(5,4,16);
%set(gcf,'Visible','off'); %prevent figures to pop up on screen
plot(tt,dF_raw_dt{i}, 'r.-'), grid;
title(['Peak: ', num2str(round(N_F_avg_peaks{i}/2))],'Interpreter', 'none');
xlim([t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))-100,...
    t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))+150]);
ylim([dF_raw_dt_min{i}-0.1*dF_raw_dt_range{i}, dF_raw_dt_max{i}+0.1*dF_raw_dt_range{i}]);
xlabel 'Time (ms)', ylabel 'dF/dt (AU/ms)';

subplot(5,4,20);
%set(gcf,'Visible','off'); %prevent figures to pop up on screen
plot(tt,dF_avg_dt{i}, 'r.-'), grid;
title(['Peak: ', num2str(round(N_F_avg_peaks{i}/2))],'Interpreter', 'none');
xlim([t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))-100,...
    t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))+150]);
ylim([dF_avg_dt_min{i}-0.1*dF_avg_dt_range{i}, dF_avg_dt_max{i}+0.1*dF_avg_dt_range{i}]);
xlabel 'Time (ms)', ylabel 'dF/dt (AU/ms)';

% Save figure

%saveas(signal_fig(i), [file_name2, '_ROI', num2str(i), '_signal'], 'fig'); 

%print('-djpeg', '-r300', [file_name2,'_signal.jpg']); %simple jpg export
r = 300; % pixels per inch
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 3000 6000]/r); % ->resolution: 1500x3000
%print(gcf, '-djpeg', [file_name2, '_ROI', num2str(i), '_signal.jpg']); %jpg export

end



%%%%%%%%%%%%%%%%%%%%%%%%%%
% Electrophysiological variables

% Signal action potential plots

%Vertical plots
% for i=1:n_rois
%    AP_fig(i) = figure(n_rois+i);
%    
%    for f=1:N_F_avg_peaks{1}
%     
%     subplot(N_F_avg_peaks{i},1,f);
%     %set(gcf,'Visible','off'); %prevent figures to pop up on screen
%     plot(t,F_avg{i}, 'b.-'), grid;
%     %plot(t,F_avg{i}, 'b.-', t(F_avg_peaks_ind{i}), F_avg_peaks{i}, 'oc'),
%     %grid; plot the circle on top of the peak
%     
%         for f_bc=1:N_F_avg_peaks{i}
%             line([t(baseline_start_ind{i}{f_bc}), t(baseline_end_ind{i}{f_bc})],...
%             [baseline_mean{i}{f_bc}, baseline_mean{i}{f_bc}],...
%             'LineWidth', 2, 'Color', [0 0.6 0]) %baseline detection plot
%         end
%     
%     title([file_name2, '_ROI', num2str(i), ' - AP number: ', ...
%         num2str(f)],'Interpreter', 'none');
%     xlim([t(F_avg_peaks_ind{i}(f))-50,...
%         t(F_avg_peaks_ind{i}(f))+100]);
%     ylim([F_avg_min{i}-0.1*F_avg_range{i}, F_avg_max{i}+0.1*F_avg_range{i}]);
%     xlabel 'Time (ms)', ylabel 'F (AU)';
%     
%     end 
%     
% end


%Horizontal plots

for i=1:n_rois
   AP_fig(i) = figure(n_rois+i);
   
   for f=1:N_F_avg_peaks{i}
    
    % Main AP plots
       
    subplot(2,N_F_avg_peaks{i},f);
    %set(gcf,'Visible','off'); %prevent figures to pop up on screen
    plot(t,F_avg{i}, 'b.-'), grid;
    %plot(t,F_avg{i}, 'b.-', t(F_avg_peaks_ind{i}), F_avg_peaks{i}, 'oc'),
    %grid; plot the circle on top of the peak
    
    %Plot different lines to the plot
    
        for f_lines=1:N_F_avg_peaks{i}
            
            line([t(baseline_start_ind{i}{f_lines}), t(baseline_end_ind{i}{f_lines})],...
            [baseline_mean{i}{f_lines}, baseline_mean{i}{f_lines}],...
            'LineWidth', 2, 'Color', [0 0.6 0]); %baseline detection plot
        
            line([t_act{i}(f_lines), t_act{i}(f_lines)],...
            [baseline_mean{i}{f_lines}, F_max{i}(f_lines)],...
            'LineWidth', 0.75, 'Color', [0 0 0],...
            'LineStyle', '--'); %Activation time vertical line
        
            line([t_act{i}(f_lines), t(F_avg_peaks_ind{i}(f_lines))],...
            [F_max{i}(f_lines), F_max{i}(f_lines)],...
            'LineWidth', 2, 'Color', [1 0 0]); %Rise time line
        
            line([t_act{i}(f_lines), t(APD30_ind{i}(f_lines))],...
            [F_AP_30{i}(f_lines), F_AP_30{i}(f_lines)],...
            'LineWidth', 2, 'Color', [1 0 0]); %APD30 line
        
            line([t_act{i}(f_lines), t(APD50_ind{i}(f_lines))],...
            [F_AP_50{i}(f_lines), F_AP_50{i}(f_lines)],...
            'LineWidth', 2, 'Color', [1 0 0]); %APD50 line
        
            line([t_act{i}(f_lines), t(APD80_ind{i}(f_lines))],...
            [F_AP_80{i}(f_lines), F_AP_80{i}(f_lines)],...
            'LineWidth', 2, 'Color', [1 0 0]); %APD80 line
        
        end
    
    title([file_name2, '_ROI', num2str(i), ' - AP', ...
        num2str(f)],'Interpreter', 'none');
    xlim([t(F_avg_peaks_ind{i}(f))-100,...
        t(F_avg_peaks_ind{i}(f))+150]);
    ylim([F_avg_min{i}-0.1*F_avg_range{i}, F_avg_max{i}+0.1*F_avg_range{i}]);
    xlabel 'Time (ms)', ylabel 'F (AU)';
    
    
    
    % First derivative plot
    
    subplot(2,N_F_avg_peaks{i},N_F_avg_peaks{i}+f);
    %set(gcf,'Visible','off'); %prevent figures to pop up on screen
    plot(tt,dF_avg_dt{i}, 'r.-', t(t_act_ind{i}(1:N_F_avg_peaks{i})), ...
    dF_dt_max{i}(1:N_F_avg_peaks{i}), 'oc'), grid;
    
    for f_lines=1:N_F_avg_peaks{i}
        line([t_act{i}(f_lines), t_act{i}(f_lines)],...
            [dF_avg_dt_min{i}-0.1*dF_avg_dt_range{i}, ...
            dF_avg_dt_max{i}+0.1*dF_avg_dt_range{i}],...
            'LineWidth', 0.75, 'Color', [0 0 0],...
            'LineStyle', '--'); %Activation time vertical line
    end

    title([file_name2, '_ROI', num2str(i), ' - AP', ...
        num2str(f)],'Interpreter', 'none');
    xlim([t(F_avg_peaks_ind{i}(f))-100,...
        t(F_avg_peaks_ind{i}(f))+150]);
    ylim([dF_avg_dt_min{i}-0.1*dF_avg_dt_range{i}, dF_avg_dt_max{i}+0.1*dF_avg_dt_range{i}]);
    xlabel 'Time (ms)', ylabel 'dF/dt (AU/ms)';
    
    
   end 
    
   % Save figure

    %saveas(signal_fig(i), [file_name2, '_ROI', num2str(i), '_AP'], 'fig'); 

    %print('-djpeg', '-r300', [file_name2,'_signal.jpg']); %simple jpg export
    r = 300; % pixels per inch
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 900*N_F_avg_peaks{i} 2160]/r); % ->resolution: N_AP*450 x 1080
    %print(gcf, '-djpeg', [file_name2, '_ROI', num2str(i), '_AP.jpg']); %jpg export

   
end
