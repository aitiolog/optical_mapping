function APD_analysis_function(file_folder, file_name, output_folder, save_figure, AP_figure)
%APD_analysis_function Analyses action potentials for several ROIs
%   Every file contains several ROIs
%   Outputs APD30, APD50, APD80, AP intervals, AP amplitudes in csv files
%   Optional: saves figures as jpg files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if save_figure argument exists
if (~exist('save_figure', 'var'))
        save_figure = '0';
end

if (~exist('AP_figure', 'var'))
        AP_figure = '0';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optical mapping APD analysis script
%
% by Klemen Ziberna and Alexandra S. Mighiu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% v2 - updates: 
% - better logic to handle APs close to the border
% - added feature to save results in csv file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear;
%close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data

file_name2 = file_name(1:end-4);
raw_data = importdata([file_folder,file_name]);

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
        'MinPeakHeight', (F_avg_min{i}+F_avg_range{i}*F_peak_MinHeight_set), ...
        'MinPeakDistance', F_peak_MinDistance);
        
    % Cutoff at the beginning
    while t(F_avg_peaks_ind{i}(1)) < T_BEGIN_CUTOFF;
        F_avg_peaks{i} = F_avg_peaks{i}(2:end);
        F_avg_peaks_ind{i} = F_avg_peaks_ind{i}(2:end);    
    end
    
    % Cutoff at the end
    while t_max - t(F_avg_peaks_ind{i}(end)) < T_END_CUTOFF;
        F_avg_peaks{i} = F_avg_peaks{i}(1:end-1);
        F_avg_peaks_ind{i} = F_avg_peaks_ind{i}(1:end-1);    
    end
    
    mean_F_interval{i} = mean(diff(t(F_avg_peaks_ind{i})));
    N_F_avg_peaks{i} = length(F_avg_peaks{i});
    
end


% Intervals between the signals

for i=1:n_rois
    F_peak_interval{i} = diff(t(F_avg_peaks_ind{i}(:)));
    N_F_peak_interval{i} = length(F_peak_interval{i});
end



% Determine the baseline before the peak

baseline_det_rel_offset = 0.15; %relative (to mean PP) time before the peak
baseline_det_abs_offset = 25; %absolute offset before the peak

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
% SAVE RESULTS - CSV FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Format:
% - one results file per loaded ROI file
% - one line for one variable (all results horizontaly)
% - naming: file name, ROI number, variable, N measurements, measurements
% - variables: APD30, APD50, APD80, peak_intervals, peak_amplitudes


% Format the output row

for i=1:n_rois  
           
    % APD30
    outputAPD30Row{i} = [{file_name2, i, 'APD30', N_F_avg_peaks{i}},...
        num2cell(APD30{i})];
    outputAPD30RowString{i} = cellfun(@num2str, outputAPD30Row{i},...
        'UniformOutput', false);
    outputAPD30RowSingleString{i} = strjoin(outputAPD30RowString{i}, ',');

    % APD50
    outputAPD50Row{i} = [{file_name2, i, 'APD50', N_F_avg_peaks{i}},...
        num2cell(APD50{i})];
    outputAPD50RowString{i} = cellfun(@num2str, outputAPD50Row{i},...
        'UniformOutput', false);
    outputAPD50RowSingleString{i} = strjoin(outputAPD50RowString{i}, ',');

    % APD80
    outputAPD80Row{i} = [{file_name2, i, 'APD80', N_F_avg_peaks{i}},...
        num2cell(APD80{i})];
    outputAPD80RowString{i} = cellfun(@num2str, outputAPD80Row{i},...
        'UniformOutput', false);
    outputAPD80RowSingleString{i} = strjoin(outputAPD80RowString{i}, ',');

    % AP intervals
    outputF_peak_intervalRow{i} = [{file_name2, i, 'Peak_Intervals',...
        N_F_peak_interval{i}}, num2cell(F_peak_interval{i})];
    outputF_peak_intervalRowString{i} = cellfun(@num2str, ...
        outputF_peak_intervalRow{i}, 'UniformOutput', false);
    outputF_peak_intervalRowSingleString{i} = ...
        strjoin(outputF_peak_intervalRowString{i}, ',');

    % AP amplitude
    outputAP_amplitudeRow{i} = [{file_name2, i, 'AP_amplitude',...
        N_F_avg_peaks{i}}, num2cell(AP_amplitude{i})];
    outputAP_amplitudeRowString{i} = cellfun(@num2str, outputAP_amplitudeRow{i},...
        'UniformOutput', false);
    outputAP_amplitudeRowSingleString{i} = strjoin(outputAP_amplitudeRowString{i}, ',');
        
end


% Merge rows

outputTableRows = [outputAPD30RowSingleString, outputAPD50RowSingleString,...
    outputAPD80RowSingleString, outputF_peak_intervalRowSingleString,...
    outputAP_amplitudeRowSingleString];


%Save file to the csv file

%outputTableFilename = [file_name2, '-table.csv'];
outputTableFilename = [output_folder, file_name2, '-table.csv'];
N_rows = length(outputTableRows);

fileID = fopen(outputTableFilename,'w');
fprintf(fileID,'%s,\n',...
    'file_name,ROI_N,Variable,N_measurements,Measurements>>>>>>'); % header line
for i=1:N_rows
    fprintf(fileID, '%s\n', outputTableRows{i});
end
fclose(fileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY OF RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%

if n_rois == 4;
    % OPTION A: Zoomed in left or right atria (ROI=4)
    disp('Zoomed atria summary results');
    
    % Convert the cell arrays to matrices for easier manipulation
    % Check if the number of APs is the same in all ROIs (-1 if not)

    APD30_size = cellfun(@length, APD30);
    if all(APD30_size == APD30_size(1));
        APD30_mat = vertcat(APD30{:});
    else
        disp('Number of APs in APD30 ROIs do not match - no summary analysis');
        APD30_mat = -1;
    end

    APD50_size = cellfun(@length, APD50);
    if all(APD50_size == APD50_size(1));
        APD50_mat = vertcat(APD50{:});
    else
        disp('Number of APs in APD50 ROIs do not match - no summary analysis');
        APD50_mat = -1;
    end

    APD80_size = cellfun(@length, APD80);
    if all(APD80_size == APD80_size(1));
        APD80_mat = vertcat(APD80{:});
    else
        disp('Number of APs in APD80 ROIs do not match - no summary analysis');
        APD80_mat = -1;
    end

    F_peak_interval_size = cellfun(@length, F_peak_interval);
    if all(F_peak_interval_size == F_peak_interval_size(1));
        F_peak_interval_mat = vertcat(F_peak_interval{:});
    else
        disp('Number of APs in F_peak_interval ROIs do not match - no summary analysis');
        F_peak_interval_mat = -1;
    end

    AP_amplitude_size = cellfun(@length, AP_amplitude);
    if all(AP_amplitude_size == AP_amplitude_size(1));
        AP_amplitude_mat = vertcat(AP_amplitude{:});
    else
        disp('Number of APs in AP_amplitude ROIs do not match - no summary analysis');
        AP_amplitude_mat = -1;
    end

    % Calculations of summary variables

    APD30_avg = mean2(APD30_mat);
    APD30_STD = std2(APD30_mat);
    Loc_APD30_STD = std(mean(APD30_mat,2)); %STD of rows' means

    APD50_avg = mean2(APD50_mat);
    APD50_STD = std2(APD50_mat);
    Loc_APD50_STD = std(mean(APD50_mat,2)); %STD of rows' means

    APD80_avg = mean2(APD80_mat);
    APD80_STD = std2(APD80_mat);
    Loc_APD80_STD = std(mean(APD80_mat,2)); %STD of rows' means

    AP_interval_avg = mean(mean(F_peak_interval_mat));
    AP_interval_STD_avg = std2(mean(F_peak_interval_mat,1));

    AP_amp_CoefVar = mean(std(AP_amplitude_mat, 0, 2)./...
        mean(AP_amplitude_mat,2)); %mean AP amplitude CV of every row

    % Create summary variables csv string

    outputSummaryRow = [{file_name2, 'zoomed atria', ...
        AP_interval_avg, AP_interval_STD_avg, ...
        APD30_avg, APD50_avg, APD80_avg,...
        APD30_STD, APD50_STD, APD80_STD,...
        Loc_APD30_STD, Loc_APD50_STD, Loc_APD80_STD,...
        AP_amp_CoefVar}];
    outputSummaryRowString = cellfun(@num2str, outputSummaryRow,...
        'UniformOutput', false);
    outputSummaryRowSingleString = strjoin(outputSummaryRowString, ',');

    % Save summary variables in csv file

    outputSummaryFilename = [file_name2, '-summary.csv'];

    fileID = fopen(outputSummaryFilename,'w');
    % Header line
    fprintf(fileID,'%s,\n',...
        'file_name,Region,AP_interval_avg,AP_interval_STD_avg,APD30_avg,APD50_avg,APD80_avg,APD30_STD,APD50_STD,APD80_STD,Loc_APD30_STD,Loc_APD50_STD,Loc_APD80_STD,AP_amp_CoefVar');
    % Results line
    fprintf(fileID, '%s\n', outputSummaryRowSingleString);
    fclose(fileID);
    
    
    
    
elseif n_rois == 6;
    % OPTION B: Biatrial analysis (ROI=6)
    % LA - ROI1-2
    % RA - ROI3-4
    % SAN1 - ROI5
    % SAN2 - ROI6
    disp('Biatrial summary results');
    
    % Obtain cell arrays for a particalar ROI region(s)
    
    % LA = ROI1 + ROI2
    LA_APD30 = {APD30{1:2}};
    LA_APD50 = {APD50{1:2}};
    LA_APD80 = {APD80{1:2}};
    LA_F_peak_interval = {F_peak_interval{1:2}};
    LA_AP_amplitude = {AP_amplitude{1:2}};
    
    % RA = ROI3 + ROI4
    RA_APD30 = {APD30{3:4}};
    RA_APD50 = {APD50{3:4}};
    RA_APD80 = {APD80{3:4}};
    RA_F_peak_interval = {F_peak_interval{3:4}};
    RA_AP_amplitude = {AP_amplitude{3:4}};
    
    % SAN1 = ROI5
    SAN1_APD30 = {APD30{5}};
    SAN1_APD50 = {APD50{5}};
    SAN1_APD80 = {APD80{5}};
    SAN1_F_peak_interval = {F_peak_interval{5}};
    SAN1_AP_amplitude = {AP_amplitude{5}};
    
    % SAN2 = ROI6
    SAN2_APD30 = {APD30{6}};
    SAN2_APD50 = {APD50{6}};
    SAN2_APD80 = {APD80{6}};
    SAN2_F_peak_interval = {F_peak_interval{6}};
    SAN2_AP_amplitude = {AP_amplitude{6}};
    
    % Convert the cell arrays to matrices for easier manipulation
    % Need 'APD_ROI_mat_conv.m' function
    
    LA_APD30_mat = APD_ROI_mat_conv(LA_APD30);
    LA_APD50_mat = APD_ROI_mat_conv(LA_APD50);
    LA_APD80_mat = APD_ROI_mat_conv(LA_APD80);
    LA_F_peak_interval_mat = APD_ROI_mat_conv(LA_F_peak_interval);
    LA_AP_amplitude_mat = APD_ROI_mat_conv(LA_AP_amplitude);
    
    RA_APD30_mat = APD_ROI_mat_conv(RA_APD30);
    RA_APD50_mat = APD_ROI_mat_conv(RA_APD50);
    RA_APD80_mat = APD_ROI_mat_conv(RA_APD80);
    RA_F_peak_interval_mat = APD_ROI_mat_conv(RA_F_peak_interval);
    RA_AP_amplitude_mat = APD_ROI_mat_conv(RA_AP_amplitude);
    
    SAN1_APD30_mat = APD_ROI_mat_conv(SAN1_APD30);
    SAN1_APD50_mat = APD_ROI_mat_conv(SAN1_APD50);
    SAN1_APD80_mat = APD_ROI_mat_conv(SAN1_APD80);
    SAN1_F_peak_interval_mat = APD_ROI_mat_conv(SAN1_F_peak_interval);
    SAN1_AP_amplitude_mat = APD_ROI_mat_conv(SAN1_AP_amplitude);
    
    SAN2_APD30_mat = APD_ROI_mat_conv(SAN2_APD30);
    SAN2_APD50_mat = APD_ROI_mat_conv(SAN2_APD50);
    SAN2_APD80_mat = APD_ROI_mat_conv(SAN2_APD80);
    SAN2_F_peak_interval_mat = APD_ROI_mat_conv(SAN2_F_peak_interval);
    SAN2_AP_amplitude_mat = APD_ROI_mat_conv(SAN2_AP_amplitude);
    
    % Calculations of summary variables

    LA_APD30_avg = mean2(LA_APD30_mat);
    LA_APD30_STD = std2(LA_APD30_mat);
    Loc_LA_APD30_STD = std(mean(LA_APD30_mat,2)); %STD of rows' means
    LA_APD50_avg = mean2(LA_APD50_mat);
    LA_APD50_STD = std2(LA_APD50_mat);
    Loc_LA_APD50_STD = std(mean(LA_APD50_mat,2)); %STD of rows' means
    LA_APD80_avg = mean2(LA_APD80_mat);
    LA_APD80_STD = std2(LA_APD80_mat);
    Loc_LA_APD80_STD = std(mean(LA_APD80_mat,2)); %STD of rows' means
    LA_AP_interval_avg = mean(mean(LA_F_peak_interval_mat));
    LA_AP_interval_STD_avg = std2(mean(LA_F_peak_interval_mat,1));
    LA_AP_amp_CoefVar = mean(std(LA_AP_amplitude_mat, 0, 2)./...
        mean(LA_AP_amplitude_mat,2)); %mean AP amplitude CV of every row

    RA_APD30_avg = mean2(RA_APD30_mat);
    RA_APD30_STD = std2(RA_APD30_mat);
    Loc_RA_APD30_STD = std(mean(RA_APD30_mat,2)); %STD of rows' means
    RA_APD50_avg = mean2(RA_APD50_mat);
    RA_APD50_STD = std2(RA_APD50_mat);
    Loc_RA_APD50_STD = std(mean(RA_APD50_mat,2)); %STD of rows' means
    RA_APD80_avg = mean2(RA_APD80_mat);
    RA_APD80_STD = std2(RA_APD80_mat);
    Loc_RA_APD80_STD = std(mean(RA_APD80_mat,2)); %STD of rows' means
    RA_AP_interval_avg = mean(mean(RA_F_peak_interval_mat));
    RA_AP_interval_STD_avg = std2(mean(RA_F_peak_interval_mat,1));
    RA_AP_amp_CoefVar = mean(std(RA_AP_amplitude_mat, 0, 2)./...
        mean(RA_AP_amplitude_mat,2)); %mean AP amplitude CV of every row
    
    SAN1_APD30_avg = mean2(SAN1_APD30_mat);
    SAN1_APD30_STD = std2(SAN1_APD30_mat);
    Loc_SAN1_APD30_STD = std(mean(SAN1_APD30_mat,2)); %STD of rows' means
    SAN1_APD50_avg = mean2(SAN1_APD50_mat);
    SAN1_APD50_STD = std2(SAN1_APD50_mat);
    Loc_SAN1_APD50_STD = std(mean(SAN1_APD50_mat,2)); %STD of rows' means
    SAN1_APD80_avg = mean2(SAN1_APD80_mat);
    SAN1_APD80_STD = std2(SAN1_APD80_mat);
    Loc_SAN1_APD80_STD = std(mean(SAN1_APD80_mat,2)); %STD of rows' means
    SAN1_AP_interval_avg = mean(mean(SAN1_F_peak_interval_mat));
    SAN1_AP_interval_STD_avg = std2(mean(SAN1_F_peak_interval_mat,1));
    SAN1_AP_amp_CoefVar = mean(std(SAN1_AP_amplitude_mat, 0, 2)./...
        mean(SAN1_AP_amplitude_mat,2)); %mean AP amplitude CV of every row
    
    SAN2_APD30_avg = mean2(SAN2_APD30_mat);
    SAN2_APD30_STD = std2(SAN2_APD30_mat);
    Loc_SAN2_APD30_STD = std(mean(SAN2_APD30_mat,2)); %STD of rows' means
    SAN2_APD50_avg = mean2(SAN2_APD50_mat);
    SAN2_APD50_STD = std2(SAN2_APD50_mat);
    Loc_SAN2_APD50_STD = std(mean(SAN2_APD50_mat,2)); %STD of rows' means
    SAN2_APD80_avg = mean2(SAN2_APD80_mat);
    SAN2_APD80_STD = std2(SAN2_APD80_mat);
    Loc_SAN2_APD80_STD = std(mean(SAN2_APD80_mat,2)); %STD of rows' means
    SAN2_AP_interval_avg = mean(mean(SAN2_F_peak_interval_mat));
    SAN2_AP_interval_STD_avg = std2(mean(SAN2_F_peak_interval_mat,1));
    SAN2_AP_amp_CoefVar = mean(std(SAN2_AP_amplitude_mat, 0, 2)./...
        mean(SAN2_AP_amplitude_mat,2)); %mean AP amplitude CV of every row
    
    % Create summary variables csv string

    LA_outputSummaryRow = [{file_name2, 'LA_biatrial', ...
        LA_AP_interval_avg, LA_AP_interval_STD_avg, ...
        LA_APD30_avg, LA_APD50_avg, LA_APD80_avg,...
        LA_APD30_STD, LA_APD50_STD, LA_APD80_STD,...
        Loc_LA_APD30_STD, Loc_LA_APD50_STD, Loc_LA_APD80_STD,...
        LA_AP_amp_CoefVar}];
    LA_outputSummaryRowString = cellfun(@num2str, LA_outputSummaryRow,...
        'UniformOutput', false);
    LA_outputSummaryRowSingleString = strjoin(LA_outputSummaryRowString, ',');
    
    RA_outputSummaryRow = [{file_name2, 'RA_biatrial', ...
        RA_AP_interval_avg, RA_AP_interval_STD_avg, ...
        RA_APD30_avg, RA_APD50_avg, RA_APD80_avg,...
        RA_APD30_STD, RA_APD50_STD, RA_APD80_STD,...
        Loc_RA_APD30_STD, Loc_RA_APD50_STD, Loc_RA_APD80_STD,...
        RA_AP_amp_CoefVar}];
    RA_outputSummaryRowString = cellfun(@num2str, RA_outputSummaryRow,...
        'UniformOutput', false);
    RA_outputSummaryRowSingleString = strjoin(RA_outputSummaryRowString, ',');
    
    SAN1_outputSummaryRow = [{file_name2, 'SAN1_biatrial', ...
        SAN1_AP_interval_avg, SAN1_AP_interval_STD_avg, ...
        SAN1_APD30_avg, SAN1_APD50_avg, SAN1_APD80_avg,...
        SAN1_APD30_STD, SAN1_APD50_STD, SAN1_APD80_STD,...
        Loc_SAN1_APD30_STD, Loc_SAN1_APD50_STD, Loc_SAN1_APD80_STD,...
        SAN1_AP_amp_CoefVar}];
    SAN1_outputSummaryRowString = cellfun(@num2str, SAN1_outputSummaryRow,...
        'UniformOutput', false);
    SAN1_outputSummaryRowSingleString = strjoin(SAN1_outputSummaryRowString, ',');
    
    SAN2_outputSummaryRow = [{file_name2, 'SAN2_biatrial', ...
        SAN2_AP_interval_avg, SAN2_AP_interval_STD_avg, ...
        SAN2_APD30_avg, SAN2_APD50_avg, SAN2_APD80_avg,...
        SAN2_APD30_STD, SAN2_APD50_STD, SAN2_APD80_STD,...
        Loc_SAN2_APD30_STD, Loc_SAN2_APD50_STD, Loc_SAN2_APD80_STD,...
        SAN2_AP_amp_CoefVar}];
    SAN2_outputSummaryRowString = cellfun(@num2str, SAN2_outputSummaryRow,...
        'UniformOutput', false);
    SAN2_outputSummaryRowSingleString = strjoin(SAN2_outputSummaryRowString, ',');
    
    % Save summary variables in csv file

    outputSummaryFilename = [file_name2, '-summary.csv'];

    fileID = fopen(outputSummaryFilename,'w');
    % Header line
    fprintf(fileID,'%s,\n',...
        'file_name,Region,AP_interval_avg,AP_interval_STD_avg,APD30_avg,APD50_avg,APD80_avg,APD30_STD,APD50_STD,APD80_STD,Loc_APD30_STD,Loc_APD50_STD,Loc_APD80_STD,AP_amp_CoefVar');
    % Results line
    fprintf(fileID, '%s\n', LA_outputSummaryRowSingleString);
    fprintf(fileID, '%s\n', RA_outputSummaryRowSingleString);
    fprintf(fileID, '%s\n', SAN1_outputSummaryRowSingleString);
    fprintf(fileID, '%s\n', SAN2_outputSummaryRowSingleString);
    fclose(fileID);
    
    
    
    
else
    disp('N_ROI not recognized - no summary results generated');
end





%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAPHS
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Switch for figure plots

if strcmp('jpg', save_figure);
    
    % Fontsize for plots
    fontSize = 4;


    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Signal smoothing plots

    % Loop for generating figures

    for i=1:n_rois

    signal_fig(i) = figure(i);

    subplot(5,4,1:3);
    set(gcf,'Visible','off'); %prevent figures to pop up on screen
    plot(t,F_raw{i}, 'b.-'), grid;
    set(gca,'FontSize',fontSize);
    title([file_name2, '_ROI', num2str(i), ' - Optical Mapping - Raw Signal'],...
        'Interpreter', 'none');
    xlim([0,t_max]);
    ylim([F_raw_min{i}-0.1*F_raw_range{i}, F_raw_max{i}+0.1*F_raw_range{i}]);
    xlabel 'Time (ms)', ylabel 'F (AU)';

    subplot(5,4,5:7);
    set(gcf,'Visible','off'); %prevent figures to pop up on screen
    plot(t,f_y{i}, 'g.-'), grid;
    set(gca,'FontSize',fontSize);
    title([file_name2, '_ROI', num2str(i), ' - Optical Mapping - Baseline Correction'],...
        'Interpreter', 'none');
    xlim([0,t_max]);
    ylim([F_raw_min{i}-0.1*F_raw_range{i}, F_raw_max{i}+0.1*F_raw_range{i}]);
    xlabel 'Time (ms)', ylabel 'F (AU)';

    subplot(5,4,9:11);
    set(gcf,'Visible','off'); %prevent figures to pop up on screen
    %plot(t,F_avg{i}, 'b.-'), grid;
    plot(t,F_avg{i}, 'b.-', t(F_avg_peaks_ind{i}(1:N_F_avg_peaks{i})), ...
        F_avg_peaks{i}(1:N_F_avg_peaks{i}), 'oc'), grid;
    set(gca,'FontSize',fontSize);
    
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
    set(gcf,'Visible','off'); %prevent figures to pop up on screen
    plot(tt,dF_raw_dt{i}, 'r.-'), grid;
    set(gca,'FontSize',fontSize);
    title([file_name2, '_ROI', num2str(i), ' - Derivative - Raw Signal'],'Interpreter', 'none');
    xlim([0,t_max]);
    ylim([dF_raw_dt_min{i}-0.1*dF_raw_dt_range{i}, dF_raw_dt_max{i}+0.1*dF_raw_dt_range{i}]);
    xlabel 'Time (ms)', ylabel 'dF/dt (AU/ms)';

    subplot(5,4,17:19);
    set(gcf,'Visible','off'); %prevent figures to pop up on screen
    %plot(tt,dF_avg_dt{i}, 'r.-'), grid;
    plot(tt,dF_avg_dt{i}, 'r.-', t(t_act_ind{i}(1:N_F_avg_peaks{i})), ...
        dF_dt_max{i}(1:N_F_avg_peaks{i}), 'oc'), grid;
    set(gca,'FontSize',fontSize);
    title([file_name2, '_ROI', num2str(i), ' - Derivative - Baseline Correction + Smoothing'],...
        'Interpreter', 'none');
    xlim([0,t_max]);
    ylim([dF_avg_dt_min{i}-0.1*dF_avg_dt_range{i}, dF_avg_dt_max{i}+0.1*dF_avg_dt_range{i}]);
    xlabel 'Time (ms)', ylabel 'dF/dt (AU/ms)';

    % Zoomed in plots

    subplot(5,4,4);
    set(gcf,'Visible','off'); %prevent figures to pop up on screen
    plot(t,F_raw{i}, 'b.-'), grid;
    set(gca,'FontSize',fontSize);
    title(['Peak: ', num2str(round(N_F_avg_peaks{i}/2))],'Interpreter', 'none');
    xlim([t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))-100,...
        t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))+150]);
    ylim([F_raw_min{i}-0.1*F_raw_range{i}, F_raw_max{i}+0.1*F_raw_range{i}]);
    xlabel 'Time (ms)', ylabel 'F (AU)';

    subplot(5,4,8);
    set(gcf,'Visible','off'); %prevent figures to pop up on screen
    plot(t,f_y{i}, 'g.-'), grid;
    set(gca,'FontSize',fontSize);
    title(['Peak: ', num2str(round(N_F_avg_peaks{i}/2))],'Interpreter', 'none');
    xlim([t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))-100,...
        t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))+150]);
    ylim([F_raw_min{i}-0.1*F_raw_range{i}, F_raw_max{i}+0.1*F_raw_range{i}]);
    xlabel 'Time (ms)', ylabel 'F (AU)';

    subplot(5,4,12);
    set(gcf,'Visible','off'); %prevent figures to pop up on screen
    %plot(t,F_avg{i}, 'b.-'), grid;
    plot(t,F_avg{i}, 'b.-', t(F_avg_peaks_ind{i}), F_avg_peaks{i}, 'oc'), grid;
    set(gca,'FontSize',fontSize);
    title(['Peak: ', num2str(round(N_F_avg_peaks{i}/2))],'Interpreter', 'none');
    xlim([t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))-100,...
        t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))+150]);
    ylim([F_avg_min{i}-0.1*F_avg_range{i}, F_avg_max{i}+0.1*F_avg_range{i}]);
    xlabel 'Time (ms)', ylabel 'F (AU)';

    subplot(5,4,16);
    set(gcf,'Visible','off'); %prevent figures to pop up on screen
    plot(tt,dF_raw_dt{i}, 'r.-'), grid;
    set(gca,'FontSize',fontSize);
    title(['Peak: ', num2str(round(N_F_avg_peaks{i}/2))],'Interpreter', 'none');
    xlim([t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))-100,...
        t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))+150]);
    ylim([dF_raw_dt_min{i}-0.1*dF_raw_dt_range{i}, dF_raw_dt_max{i}+0.1*dF_raw_dt_range{i}]);
    xlabel 'Time (ms)', ylabel 'dF/dt (AU/ms)';

    subplot(5,4,20);
    set(gcf,'Visible','off'); %prevent figures to pop up on screen
    plot(tt,dF_avg_dt{i}, 'r.-'), grid;
    set(gca,'FontSize',fontSize);
    title(['Peak: ', num2str(round(N_F_avg_peaks{i}/2))],'Interpreter', 'none');
    xlim([t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))-100,...
        t(F_avg_peaks_ind{i}(round(N_F_avg_peaks{i}/2)))+150]);
    ylim([dF_avg_dt_min{i}-0.1*dF_avg_dt_range{i}, dF_avg_dt_max{i}+0.1*dF_avg_dt_range{i}]);
    xlabel 'Time (ms)', ylabel 'dF/dt (AU/ms)';

    % Save figure

    %saveas(signal_fig(i), [file_name2, '_ROI', num2str(i), '_signal'], 'fig'); 

    %print('-djpeg', '-r300', [file_name2,'_signal.jpg']); %simple jpg export
    r = 300; % pixels per inch
    %set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 3000 6000]/r); % ->resolution: 1500x3000
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 1000 2000]/r); % ->resolution: 500x1000
    print(gcf, '-djpeg', [output_folder, file_name2, '_ROI', num2str(i), '_signal.jpg']); %jpg export

    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Electrophysiological variables

    if strcmp('true', AP_figure);

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
            set(gcf,'Visible','off'); %prevent figures to pop up on screen
            plot(t,F_avg{i}, 'b.-'), grid;
            set(gca,'FontSize',fontSize);
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
            set(gcf,'Visible','off'); %prevent figures to pop up on screen
            plot(tt,dF_avg_dt{i}, 'r.-', t(t_act_ind{i}(1:N_F_avg_peaks{i})), ...
            dF_dt_max{i}(1:N_F_avg_peaks{i}), 'oc'), grid;
            set(gca,'FontSize',fontSize);

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
            %set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 900*N_F_avg_peaks{i} 2160]/r); % ->resolution: N_AP*450 x 1080
            set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 300*N_F_avg_peaks{i} 720]/r); % ->resolution: N_AP*150 x 360
            print(gcf, '-djpeg', [output_folder, file_name2, '_ROI', num2str(i), '_AP.jpg']); %jpg export


        end
    else
        disp('Plotting single AP figures: OFF');
    end
    
else
    disp('Plotting figures: OFF');
    
end

end

