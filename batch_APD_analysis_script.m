%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BATCH ANALYSIS OF OPTICAL MAPPING APDs
%
% by Klemen Ziberna and Alexandra S. Mighiu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL VARIABLES

%ROI_FILE_FOLDER = 'ROI_Files/'; %it must contain the last '/' or '\'(Win)
ROI_FILE_FOLDER = '/Users/klemen/Repositories/optical_mapping/ROI_Files/';
ROI_EXTENSION = '*.roi';

%OUTPUT_FOLDER = 'Output/'; %it must contain the last '/' or '\'(Win)
OUTPUT_FOLDER = '/Users/klemen/Repositories/optical_mapping/Output/';

SAVE_FIGURES = 'jpg';
%SAVE_FIGURES = '0'; %switch off figure saving

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a list of files
ROI_files = dir([ROI_FILE_FOLDER,...
    ROI_EXTENSION]); %struct - filename accessed using ROI_files(i).name
N_ROI_files = length(ROI_files);

disp(['Number of files for analysis: ', num2str(N_ROI_files)]);

% Create a list of all filenames (without the extension)
for i=1:N_ROI_files
    ROI_file_name{i} = ROI_files(i).name(1:end-4);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the APD analysis

% Batch processing
tic

for i=1:N_ROI_files
    toc
    disp(['Analysis of the file (',num2str(i),'/',num2str(N_ROI_files),...
        '): ',ROI_files(i).name,...
        ' - (',num2str(round(10*100*i/N_ROI_files)/10), '%)']);
        
    % Run the function that analyses the ROI file
    APD_analysis_function(ROI_FILE_FOLDER, ROI_files(i).name, ...
        OUTPUT_FOLDER, SAVE_FIGURES);
    close all;
    
end
toc

