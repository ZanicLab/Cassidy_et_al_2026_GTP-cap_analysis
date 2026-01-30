%% Functions (if they have been changed or updated (ie different function name), then they need to be updated below
% line 51 - CometAnalysis_Align_20240709(pathtofile,outputfoldersuffix,filename,windowsize,pospix,timepix)
% line 65 - CometAnalysis_SubtractBGFixedPixel_20240709-noBG(foldername,solpix,latpix)
% line 91 - CometAnalysis_ExponentialFitAverageLinescan_20250722_65nm(path,subfolder,filename,outputfilename,fitstartpix,fitendpix)

%% Before proceeding, the following must be updated:
% line 14 - base folder directory (base_dir)
% line 18 - datasets = {folder within the directory, output folder suffix, points file name}
% line 30 - align_params = {windowsize,pospix,timepix}
% line 32 - bg_subtract_params = {solpix, latpix}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up the base directory
base_dir = "C:\Users\annao\OneDrive - Vanderbilt\Desktop\To_be_analyzed\Interphase_Test_1-25_fit";

% Define datasets to be processed
% datasets = {folder within the directory, output folder suffix, points file name}
datasets = {
    {'100X_20p_LLCPK1-EB1-GFP-H2B-mCherry025\RESULTS', '4s_interphase25', 'Points_file_4sORmore.xlsx'}; %% 0.1615
    {'100X_20p_LLCPK1-EB1-GFP-H2B-mCherry010\RESULTS', '4s_interphase10', 'Points_file_4sORmore.xlsx'}; %% 0.1613
    {'100X_20p_LLCPK1-EB1-GFP-H2B-mCherry012\RESULTS', '4s_interphase12', 'Points_file_4sORmore.xlsx'}; %% 0.161
    {'100x_15p-LLCPK1-EB1GFP-0021\RESULTS', '4s_interphase21', 'Points_file_4sORmore.xlsx'}; %% 0.1501
    {'100x_15p-LLCPK1-EB1GFP-0023\RESULTS', '4s_interphase23', 'Points_file_4sORmore.xlsx'}; %% 0.1491
    {'100x_15p-LLCPK1-EB1GFP-0034\RESULTS', '4s_interphase34', 'Points_file_4sORmore.xlsx'}; %% 0.1467    
    {'100X_15p_LLCPK1-EB1GFP-H2B-mCherry_0058\RESULTS', '4s_interphase58', 'Points_file_4sORmore.xlsx'}; %% 0.1468
    {'100X_15p_LLCPK1-EB1GFP-H2B-mCherry_0059\RESULTS', '4s_interphase59', 'Points_file_4sORmore.xlsx'}; %% 0.1466
    {'100X_20p_LLCPK1-EB1GFP-H2B-mCherry_0016\RESULTS', '4s_interphase16', 'Points_file_4sORmore.xlsx'}; %% 0.1466    
};

% Parameters for alignment and background subtraction
% align_params = {windowsize,pospix,timepix}
% do not change windowsize = 5
%% SDC pixel size 110 nm at 5 frames/s (0.2 s/frame)
%%Burnette lab SDC pixel size 65 nm with variable frame rates
align_params = {5, 65, 0.15}; 
% bg_subtract_params = {solpix, latpix}
bg_subtract_params = {30, 50};

% % Perform alignment for each dataset
% for i = 1:length(datasets)
%     disp(['Processing dataset ' num2str(i) ' of ' num2str(length(datasets))]);  % Debug output
% 
%     data = datasets{i};
%     date = data{1};  % Extract the date
%     sample = data{2};  % Extract the sample name
%     points_file = data{3};  % Extract the points file name
%     
%     % Construct the full file path for the points file
%     points_file_path = fullfile(base_dir, date, points_file);
%     
%     % Check if the file exists at the constructed path
%     if exist(points_file_path, 'file') == 2
%         % If the file exists, perform the alignment
%         disp(['Aligning dataset: ' points_file]);  % Debug output
%         CometAnalysis_Align_20240709(fullfile(base_dir, date), sample, points_file, align_params{:});
%     else
%         % If the file doesn't exist, display an error message
%         disp(['Error: Points file not found: ' points_file_path]);  % Debug output
%     end
% end
% 
% 
% % Perform background subtraction for each dataset
% for i = 1:length(datasets)
%     data = datasets{i};
%     date = data{1};
%     sample = data{2};
%     foldername = fullfile(base_dir, date, ['individualLS_', sample]);
%     CometAnalysis_SubtractBGFixedPixel_20240709(foldername, bg_subtract_params{:});
% end




% Perform exponential fit analysis for each dataset
for i = 1:length(datasets)
    data = datasets{i};
    date = data{1};  % The date of the dataset (e.g., '20211116')
    sample = data{2};  % The sample name (e.g., '5s_Cells')

    % Create the full folder path for the dataset
    this_folder = fullfile(base_dir, date, ['individualLS_', sample]);

    % Get the list of files to process (files with '-Average.dat' extension)
    filesInFolder = dir(fullfile(this_folder, '*-Average.dat'));
    NumofFiles = numel(filesInFolder);

    % Perform exponential fit analysis on each file found in the folder
    for ii = 1:NumofFiles
        close all  % Close any existing figures
        iiname = filesInFolder(ii).name;  % Get the file name (e.g., 'K0119-BGsubstracted-Average.dat')
        
        % Define the output filename based on the sample
        fit_filename = sprintf('Fit_p1-25-NoYOffset_%s', sample);  % Example: 'Fit_p1-25_5s_Cells'
        
        % Call the function to perform the exponential fit analysis
        CometAnalysis_ExponentialFitAverageLinescan_20250722_65nm(...
            this_folder, ...   % Path to the folder containing the file
            '', ...            % Empty string since we're passing the full path already
            iiname, ...        % Name of the current file
            fit_filename, ...  % Output filename
            1, ...             % Parameter for fitting (e.g., starting frame)
            25 ...             % Parameter for fitting (e.g., number of frames)
        );
    end
end
