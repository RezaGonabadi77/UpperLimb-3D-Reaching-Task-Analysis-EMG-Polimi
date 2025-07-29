clear 
close all
clc

%% Load Data

% Define subject number
subject_number = 13; % Change this number as needed

cell_number = 10; % 1---> TRANS_L1, 2----> TRANS_L2 and ...

% Create the file name based on subject number

file_name = ['table_sub', num2str(subject_number), '.mat'];

% Check if the file exists in the current directory
if isfile(file_name)
    % If the file exists, load it
    loaded_data = load(file_name);
    table_var_name = fieldnames(loaded_data);
    loaded_table = loaded_data.(table_var_name{1});
else
    % If the file doesn't exist, create a new table with dimensions 3x9
    loaded_table = cell(2, 9);
    save(file_name, 'loaded_table');
end


% For debugging, display the variable (optional)
disp(['Variable from ', file_name, ' is now loaded or created.']);
disp(loaded_table); % Display the table


% add to path "SAGA_interface" with all subfolders
addpath("SAGA_interface")                               %add path for functions

home = cd;

[FileName1,PathName1] = uigetfile('*.Poly5');         %poly5 data are the ones recorded from the TMSi

cd(PathName1);

data_mot_emg = TMSiSAGA.Poly5.read(FileName1);

fs=data_mot_emg.sample_rate;                 %extract sample frequency

cd(home)

[FileName2,PathName2] = uigetfile('*.mat');            %mat data are the ones extracted from AGREE
cd(PathName2);
load(FileName2)

cd(home)

agree_trajectory = vector;

clear vector;                                                                                       %vector is the name of the variable in the mat file. since it is useless, i cleared it

%% Extract and Define Variables
% Extract trigger and EMG data, and define related variables
% This section isolates trigger and EMG channels, creates a time vector for the data, and lists the muscle names.

trigger=data_mot_emg.samples(21,:);                                                                 %isolate the trigger channel (bipolar electrode)

EMG_mono = data_mot_emg.samples(1:20,:);                                                            %group the signal channels (monopolar electrodes)

time=linspace( 0 , max(size(trigger(1,:)))/fs, max(size(trigger(1,:))) )';                          %time vector for TMSi acquisitions

muscles = {'Brachi', 'Bic-LH', 'Bic-SH', 'Tri', 'Pec', 'AntD', 'MedD', 'PosD', 'Teres', 'Trap'};    %muscle list

%% Notch Filtering and Bipolar EMG Calculation
% Apply notch filter and compute bipolar EMG signals
% This section designs and applies a notch filter to remove 50 Hz noise
% from EMG data and then computes bipolar EMG signals from the filtered monopolar signals.

% Notch filter design
f0 = 50;  % Frequency to notch out
bw = 0.5;  % Bandwidth around the notch frequency
[bn, an] = iirnotch(f0/(fs/2), bw/(fs/2));  % Notch filter coefficients
EMG_mono_f = filtfilt(bn, an, EMG_mono');
% Now you can use smoothed_signal for further analysis or visualization
%%%from monopolar signals bandpass filtered, it is possible to make the
%%%difference between the monopolar electrodes pair used for each muscle
n_channels= 10;                                                                                     % insert the number of channels used, rembember that each muscle activity is acquired with 2 mono channels
for n=1:n_channels
    EMG_bip(:,n) = EMG_mono_f(:,2*(n-1)+1) - EMG_mono_f(:,2*n);
end

%% Removal of the ECG artifacts affecting EMG signals of pectoralis and teres major
% We are using https://github.com/ime-luebeck/ecg-removal so we need to
% change path to access to the already implemented functions 
original_directory = cd;
cd('ecg-removal-main/code');  % Replace with the actual path to your parent directory
% Add the subfolders to the MATLAB path
addpath('filters');
addpath('ecg_utils');
addpath('template_subtraction');
addpath('swt_denoising');

% Change back to the original directory
cd(original_directory);

use_filtfilt = true;

% Mild high-pass filtering to remove baseline and ECG interference 
signalhp05_Pec = butter_filt_stabilized(EMG_bip(:,5), 5, fs, 'high', use_filtfilt, 6);
signalhp20_Pec = butter_filt_stabilized(EMG_bip(:,5), 20, fs, 'high', use_filtfilt, 6);
signalhp20_Teres = butter_filt_stabilized(EMG_bip(:,9), 20, fs, 'high', use_filtfilt, 6);
% We find the rpeaks position
rpeaks_Pec = peak_detection(signalhp05_Pec, fs);

cleaned_ats_Pec = adaptive_template_subtraction(signalhp20_Pec, rpeaks_Pec, fs);
cleaned_ats_Teres = adaptive_template_subtraction(signalhp20_Teres, rpeaks_Pec, fs);
% Replace the original signals with the cleaned signals
EMG_bip(:,5) = cleaned_ats_Pec;
EMG_bip(:,9) = cleaned_ats_Teres;

% Define sampling frequency
fs = 1000; % Example value, replace with your actual sampling frequency

% Define start and end times
startTime = 3; % in seconds
endTime = 8;   % in seconds

% Calculate indices for the time window
startIndex = startTime * fs + 1;
endIndex = endTime * fs;

% Create time vector
timeVector = linspace(startTime, endTime, endIndex - startIndex + 1);

% Plot the original signal
figure
plot(timeVector, signalhp20_Pec(startIndex:endIndex), 'b', 'DisplayName', 'Original');
hold on;

% Plot the cleaned signal
plot(timeVector, cleaned_ats_Pec(startIndex:endIndex), 'r', 'DisplayName', 'Cleaned');

% Add title and labels
title('Overlay of Original and Cleaned EMG Pectoralis','FontSize', 13, 'FontWeight', 'bold');
xlabel('Time (s)','FontSize', 12);
ylabel('Amplitude (mV)','FontSize', 12);

% Add legend and grid
legend('show', 'FontSize', 10);
grid on;

set(gca);
% Adjust figure properties
set(gcf, 'PaperPositionMode', 'auto');

% Save as high-quality PNG
% exportgraphics(gcf, 'high_quality_plot.png', 'Resolution', 300);
%% Convert to Monopolar and Apply Band-Pass Filter
% Convert bipolar EMG to monopolar and apply band-pass filter
% This section converts bipolar EMG signals to monopolar for band-pass filtering,
% and then recomputes the bipolar signals from the filtered data.

n_channels = 10;  % Number of bipolar channels

% Assuming num_samples is the number of samples in your data
num_samples = size(EMG_bip, 1);

% Initialize EMG_mono_f_new with zeros
EMG_mono_f_new = zeros(num_samples, 2 * n_channels);

for n = 1:n_channels-1
    EMG_mono_f_new(:, 2*(n-1)+1) = (EMG_bip(:, n) + EMG_bip(:, n+1)) / 2;  % First monopolar channel
    EMG_mono_f_new(:, 2*n) = -(EMG_bip(:, n) - EMG_bip(:, n+1)) / 2;        % Second monopolar channel
end

% Special case for the last channel
EMG_mono_f_new(:, 2*n_channels-1) = EMG_bip(:, n_channels);  % Last monopolar channel
EMG_mono_f_new(:, 2*n_channels) = -EMG_bip(:, n_channels);    % Last monopolar channel

fcutlow=20;                                             %[Hz], lower limit for bandpass filter
fcuthigh=400;                                                                                       %[Hz], higher limit for bandpass filter
[B,A]=butter(5,[fcutlow,fcuthigh]/(fs/2),'bandpass');                                               %coefficients extractions for the bandpass filter
EMG_mono_f = filtfilt(B,A,EMG_mono_f_new);   

n_channels= 10;                                                                                     % insert the number of channels used, rembember that each muscle activity is acquired with 2 mono channels
for n=1:n_channels
    EMG_bip(:,n) = EMG_mono_f(:,2*(n-1)+1) - EMG_mono_f(:,2*n);
end
%% Trigger

[pks,locs]=findpeaks(trigger,'MinPeakDistance',300,'MinPeakHeight',1e4);                            %find the trigger peaks
                                                                                                    %pks is the value, locs the position in the array
% figure
% plot(time,trigger);
% hold on;
% plot(time(locs),pks,'o');

%% Signal Processing Steps
% Rectify and envelope EMG signals, then apply smoothing
% This section performs rectification, low-pass filtering to create the envelope, and applies a moving average for further smoothing of the EMG signals.

%Rectification
EMG_abs = abs(EMG_bip);

%Envelope
%low pass 10 Hz - 5 order
[B1,A1] = butter(5,10/(fs/2),'low');
EMG_en1 = filtfilt(B1,A1,EMG_abs);

% Apply moving Average for smoothing the signal
window_size = 100;  % Set the window size for the moving average
type = 'exponential';  % Choose the type of moving average ('simple', 'weighted', 'exponential')
% Ensure the moving_average function is in your MATLAB path
EMG_en2 = moving_average(EMG_en1, window_size, type);

% Create subplots
figure;

% Plot EMG_abs
subplot(2, 1, 1);
plot(time, EMG_en1(:,6));
title('EMG Signal after Low-pass : Anterior Deltoid','FontSize', 12);
xlabel('Time (s)','FontSize', 10);
ylabel('Amplitude (mV)','FontSize', 10);
grid on;

% Plot EMG_en2
subplot(2, 1, 2);
plot(time, EMG_en2(:,6));
title('EMG Signal after Exponential MA : Anterior Deltoid','FontSize', 12);
xlabel('Time (s)','FontSize', 10);
ylabel('Amplitude (mV)','FontSize', 10);
grid on;
set(gca);
% Adjust figure properties
set(gcf, 'PaperPositionMode', 'auto');

% Save as high-quality PNG
% exportgraphics(gcf, 'high_quality_plot3.png', 'Resolution', 300);
%% kinematics data from AGREE
% Extract and synchronize kinematics data from AGREE
% This section finds the trigger point in the kinematics data,
% normalizes the time vector, and creates a common time vector for synchronization with EMG data.

agree_tr = find(agree_trajectory(:,1)==9,1,'last');
time_ag = agree_trajectory(agree_tr:end,2)/1000 - agree_trajectory(agree_tr,2)/1000;                %time from trigger to end of trajectory mod (just before weight comand) (beginning set to zero)

time_tr = time(locs(end):end)-time(locs(end));                                                      %EMG time from 2nd trigger to the end (beginning set to zero)

                                                                                   %link x axis of all subplots
if size(time_ag,1) < size(time_tr,1)                                                                %create a common time vector (the smallest one)
    time_comm = time_ag;
else
    time_comm = time_tr;
end

%% Identification of on-sets and off-sets of the movements
% Identify movement onsets and offsets using kinematic data
% This section creates a matrix with shared time instants and joint angles,
% identifies movement onsets and offsets, and plots the flexion/extension of elbow movement with markers for these events.

% I create a matrix containing shared time instants and the relevant angles
% of the interested joints
joints_displ = [time_comm, rad2deg(agree_trajectory(agree_tr:(size(time_comm,1)+agree_tr-1),3)),...
    rad2deg(agree_trajectory(agree_tr:(size(time_comm,1)+agree_tr-1),4)), rad2deg(agree_trajectory(agree_tr:(size(time_comm,1)+agree_tr-1),6))];
[onset, offset] = movement_identification(joints_displ,fs,0.01,5);

figure
plot(time_comm,joints_displ(:,4),"b")
hold on
% Add vertical lines at offset and onset
xline(time_comm(offset(1:end)), '-.r', 'LineWidth', 1); % Red vertical line at offset
xline(time_comm(onset(1:end)), '-.g', 'LineWidth', 1);  % Red vertical line at onset
% Add title and axis labels
title('Flexion/Extension of Elbow Movement','FontSize', 12)
xlabel('Time (s)','FontSize', 10)
ylabel('Angle (Degree)','FontSize', 10)
set(gca);
% Adjust figure properties
set(gcf, 'PaperPositionMode', 'auto');

% Save as high-quality PNG
% exportgraphics(gcf, 'high_quality_plot2.png', 'Resolution', 300);

%we have detected onset and offset points :))

%% Extract EMG Signals During Identified Activities
% Extract EMG signals during movement activity windows
% This section extracts EMG data within the identified onset and offset time windows,
% synchronizes it, and stores the results in a table for further analysis.

EMG_patch_with_int = cell(length(onset),size(EMG_en2,2));

for i=1:size(EMG_en2,2)
    for j=1:length(onset)
        if (j < 9)
            EMG_start = onset(j)+locs(end)-2500;
            EMG_stop  = offset(j)+locs(end)+2500;
            EMG_patch_with_int{j,i} = EMG_en2(EMG_start:EMG_stop,i);
        end
    end
end

row = ceil(cell_number / 9);
col = mod(cell_number - 1, 9) + 1;
loaded_table{row, col} = EMG_patch_with_int;
save(file_name, 'loaded_table');






