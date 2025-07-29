clc
clear
close all

%% Select the subject

sub_number = 3; % You can change this number to whatever you want

file_name = sprintf('table_sub%d.mat', sub_number);
load(file_name);
EMG_patches=loaded_table;

%% Time Normalization
EMG_t_norm = cell(2,9);


for i = 1:size(EMG_patches,1) % for all the modalities
    for j=1:size(EMG_patches,2) % for all the tasks
        for k=2:8 % for all the repetitions
            for y=1:10 % for all the muscles
                N = length(EMG_patches{i,j}{k,y});  % Numero di campioni originali
                t_original = linspace(1, N, N);  % Vettore delle ascisse originali
                t_N = linspace(1,N,7500);
                %EMG_t_norm{i,j}{k-1,y} = interp1(EMG_patches{i,j}{k,y}, t_N, 'spline');
                EMG_t_norm{i,j}{k-1,y} = max(interp1(t_original, EMG_patches{i,j}{k,y}, t_N, 'spline'), 0); %interpolation and clipping operation

            end
        end
    end
end

%% Amplitude Normalization

% It has to be done with the data acquired in the transparent mode
% It has to be done separately for each subject
medians = [];
for i = 1:size(EMG_t_norm,2) % for each task
    for j=1:10 % for each muscle
        peak_values = [];
        for k=1:7 % for each repetition
            peak_values(k) = max(EMG_t_norm{1,i}{k,j});
        end
        medians(i,j) = median(peak_values);
    end
end

%% Compute the maximum for each muscle

max_values = [];
for i=1:10
    max_values(i)= max(medians(:,i));
end

%% NORMALIZATION OF ALL THE DATA OF THE SUBJECT
% We divide the value of each rep of a muscle for the maximum computed for
% that muscle

EMG_amp_norm=cell(2,9);
for i = 1:size(EMG_t_norm,1) % for all the modalities
    for j=1:size(EMG_t_norm,2) % for all the tasks
        for k=1:7 % for all the repetitions
            for y=1:10 % for all the muscles
                EMG_amp_norm{i,j}{k,y} = EMG_t_norm{i,j}{k,y}/max_values(y);
            end
        end
    end
end
%% Phasic extracion
EMG_phasic = cell(2,9);
EMG_phasic_OK = cell(2,9);
%count = 0;

for i = 1:size(EMG_amp_norm,1) % for all the modalities
    for j=1:size(EMG_amp_norm,2) % for all the tasks
        for k=1:7 % for all the repetitions
            for y=1:10 % for all the muscles
              mean_on = mean(EMG_amp_norm{i,j}{k,y}(1200:1400));
              mean_off = mean(EMG_amp_norm{i,j}{k,y}(5900:6100));
              ramp = linspace(mean_on, mean_off, 2700);
              EMG_phasic{i,j}{k,y} = EMG_amp_norm{i,j}{k,y}(2301:5000) - ramp;

                for n = 1:length(EMG_phasic{i,j}{k,y})
                    if EMG_phasic{i,j}{k,y}(n) <= 0
                       EMG_phasic_OK{i,j}{k,y}(n)=0;
                       %count = count +1;
                    else
                       EMG_phasic_OK{i,j}{k,y}(n) = EMG_phasic{i,j}{k,y}(n);
                    end
                end
            end
        end
    end
end

%% M matrix organization
% one M matrix for each subject
% one M matrix for each modality
% one line for each muscle
% in each line all the repetitions of a task one after the other and all the
% task one after the other
M_trans = [];
M_imp = [];
count=0;
total_muscles = 10;
total_tasks = 7;
for i=1:size(EMG_phasic_OK,1) %generation of 2 M matrices
    for j =1:size(EMG_phasic_OK,2) %all the task in a given mode must be included in the same line of a matrix
        for m=1:total_muscles %all the muscles
            for r=1:total_tasks %all the repetitions
                if i==1
                    M_trans(m, count+(2700*(r-1)+1): count+(2700*r)) = EMG_phasic_OK{i,j}{r,m};
                end
                if i==2
                    M_imp(m, count+(2700*(r-1)+1): count+(2700*r)) = EMG_phasic_OK{i,j}{r,m};
                end

            end
        end
        count= count + 2700*7;
    end
    count=0;
end

%% Synergic with nnmf
% Assume your sEMG data matrix is called `data` with size 10 x samples
data = M_imp; % Replace this with your actual data
data2 = M_trans;
% Define the number of synergies (factors) you want to extract
num_synergies = 4; % You can change this based on your requirements

% Perform Nonnegative Matrix Factorization
[W, H] = nnmf(data, num_synergies);

%normalization of weight vector
for k = 1:num_synergies    
    nn = norm(W(:,k));
    W(:,k) = W(:,k)./nn;
    H(k,:) = H(k,:).*nn;
end


% Perform Nonnegative Matrix Factorization
[W2, H2] = nnmf(data2, num_synergies);

%normalization of weight vector
for k = 1:num_synergies   
    nn = norm(W2(:,k));
    W2(:,k) = W2(:,k)./nn;
    H2(k,:) = H2(k,:).*nn;
end

% W contains the muscle synergies (basis matrix)
% H contains the activation coefficients (coefficient matrix)

% Display the results
disp('Muscle Synergies (W):');
disp(W);

disp('Activation Coefficients (H):');
disp(H);

% Define a set of colors for plotting
colors = lines(num_synergies);

% Plot the muscle synergies
figure;
for i = 1:num_synergies
    subplot(num_synergies, 1, i);
    plot(W(:, i), 'Color', colors(i, :), 'LineWidth', 1.5);
    title(['Muscle Synergy ', num2str(i)]);
    xlabel('Muscle Index');
    ylabel('Activation Level');
    grid on;
end

% Plot the activation coefficients
figure;
for i = 1:num_synergies
    subplot(num_synergies, 1, i);
    plot(H(i, :), 'Color', colors(i, :), 'LineWidth', 1.5);
    title(['Activation Coefficient ', num2str(i)]);
    xlabel('Time Index');
    ylabel('Activation Level');
    grid on;
end

% Reconstruct the data matrix using W and H
data_reconstructed = W * H;
data_reconstructed2 = W2 * H2;

% Calculate the Variance Accounted For (VAF) for each synergy
vaf_synergies = zeros(1, num_synergies);
for i = 1:num_synergies
    vaf_synergies(i) = (1 - sum((data(:) - data_reconstructed(:)).^2) / sum(data(:).^2)) * 100;
    disp(['VAF for Synergy ', num2str(i), ': ', num2str(vaf_synergies(i)), '%']);
end

% Calculate the Variance Accounted For (VAF)
vaf = (1 - sum((data(:) - data_reconstructed(:)).^2) / sum(data(:).^2)) * 100;

disp(['Variance Accounted For (VAF) - Impedence: ', num2str(vaf), '%']);

vaf_synergies2 = zeros(1, num_synergies);
for i = 1:num_synergies
    vaf_synergies2(i) = (1 - sum((data2(:) - data_reconstructed2(:)).^2) / sum(data2(:).^2)) * 100;
    disp(['VAF for Synergy ', num2str(i), ': ', num2str(vaf_synergies2(i)), '%']);
end

% Calculate the Variance Accounted For (VAF)
vaf = (1 - sum((data2(:) - data_reconstructed2(:)).^2) / sum(data2(:).^2)) * 100;

disp(['Variance Accounted For (VAF) - Transparent: ', num2str(vaf), '%']);
% % Plot original data vs. reconstructed data for comparison
% figure;
% subplot(2, 1, 1);
% imagesc(data);
% colorbar;
% title('Original Data');
% xlabel('Time Index');
% ylabel('Muscle Index');
% 
% subplot(2, 1, 2);
% imagesc(data_reconstructed);
% colorbar;
% title('Reconstructed Data');
% xlabel('Time Index');
% ylabel('Muscle Index');

% % Plot residuals (difference between original and reconstructed data)
% residuals = data - data_reconstructed;
% figure;
% imagesc(residuals);
% colorbar;
% title('Residuals (Original - Reconstructed)');
% xlabel('Time Index');
% ylabel('Muscle Index');

% % Plot histogram of residuals
% figure;
% histogram(residuals(:), 50);
% title('Histogram of Residuals');
% xlabel('Residual Value');
% ylabel('Frequency');

% Plot original data and reconstructed data on the same plot for muscle 5
% muscle_index = 4; % Index of the muscle to plot
% 
% figure;
% plot(data(muscle_index, :), 'b', 'LineWidth', 1.5); % Plot original data in blue
% hold on;
% plot(data_reconstructed(muscle_index, :), 'r.', 'MarkerSize', 5); % Plot reconstructed data in red dots
% title(['Muscle ', num2str(muscle_index)]);
% xlabel('Time Index');
% ylabel('Activation Level');
% grid on;
% legend('Original Data', 'Reconstructed Data');
% hold off;

%% Save the W and H data into a cells
% Uncomment this section if you wanna create tables for W and H
% Load existing tables or create new ones if they don't exist
if isfile('W_table.mat')
    load('W_table.mat', 'W_table');
else
    W_table = {}; 
end

if isfile('H_table.mat')
    load('H_table.mat', 'H_table');
else
    H_table = {};
end

% Create the subject identifier string
subject_id = sprintf('Subject_%d', sub_number);

% Check if the tables are empty and set the initial index
subject_index = [];

if ~isempty(W_table) && ~isempty(W_table{1})
    % Find if the subject already exists
    subject_index = find(strcmp(W_table(:, 1), subject_id));
end

if isempty(subject_index)
    % Append new data if the subject does not exist
    W_table{end+1, 1} = subject_id; % Subject identifier
    W_table{end, 2} = W; % W matrix

    H_table{end+1, 1} = subject_id; % Subject identifier
    H_table{end, 2} = H; % H matrix
else
    % Update existing data if the subject exists
    W_table{subject_index, 2} = W;
    H_table{subject_index, 2} = H;
end

% Save updated tables to .mat files
save('W_table.mat', 'W_table');
save('H_table.mat', 'H_table');


%% 
% Ensure the dimensions are correct
[num_muscles, num_samples] = size(M_imp);
[~, num_synergies] = size(W);

% Normalize each column of W to unit norm
W_normalized = normalize(W, 'norm', 2);

% Initialize H with random values
H = rand(num_synergies, num_samples);

% Parameters for the multiplicative update rule
max_iterations = 1000;
tolerance = 1e-5;

% Perform the NNR algorithm
for iter = 1:max_iterations
    H_old = H;
    
    % Update rule: Hrc = Hrc .* (W' * M_imp) ./ (W' * W * H)
    numerator = W_normalized' * M_imp;
    denominator = (W_normalized' * W_normalized) * H;
    H = H .* (numerator ./ (denominator + eps));
    
    % Check for convergence
    if norm(H - H_old, 'fro') / norm(H_old, 'fro') < tolerance
        disp(['Converged at iteration: ', num2str(iter)]);
        break;
    end
end

% Compute VAF (Variance Accounted For) to evaluate the reconstruction
M_reconstructed = W_normalized * H;
VAF = 1 - sum((M_imp - M_reconstructed).^2, 'all') / sum(M_imp.^2, 'all');

disp(['VAF: ', num2str(VAF)]);



