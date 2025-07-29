%% Loading Data Tables
% This section clears the workspace, closes all figures, and loads necessary data tables from .mat files, renaming them for further use.
clear
close all
clc
load('W_table_IMP.mat');
W_table_IMP = W_table;
clear W_table;
load('W_table_TRANSP.mat');
W_table_TRANSP = W_table;
clear W_table;
load('H_table_IMP.mat');
H_table_IMP = H_table;
clear H_table;
load('H_table_TRANSP.mat');
H_table_TRANSP = H_table;
clear H_table;
%% Reordering the W and H th of TRANSP based on similarity to its IMP mode for each subjects using Cosine metric
% This section reorders the W and H matrices of the TRANSP data to match the IMP mode
% for each subject using the cosine similarity metric. The code processes multiple subjects,
% aligns the matrices, and saves the reordered tables.
% List of subjects
subjects = {"subject3", "subject5", "subject6", "subject7", "subject9", "subject10", "subject11", "subject13"};
num_subjects = length(subjects);

% Initialize the reordered TRANSP tables
W_table_TRANSP_reordered = cell(size(W_table_TRANSP));
H_table_TRANSP_reordered = cell(size(H_table_TRANSP));

% Loop through each subject
for i = 1:length(subjects)
    % Extract the IMP and TRANSP data for the current subject
    W_imp = W_table_IMP{i, 2};
    W_transp = W_table_TRANSP{i, 2};
    H_transp = H_table_TRANSP{i, 2};
    
    % Initialize reordered W_transp and H_transp
    W_transp_reordered = zeros(size(W_transp));
    H_transp_reordered = zeros(size(H_transp));
    used_indices = false(1, 4); % Track which W matrices in TRANSP have been used
    
    % Loop over each W (W1, W2, W3, W4) in IMP and find the best match in TRANSP
    for w_imp_idx = 1:4
        max_cosine_sim = -Inf;
        best_match_idx = -1;
        
        for w_transp_idx = 1:4
            if ~used_indices(w_transp_idx)
                % Reshape the matrices into vectors for cosine similarity calculation
                vec_imp = reshape(W_imp(:, w_imp_idx), [], 1);
                vec_transp = reshape(W_transp(:, w_transp_idx), [], 1);
                
                % Compute the cosine similarity
                cosine_sim = dot(vec_imp, vec_transp) / (norm(vec_imp) * norm(vec_transp));
                
                % Find the maximum cosine similarity
                if cosine_sim > max_cosine_sim
                    max_cosine_sim = cosine_sim;
                    best_match_idx = w_transp_idx;
                end
            end
        end
        
        % Assign the best match to the reordered TRANSP matrix
        W_transp_reordered(:, w_imp_idx) = W_transp(:, best_match_idx);
        H_transp_reordered(w_imp_idx, :) = H_transp(best_match_idx, :);
        used_indices(best_match_idx) = true; % Mark this TRANSP W matrix as used
    end
    
    % Store the reordered TRANSP matrices back to the tables
    W_table_TRANSP_reordered{i, 1} = subjects{i};
    W_table_TRANSP_reordered{i, 2} = W_transp_reordered;
    
    H_table_TRANSP_reordered{i, 1} = subjects{i};
    H_table_TRANSP_reordered{i, 2} = H_transp_reordered;
    
    % Display the result in the command window
    fprintf('Subject: %s\n', subjects{i});
    disp(W_transp_reordered);
    
    % % Plotting
    % figure;
    % 
    % for w = 1:4
    %     % Plot original IMP
    %     subplot(3, 4, w);
    %     bar(W_imp(:, w));
    %     title(['IMP W' num2str(w)]);
    %     xlabel('Muscles');
    %     ylabel('W values');
    % 
    %     % Plot original TRANSP
    %     subplot(3, 4, 4 + w);
    %     bar(W_transp(:, w));
    %     title(['Original TRANSP W' num2str(w)]);
    %     xlabel('Muscles');
    %     ylabel('W values');
    % 
    %     % Plot reordered TRANSP
    %     subplot(3, 4, 8 + w);
    %     bar(W_transp_reordered(:, w));
    %     title(['Reordered TRANSP W' num2str(w)]);
    %     xlabel('Muscles');
    %     ylabel('W values');
    % end
    % 
    % % Adjust layout
    % sgtitle([subjects{i}]);
end

% Save the reordered tables if needed
save('W_table_TRANSP_reordered.mat', 'W_table_TRANSP_reordered');
save('H_table_TRANSP_reordered.mat', 'H_table_TRANSP_reordered');
%% Statistical Analysis: Shapiro-Wilk and Kruskal-Wallis Tests
% This section performs statistical analysis on the reordered TRANSP and IMP matrices
% using the Shapiro-Wilk test for normality and the Kruskal-Wallis test for differences
% between groups.The results are displayed and stored in a table.
% Initialize variables to store results
subjects = W_table_TRANSP_reordered(:, 1);
num_subjects = length(subjects);
num_synergies = 4; % Assuming each data matrix is 10x4

% Initialize result matrices
shapiro_results_TRANSP = zeros(num_subjects, num_synergies, 2); % (subject, synergy, [p-value, W-statistic])
shapiro_results_IMP = zeros(num_subjects, num_synergies, 2);

for i = 1:num_subjects
    % Extract subject name
    subject_name = subjects{i};
    
    % Extract data for the subject
    data_TRANSP = W_table_TRANSP_reordered{i, 2};
    data_IMP = W_table_IMP{i, 2};
    
    for j = 1:num_synergies
        % Extract data for the current synergy
        data_TRANSP_synergy = data_TRANSP(:, j);
        data_IMP_synergy = data_IMP(:, j);
        
        % Perform Shapiro-Wilk test for normality on the current synergy
        [h_TRANSP, p_TRANSP, w_TRANSP] = swtest(data_TRANSP_synergy);
        [h_IMP, p_IMP, w_IMP] = swtest(data_IMP_synergy);
        
        % Store results for the current synergy
        shapiro_results_TRANSP(i, j, :) = [p_TRANSP, w_TRANSP];
        shapiro_results_IMP(i, j, :) = [p_IMP, w_IMP];
        
        % Display results for the current synergy
        fprintf('Subject: %s, Synergy: W%d\n', subject_name, j);
        fprintf('TRANSP Mode - p-value: %f, W-statistic: %f\n', p_TRANSP, w_TRANSP);
        fprintf('IMP Mode - p-value: %f, W-statistic: %f\n\n', p_IMP, w_IMP);
    end
end
% Initialize variables to store results
% Initialize variables to store results
subjects = W_table_TRANSP_reordered(:, 1);
num_subjects = length(subjects);
num_synergies = 4; % Assuming each data matrix is 10x4

% Initialize result matrices
kruskal_results = zeros(num_subjects, num_synergies, 2); % (subject, synergy, [p-value, test-statistic])

% Store subject names and synergy labels for table display
subject_names = [];
synergy_labels = [];

for i = 1:num_subjects
    % Extract subject name
    subject_name = subjects{i};
    
    % Extract data for the subject
    data_TRANSP = W_table_TRANSP_reordered{i, 2};
    data_IMP = W_table_IMP{i, 2};
    
    for j = 1:num_synergies
        % Extract data for the current synergy
        data_TRANSP_synergy = data_TRANSP(:, j);
        data_IMP_synergy = data_IMP(:, j);
        
        % Combine the data for Kruskal-Wallis test
        group = [ones(size(data_TRANSP_synergy)); 2*ones(size(data_IMP_synergy))];
        data = [data_TRANSP_synergy; data_IMP_synergy];
        
        % Perform Kruskal-Wallis test
        [p, tbl, stats] = kruskalwallis(data, group, 'off'); % 'off' suppresses the display of the boxplot
        
        % Store results for the current synergy
        kruskal_results(i, j, :) = [p, tbl{2,5}]; % p-value and test statistic (Chi-square value)
        
        % Store subject name and synergy label for table
        subject_names = [subject_names; {subject_name}];
        synergy_labels = [synergy_labels; {sprintf('W%d', j)}];
    end
end

% Convert results to table format
p_values = kruskal_results(:, :, 1);
chi_square_stats = kruskal_results(:, :, 2);

% Create table
results_table = table(subject_names, synergy_labels, p_values(:), chi_square_stats(:), ...
    'VariableNames', {'Subject', 'Synergy', 'P_Value', 'Chi_Square_Stat'});

% Display table
disp(results_table);

%% Calculating Cosine Similarity
% This section computes the cosine similarity between the reordered TRANSP and IMP matrices
% for each synergy across multiple subjects. The results are organized into a table,
% sorted, and the mean cosine similarity for each subject is calculated and displayed.
% Initialize variables to store results
subjects = W_table_TRANSP_reordered(:, 1);
num_subjects = length(subjects);
num_synergies = 4; % Assuming each data matrix is 10x4

% Initialize result matrices
cosine_similarity_results = zeros(num_subjects, num_synergies); % Initialize for cosine similarity

% Store subject names and synergy labels for table display
subject_names = [];
synergy_labels = [];

for i = 1:num_subjects
    % Extract subject name
    subject_name = subjects{i};
    
    % Extract data for the subject
    data_TRANSP = W_table_TRANSP_reordered{i, 2};
    data_IMP = W_table_IMP{i, 2};
    
    for j = 1:num_synergies
        % Extract data for the current synergy
        data_TRANSP_synergy = data_TRANSP(:, j);
        data_IMP_synergy = data_IMP(:, j);
        
        % Compute cosine similarity
        cosine_similarity = dot(data_TRANSP_synergy, data_IMP_synergy) / (norm(data_TRANSP_synergy) * norm(data_IMP_synergy));
        
        % Store results for the current synergy
        cosine_similarity_results(i, j) = cosine_similarity;
        
        % Store subject name and synergy label for table
        subject_names = [subject_names; {subject_name}];
        synergy_labels = [synergy_labels; {sprintf('W%d', j)}];
    end
end

% Reshape cosine_similarity_results into a single column vector
cosine_similarity_results_reshaped = reshape(cosine_similarity_results, [], 1);

% Create table
results_table = table(subject_names, synergy_labels, cosine_similarity_results_reshaped, ...
    'VariableNames', {'Subject', 'Synergy', 'Cosine_Similarity'});


% Display table
disp(results_table);
% Sort results_table by Cosine_Similarity in ascending order
sorted_results_table = sortrows(results_table, 'Cosine_Similarity');

% Display sorted table
disp(sorted_results_table);

% Assuming sorted_results_table is already calculated and displayed
disp(sorted_results_table);

% Display sorted_results_table
disp(sorted_results_table);

% Convert Subject column to cell array of character vectors
sorted_results_table.Subject = cellstr(sorted_results_table.Subject);

% Calculate mean for each subject
unique_subjects = unique(sorted_results_table.Subject);
mean_cosine_similarity = zeros(length(unique_subjects), 1);

for i = 1:length(unique_subjects)
    subject = unique_subjects{i};
    idx = strcmp(sorted_results_table.Subject, subject);
    cosine_similarities = sorted_results_table.Cosine_Similarity(idx);
    mean_cosine_similarity(i) = mean(cosine_similarities);
end

% Display mean cosine similarity for each subject
fprintf('\nMean Cosine Similarity per Subject:\n');
for i = 1:length(unique_subjects)
    fprintf('Subject %s: %.4f\n', unique_subjects{i}, mean_cosine_similarity(i));
end

%% Calculating P-values and Cosine Similarities Between Two Modalities
% This section calculates cosine similarities between the IMP and reordered TRANSP modalities for each subject.
% It then computes and displays the average cosine similarities for each modality across subjects,
% excluding self-similarities. The results are organized into tables for easy comparison and interpretation.

% Extract the subject names
subjects_IMP = W_table_IMP(:, 1);
subjects_TRANSP = W_table_TRANSP_reordered(:, 1);

% Extract the data columns
data_IMP = W_table_IMP(:, 2);
data_TRANSP = W_table_TRANSP_reordered(:, 2);

% Initialize an array to store cosine similarities
cosine_similarities = zeros(height(W_table_IMP), 1);

% Loop over each subject
for i = 1:height(W_table_IMP)
    % Extract and reshape the data for IMP and TRANSP
    data_IMP_reshaped = reshape(data_IMP{i}, [10, 4]);
    data_TRANSP_reshaped = reshape(data_TRANSP{i}, [10, 4]);
    
    % Flatten the matrices to vectors
    vector_IMP = data_IMP_reshaped(:);
    vector_TRANSP = data_TRANSP_reshaped(:);
    
    % Compute the cosine similarity between the two vectors
    cosine_similarities(i) = dot(vector_IMP, vector_TRANSP) / (norm(vector_IMP) * norm(vector_TRANSP));
end

% Display the cosine similarities along with the subjects
results_table = table(subjects_IMP, cosine_similarities, 'VariableNames', {'Subject', 'Cosine_Similarity'});
disp(results_table);

% Extract subject names and data
subject_names = W_table_IMP(:, 1);

num_subjects = height(W_table_IMP);

% Initialize arrays to store average cosine similarities
avg_cosine_similarity_IMP = zeros(num_subjects, 1);
avg_cosine_similarity_TRANSP = zeros(num_subjects, 1);

% Calculate average cosine similarity for IMP modality
for i = 1:num_subjects
    % Initialize sum of cosine similarities for subject i
    sum_cos_sim_IMP = 0;
    sum_cos_sim_TRANSP = 0;
    
    for j = 1:num_subjects
        if i ~= j
            % Extract and reshape the data for each subject
            data_IMP_i = reshape(W_table_IMP{i, 2}, [10, 4]);
            data_IMP_j = reshape(W_table_IMP{j, 2}, [10, 4]);
            
            data_TRANSP_i = reshape(W_table_TRANSP_reordered{i, 2}, [10, 4]);
            data_TRANSP_j = reshape(W_table_TRANSP_reordered{j, 2}, [10, 4]);

            % Compute cosine similarity between subject i and subject j for both modalities
            cosine_sim_IMP = dot(data_IMP_i(:), data_IMP_j(:)) / (norm(data_IMP_i(:)) * norm(data_IMP_j(:)));
            cosine_sim_TRANSP = dot(data_TRANSP_i(:), data_TRANSP_j(:)) / (norm(data_TRANSP_i(:)) * norm(data_TRANSP_j(:)));
            
            % Accumulate cosine similarities
            sum_cos_sim_IMP = sum_cos_sim_IMP + cosine_sim_IMP;
            sum_cos_sim_TRANSP = sum_cos_sim_TRANSP + cosine_sim_TRANSP;
        end
    end
    
    % Calculate average cosine similarity for subject i (subtract 1 to exclude self-similarity)
    avg_cosine_similarity_IMP(i) = sum_cos_sim_IMP / (num_subjects - 1);
    avg_cosine_similarity_TRANSP(i) = sum_cos_sim_TRANSP / (num_subjects - 1);
end

% Display the average cosine similarities
avg_cosine_table_IMP = table(subject_names, avg_cosine_similarity_IMP, 'VariableNames', {'Subject', 'Avg_Cosine_IMP'});
avg_cosine_table_TRANSP = table(subject_names, avg_cosine_similarity_TRANSP, 'VariableNames', {'Subject', 'Avg_Cosine_TRANSP'});

disp('Average Cosine Similarities for IMP modality:');
disp(avg_cosine_table_IMP);

disp('Average Cosine Similarities for TRANSP modality:');
disp(avg_cosine_table_TRANSP);


%% Plotting W Matrices for IMP and TRANSP Modalities for Specific subject

W_imp = W_table_IMP{8, 2};
W_transp_reordered = W_table_TRANSP_reordered{8, 2};

muscles = {'Brachi', 'Bic-LH', 'Bic-SH', 'Tri', 'Pec', 'AntD', 'MedD', 'PosD', 'Teres', 'Trap'};    %muscle list

% Define the colors for each W
colors = [
    0.1, 0.6, 0.1; % Green for W1
    0.1, 0.4, 0.8; % Blue for W2
    0.8, 0.6, 0.1; % Yellow for W3
    0.7, 0.1, 0.1  % Red for W4
];

figure;

% Loop through each W component for IMP (first row)
for i = 1:4
    subplot(2, 4, i); % Adjust to fit 2 rows and 4 columns of subplots
    barh(W_imp(:, i), 'FaceColor', colors(i, :)); % Use barh for horizontal bar plot
    xlim([0 1.1]); % Adjust to fit the number of channels
    if i == 1
        set(gca, 'YTick', 1:10, 'YTickLabel', muscles,'FontSize',12); % Set muscle names as y-axis labels
    else
        set(gca, 'YTick', 1:10, 'YTickLabel', []); % No labels for other plots
    end
    if i ==1
        ylabel('Channels');
    end
end

% Add title for the first row in the middle
subplot(2, 4, 2); % Use the second subplot in the first row to set the title
title('IMP', 'FontSize', 14, 'FontWeight', 'bold', 'Units', 'normalized', 'Position', [1.1, 1.05, 0]);

% Loop through each W component for TRANSP (second row)
for i = 1:4
    subplot(2, 4, 4 + i); % Adjust to fit 2 rows and 4 columns of subplots
    barh(W_transp_reordered(:, i), 'FaceColor', colors(i, :)); % Use barh for horizontal bar plot
    xlim([0 1.1]); % Adjust to fit the number of channels
    if i == 1
        set(gca, 'YTick', 1:10, 'YTickLabel', muscles,'FontSize',12); % Set muscle names as y-axis labels
    else
        set(gca, 'YTick', 1:10, 'YTickLabel', []); % No labels for other plots
    end
    if i ==1
        ylabel('Channels');
    end
end

% Add title for the second row in the middle
subplot(2, 4, 6); % Use the second subplot in the second row to set the title
title('TRANSP', 'FontSize', 14, 'FontWeight', 'bold', 'Units', 'normalized', 'Position', [1.1, 1.05, 0]);

% Set the font size for all axes
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 9);

% Adjust figure properties
set(gcf, 'PaperPositionMode', 'auto');

% Save as high-quality PNG
% exportgraphics(gcf, 'plotIMPTRANSP_13.png', 'Resolution', 900);
%% Averaging among repetitions among tasks for both H in IMP and TRANSP  

subjects_number_for_plot = 4; % you can change this to show another subject plot
H_final_IMP     = cell(4,9);
H_final_TRANSP  = cell(4,9);
rep = 1:2700:21600;
H_average_IMP = H_table_IMP{subjects_number_for_plot,2};
H_average_TRANSP = H_table_TRANSP_reordered{subjects_number_for_plot,2};

 for i = 1:4
     for j = 0:8
         H_final_sum_IMP    = zeros;
         H_final_sum_TRANSP = zeros;
         for k = 1:7
             H_final_sum_IMP        = H_final_sum_IMP + H_average_IMP(i,(j*18900)+rep(k):(j*18900)+rep(k+1) -1);  
             H_final_sum_TRANSP     = H_final_sum_TRANSP + H_average_TRANSP(i,(j*18900)+rep(k):(j*18900)+rep(k+1) -1);  
         end
         H_final_IMP{i,j+1}     = H_final_sum_IMP/7; 
         H_final_TRANSP{i,j+1}  = H_final_sum_TRANSP/7; 
     end
 end



W_Final_IMP = W_table_IMP{subjects_number_for_plot,2};
W_Final_TRANSP = W_table_TRANSP_reordered{subjects_number_for_plot, 2};


%% Plot specific SUBJECT In two modalities : IMP & TRANSP  W,H
% visualizes the synergies (spatial and temporal) for the IMP and TRANSP modalities of Subject 13,
% highlighting differences through cosine similarity calculations and detailed bar and line plots.
% The results are displayed and saved as high-quality images for further analysis.
% Define titles for the first row of subplots
titles = {'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'T1', 'T2', 'T3'};

% Define figure size and layout
figure('Position', [100, 100, 1300, 900]);

% Define colors for the rows
colors = [
    0.1, 0.6, 0.1; % Dark Green for W1 and row 1 temporal plots
    0.1, 0.4, 0.8; % Dark Blue for W2 and row 2 temporal plots
    0.8, 0.6, 0.1; % Dark Yellow for W3 and row 3 temporal plots
    0.7, 0.1, 0.1  % Dark Red for W4 and row 4 temporal plots
];

% Define gray color for TRANSP mode
gray_color = [0.3, 0.3, 0.3]; % Medium gray

% Plot Spatial Synergies on the Left (W_Final_IMP) and W_Final_TRANSP
for i = 1:4
    % Plot IMP W
    subplot(4, 10, (i-1)*10 + 1); % Adjust to fit 4 rows of subplots
    barh(W_Final_IMP(:, i), 'FaceColor', colors(i, :)); % Set color based on the row
    hold on; % Hold on to overlay TRANSP W
    barh(W_Final_TRANSP(:, i), 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 2); % Black outline for TRANSP with increased width
    title(['W' num2str(i)]);
    ylim([0.5 10.5]); % Adjust to fit the number of channels
    set(gca, 'YDir', 'reverse'); % To match typical horizontal bar orientation
    set(gca, 'YTick', 1:10, 'YTickLabel', muscles); % Set muscle names as y-axis labels
    xlabel('Amplitude');
    ylabel('Channels', 'FontSize', 12);
    hold off; % Release hold after plotting
end

% Plot Temporal Synergies with different line styles for IMP and TRANSP
for i = 1:4
    for j = 1:9
        subplot(4, 10, (i-1)*10 + j + 1); % Shift by 1 column to accommodate spatial plots
        x = linspace(0, 2.7, length(H_final_IMP{i, j})); % Define x-axis values
        xq = linspace(0, 2.7, 10*length(H_final_IMP{i, j})); % Finer x-axis for smoothing
        y_imp = H_final_IMP{i, j};
        y_transp = H_final_TRANSP{i, j};
        
        % Interpolate IMP H with dark color
        yq_imp = interp1(x, y_imp, xq, 'spline');
        % Plot IMP H
        area(x, y_imp, 'FaceColor', colors(i, :), 'EdgeColor', 'none', 'FaceAlpha', 0.5); % Same color for each row, light opacity
        hold on;
        plot(x, y_imp, 'Color', colors(i, :), 'LineWidth', 1.5); % Overlay the line plot with the same color
        
        % Interpolate TRANSP H with black color
        yq_transp = interp1(x, y_transp, xq, 'spline');
        plot(xq, yq_transp, 'Color', gray_color, 'LineWidth', 1.3); % Black color for TRANSP
        
        hold off;
        if i == 1
            title(titles{j});
        end
        % Set y-axis limit to ensure flexibility
        ylim('auto');
        % Remove x-axis and y-axis numbers
        set(gca, 'XTick', [], 'YTick', []);
        % Label x-axis for the last row
        if i == 4
            xlabel('Time (s)');
        end
        % Label y-axis for the first column of temporal plots
        if j == 1
            ylabel(['C' num2str(i)], 'FontSize', 12);
        end
    end
end

% Set the font size for all axes
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 9);

% Adjust figure properties
set(gcf, 'PaperPositionMode', 'auto');

% Save as high-quality PNG
% exportgraphics(gcf, 'Sub7_WH.png', 'Resolution', 600);

%% Re-ordering W and H Indices of IMP and TRANSP for Subjects using Cosine Similarity
% This script reorders the synergies (W and H matrices) of different subjects
% in the IMP and TRANSP modalities by maximizing cosine similarity with a reference subject.
% The reordered synergies are then plotted and saved for further analysis.

% Load all tables
clear
clc

load('W_table_IMP.mat');
W_table_IMP = W_table;
clear W_table;
load('W_table_TRANSP_reordered.mat');
W_table_TRANSP = W_table_TRANSP_reordered;
clear W_table;
load('H_table_IMP.mat');
H_table_IMP = H_table;
clear H_table;
load('H_table_TRANSP_reordered.mat');
H_table_TRANSP = H_table_TRANSP_reordered;
clear H_table;

% List of subjects
subjects = {"subject3", "subject5", "subject6", "subject7", "subject9", "subject10", "subject11", "subject13"};
num_subjects = length(subjects);

% Initialize the reordered tables
W_table_imp_reordered = cell(num_subjects, 2);
W_table_transp_reordered = cell(num_subjects, 2);
H_table_transp_reordered = cell(num_subjects, 2);
H_table_imp_reordered = cell(num_subjects, 2);

% Subject 7 is the reference
ref_subject_idx = 4;
W_imp_ref = W_table_IMP{ref_subject_idx, 2};

% Loop through each subject excluding the reference
for i = 1:num_subjects
    if i == ref_subject_idx
        W_table_imp_reordered{i, 1} = subjects{i};
        W_table_imp_reordered{i, 2} = W_table_IMP{i, 2};
        
        W_table_transp_reordered{i, 1} = subjects{i};
        W_table_transp_reordered{i, 2} = W_table_TRANSP{i, 2};
        
        H_table_transp_reordered{i, 1} = subjects{i};
        H_table_transp_reordered{i, 2} = H_table_TRANSP{i, 2};
        
        H_table_imp_reordered{i, 1} = subjects{i};
        H_table_imp_reordered{i, 2} = H_table_IMP{i, 2};
        continue;
    end

    % Extract the IMP data for the current subject
    W_imp_cmp = W_table_IMP{i, 2};

    % Initialize the reordered indices
    reordered_indices = zeros(1, 4);
    used_indices_ref = false(1, 4); % Track which W matrices of reference have been used
    used_indices_cmp = false(1, 4); % Track which W matrices of the current subject have been used

    % Loop over each W (W1, W2, W3, W4) in the reference IMP
    for w_ref_idx = 1:4
        max_cosine_sim = -Inf;
        best_match_idx = -1;

        % Compare each W in IMP of the reference subject with each W in IMP of the current subject
        for w_cmp_idx = 1:4
            if ~used_indices_cmp(w_cmp_idx)
                % Reshape the matrices into vectors for cosine similarity calculation
                vec_imp_ref = reshape(W_imp_ref(:, w_ref_idx), [], 1);
                vec_imp_cmp = reshape(W_imp_cmp(:, w_cmp_idx), [], 1);

                % Compute the cosine similarity
                cosine_sim = dot(vec_imp_ref, vec_imp_cmp) / (norm(vec_imp_ref) * norm(vec_imp_cmp));

                % Find the maximum cosine similarity
                if cosine_sim > max_cosine_sim
                    max_cosine_sim = cosine_sim;
                    best_match_idx = w_cmp_idx;
                end
            end
        end

        % Assign the best match index
        reordered_indices(w_ref_idx) = best_match_idx;
        used_indices_cmp(best_match_idx) = true; % Mark this W matrix of the current subject as used
    end

    % Reorder the IMP matrix based on the reordered indices
    W_reordered = W_imp_cmp(:, reordered_indices);

    % Reorder the TRANSP and H matrices based on the same reordered indices
    W_transp_reordered = W_table_TRANSP{i, 2}(:, reordered_indices);
    H_transp_reordered = H_table_TRANSP{i, 2}(reordered_indices, :);
    H_imp_reordered = H_table_IMP{i, 2}(reordered_indices, :);

    % Store the reordered data in the tables
    W_table_imp_reordered{i, 1} = subjects{i};
    W_table_imp_reordered{i, 2} = W_reordered;

    W_table_transp_reordered{i, 1} = subjects{i};
    W_table_transp_reordered{i, 2} = W_transp_reordered;

    H_table_transp_reordered{i, 1} = subjects{i};
    H_table_transp_reordered{i, 2} = H_transp_reordered;

    H_table_imp_reordered{i, 1} = subjects{i};
    H_table_imp_reordered{i, 2} = H_imp_reordered;

    % Display the result in the command window
    fprintf('Subject: %s\n', subjects{i});
    disp(reordered_indices);
end

% Plot original IMPs for all subjects
figure;
sgtitle('Original IMPs for All Subjects');
for i = 1:num_subjects
    W_imp = W_table_IMP{i, 2};
    for w = 1:4
        subplot(num_subjects, 4, (i-1)*4 + w);
        bar(W_imp(:, w));
        title(['W' num2str(w)]);
        ylabel(subjects{i});
    end
end

% Plot reordered IMPs for all subjects
figure;
sgtitle('Reordered IMPs for All Subjects');
for i = 1:num_subjects
    W_imp = W_table_imp_reordered{i, 2};
    for w = 1:4
        subplot(num_subjects, 4, (i-1)*4 + w);
        bar(W_imp(:, w));
        title(['W' num2str(w)]);
        ylabel(subjects{i});
    end
end
% 
% % Plot original TRANSPs for all subjects
% figure;
% sgtitle('Original TRANSPs for All Subjects');
% for i = 1:num_subjects
%     W_transp = W_table_TRANSP{i, 2};
%     for w = 1:4
%         subplot(num_subjects, 4, (i-1)*4 + w);
%         bar(W_transp(:, w));
%         title(['W' num2str(w)]);
%         ylabel(subjects{i});
%     end
% end
% 
% % Plot reordered TRANSPs for all subjects
% figure;
% sgtitle('Reordered TRANSPs for All Subjects');
% for i = 1:num_subjects
%     W_transp = W_table_transp_reordered{i, 2};
%     for w = 1:4
%         subplot(num_subjects, 4, (i-1)*4 + w);
%         bar(W_transp(:, w));
%         title(['W' num2str(w)]);
%         ylabel(subjects{i});
%     end
% end

% Display the reordered indices of the reference subject
W_table_imp_reordered{ref_subject_idx, 1} = subjects{ref_subject_idx};
W_table_imp_reordered{ref_subject_idx, 2} = W_table_IMP{ref_subject_idx, 2}; % No reordering needed for the reference subject

% Save the reordered tables
% save('W_table_transp_reordered_2.mat', 'W_table_transp_reordered');
% save('H_table_transp_reordered_2.mat', 'H_table_transp_reordered');
% save('H_table_imp_reordered_2.mat', 'H_table_imp_reordered');
% save('W_table_imp_reordered_indices_2.mat', 'W_table_imp_reordered');
%% Cosine similarity between sub3 and sub10 among  others
% This MATLAB script calculates and displays the cosine similarities of the W matrices between
% selected subjects (subject3 and subject10) and other subjects (excluding subject13)
% in both IMP and TRANSP modalities.
% The script also calculates and prints the average cosine similarities for each W matrix.

% Calculate and display cosine similarities for specific subjects, excluding subject 13
subjects_of_interest = [1, 6]; % Indices of subject3 and subject10
num_W = 4;
excluded_subject_idx = 8; % Index of subject13

% Initialize matrices to store the cosine similarities
cosine_similarities_imp = zeros(num_W, num_subjects - 2, length(subjects_of_interest)); % -2 to exclude reference subject and subject13
cosine_similarities_transp = zeros(num_W, num_subjects - 2, length(subjects_of_interest));

for subj_interest_idx = 1:length(subjects_of_interest)
    subj_idx = subjects_of_interest(subj_interest_idx);
    for w_idx = 1:num_W
        fprintf('Comparisons for %s, W%d\n', subjects{subj_idx}, w_idx);
        count = 1;
        for other_subj_idx = 1:num_subjects
            if other_subj_idx ~= subj_idx && other_subj_idx ~= excluded_subject_idx
                % For IMP
                vec1_imp = reshape(W_table_imp_reordered{subj_idx, 2}(:, w_idx), [], 1);
                vec2_imp = reshape(W_table_imp_reordered{other_subj_idx, 2}(:, w_idx), [], 1);
                cosine_sim_imp = dot(vec1_imp, vec2_imp) / (norm(vec1_imp) * norm(vec2_imp));
                cosine_similarities_imp(w_idx, count, subj_interest_idx) = cosine_sim_imp;
                fprintf('  %s W%d (IMP): %f\n', subjects{other_subj_idx}, w_idx, cosine_sim_imp);

                % For TRANSP
                vec1_transp = reshape(W_table_transp_reordered{subj_idx, 2}(:, w_idx), [], 1);
                vec2_transp = reshape(W_table_transp_reordered{other_subj_idx, 2}(:, w_idx), [], 1);
                cosine_sim_transp = dot(vec1_transp, vec2_transp) / (norm(vec1_transp) * norm(vec2_transp));
                cosine_similarities_transp(w_idx, count, subj_interest_idx) = cosine_sim_transp;
                fprintf('  %s W%d (TRANSP): %f\n', subjects{other_subj_idx}, w_idx, cosine_sim_transp);

                count = count + 1;
            end
        end
    end
end

fprintf('\nAverage Cosine Similarities (IMP):\n');
for subj_interest_idx = 1:length(subjects_of_interest)
    subj_idx = subjects_of_interest(subj_interest_idx);
    fprintf('Averages for %s:\n', subjects{subj_idx});
    for w_idx = 1:num_W
        avg_cosine_sim_imp = mean(cosine_similarities_imp(w_idx, :, subj_interest_idx));
        fprintf('  W%d: %f\n', w_idx, avg_cosine_sim_imp);
    end
end

fprintf('\nAverage Cosine Similarities (TRANSP):\n');
for subj_interest_idx = 1:length(subjects_of_interest)
    subj_idx = subjects_of_interest(subj_interest_idx);
    fprintf('Averages for %s:\n', subjects{subj_idx});
    for w_idx = 1:num_W
        avg_cosine_sim_transp = mean(cosine_similarities_transp(w_idx, :, subj_interest_idx));
        fprintf('  W%d: %f\n', w_idx, avg_cosine_sim_transp);
    end
end

%% Averaging W ,H among subjects

% Initialize a matrix to store the sum of all data matrices
W_sum_IMP       = zeros(10, 4);
W_sum_TRANSP    = zeros(10, 4);
H_sum_IMP       = zeros(4, 170100);
H_sum_TRANSP    = zeros(4, 170100);

% Loop through each cell in the second column and accumulate the sum
subject_considered = [2,3,4,5,7,8] ;
for i = subject_considered
    W_sum_IMP = W_sum_IMP + W_table_imp_reordered{i, 2};
    W_sum_TRANSP = W_sum_TRANSP + W_table_transp_reordered{i, 2};
    H_sum_IMP = H_sum_IMP + H_table_imp_reordered{i, 2};
    H_sum_TRANSP = H_sum_TRANSP + H_table_transp_reordered{i, 2};
end

% Calculate the mean by dividing by the number of data matrices (8)
W_mean_IMP       = W_sum_IMP/numel(subject_considered);
W_mean_TRANSP    = W_sum_TRANSP/numel(subject_considered);
H_mean_IMP       = H_sum_IMP/numel(subject_considered);
H_mean_TRANSP    = H_sum_TRANSP/numel(subject_considered);
%% Calculate P-value for Averaged W, Kruskal-Wallis Test, and Cosine Similarity
% This script performs the following statistical analyses between the W_mean_IMP and W_mean_TRANSP datasets:
% 
% 1-----Paired t-test to compare the means of the two datasets.
% 2-----Kruskal-Wallis Test for each variable (W1, W2, W3, W4) to assess if the samples originate from the same distribution.
% 3-----Cosine Similarity between corresponding columns of the two datasets.


% Reshape the tables to 1D vectors
data1 = W_mean_IMP(:);
data2 = W_mean_TRANSP(:);

% Perform the paired t-test
[~, P_value] = ttest(data1, data2);

% Display the P-value
disp(['P-value: ', num2str(P_value)]);

% Generate example data (replace with your actual data)
data1 = W_mean_IMP;  % Example dataset 1 (10 rows, 4 columns)
data2 = W_mean_TRANSP;  % Example dataset 2 (10 rows, 4 columns)

% Initialize arrays to store p-values and H-statistics
p_values = zeros(1, 4);  % Array to store p-values for each W
H_stats = zeros(1, 4);   % Array to store H-statistics for each W

% Loop over each variable (W1, W2, W3, W4)
for i = 1:4
    % Extract current variable from both datasets
    variable_data1 = data1(:, i);
    variable_data2 = data2(:, i);
    
    % Combine the data and create groups vector
    all_data = [variable_data1; variable_data2];
    groups = [ones(10, 1); 2*ones(10, 1)];  % Group 1 for data1, Group 2 for data2
    
    % Perform Kruskal-Wallis test
    [p, tbl, stats] = kruskalwallis(all_data, groups, 'off');
    
    % Store p-value and H-statistic
    p_values(i) = p;
    H_stats(i) = tbl{2, 5};  % Extract the H-statistic from the table
end

% Display results
disp('Kruskal-Wallis Test Results:');
disp(['Variable   p-value   H-statistic']);
for i = 1:4
    fprintf('W%d         %.4f     %.2f\n', i, p_values(i), H_stats(i));
end
% Example data (replace with your actual data)
data1 = W_mean_IMP;  % Example dataset 1 (10 rows, 4 columns)
data2 = W_mean_TRANSP;  % Example dataset 2 (10 rows, 4 columns)

% Initialize arrays to store cosine similarities
cosine_similarities = zeros(1, size(data1, 2));  % Array to store cosine similarities

% Loop over each variable (W1, W2, W3, W4)
for i = 1:size(data1, 2)
    % Extract current variable (column) from both datasets
    variable_data1 = data1(:, i);
    variable_data2 = data2(:, i);
    
    % Compute cosine similarity
    cosine_sim = dot(variable_data1, variable_data2) / (norm(variable_data1) * norm(variable_data2));
    
    % Store cosine similarity
    cosine_similarities(i) = cosine_sim;
end

% Display results
disp('Cosine Similarities:');
disp(['Variable   Cosine Similarity']);
for i = 1:size(data1, 2)
    fprintf('W%d         %.4f\n', i, cosine_similarities(i));
end



%% Averaging H Among Repetitions for Both IMP and TRANSP Tasks 
% This script calculates the average of H matrices
% for both IMP and TRANSP tasks across different repetitions and tasks.
H_final_IMP     = cell(4,9);
H_final_TRANSP  = cell(4,9);
rep = 1:2700:21600;
H_average_IMP = H_mean_IMP;
H_average_TRANSP = H_mean_TRANSP;

 for i = 1:4
     for j = 0:8
         H_final_sum_IMP    = zeros;
         H_final_sum_TRANSP = zeros;
         for k = 1:7
             H_final_sum_IMP        = H_final_sum_IMP + H_average_IMP(i,(j*18900)+rep(k):(j*18900)+rep(k+1) -1);  
             H_final_sum_TRANSP     = H_final_sum_TRANSP + H_average_TRANSP(i,(j*18900)+rep(k):(j*18900)+rep(k+1) -1);  
         end
         H_final_IMP{i,j+1}     = H_final_sum_IMP/7; 
         H_final_TRANSP{i,j+1}  = H_final_sum_TRANSP/7; 
     end
 end

W_final_IMP = W_mean_IMP;
W_final_TRANSP = W_mean_TRANSP;

%% Plot Average among subjects In two modalities : IMP & TRANSP  W,H
% Define the angles in degrees for the specified tasks
angles_deg = [220, 180, 140, 320, 0, 40];

% Convert the angles to radians for plotting
angles_rad = deg2rad(angles_deg);

% Define the common length for all tasks
length_common = 3.5;

% Define the task parameters as a cell array
tasks = {
    'L1', angles_rad(1), length_common;
    'M1', angles_rad(2), length_common;
    'T1', angles_rad(3), length_common;
    'L3', angles_rad(4), length_common;
    'M3', angles_rad(5), length_common;
    'T3', angles_rad(6), length_common;
};

% Define the colors for the tasks
colors = [
    0.1, 0.6, 0.1; % Dark Green for first plot
    0.1, 0.4, 0.8; % Dark Blue for second plot
    0.8, 0.6, 0.1; % Dark Yellow for third plot
    0.7, 0.1, 0.1  % Dark Red for fourth plot
];

% Define plot settings for each plot
plot_settings = {
    {'L3', 'M3', 'T3'}, length_common - 1.5, colors(1, :); % First plot
    {'L1', 'M1', 'T1'}, length_common, colors(2, :);       % Second plot
    {'L3', 'M3', 'T3'}, length_common, colors(3, :);       % Third plot
    {'L1', 'M1', 'T1'}, length_common - 1.5, colors(4, :); % Fourth plot
};

% Define titles for the first row of subplots
muscles = {'Brachi', 'Bic-LH', 'Bic-SH', 'Tri', 'Pec', 'AntD', 'MedD', 'PosD', 'Teres', 'Trap'};    % Muscle list
titles = {'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'T1', 'T2', 'T3'};

% Define figure size and layout
figure('Position', [100, 100, 1600, 900]); % Adjusted figure size

% Define colors for the rows
colors_original = [
    0.1, 0.6, 0.1; % Dark Green for W1 and row 1 temporal plots
    0.1, 0.4, 0.8; % Dark Blue for W2 and row 2 temporal plots
    0.8, 0.6, 0.1; % Dark Yellow for W3 and row 3 temporal plots
    0.7, 0.1, 0.1  % Dark Red for W4 and row 4 temporal plots
];

% Define gray color for TRANSP mode
gray_color = [0.3, 0.3, 0.3]; % Medium gray

% Plot Spatial Synergies on the Left (W_Final_IMP) and W_Final_TRANSP
for i = 1:4
    % Plot IMP W
    subplot(4, 11, (i-1)*11 + 1); % Adjust to fit 4 rows of subplots
    barh(W_final_IMP(:, i), 'FaceColor', colors_original(i, :)); % Set color based on the row
    hold on; % Hold on to overlay TRANSP W
    barh(W_final_TRANSP(:, i), 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 2); % Black outline for TRANSP with increased width
    title(['W' num2str(i)]);
    ylim([0.5 10.5]); % Adjust to fit the number of channels
    set(gca, 'YDir', 'reverse'); % To match typical horizontal bar orientation
    set(gca, 'YTick', 1:10, 'YTickLabel', muscles); % Set muscle names as y-axis labels
    xlabel('Amplitude');
    ylabel('Channels', 'FontSize', 12);
    hold off; % Release hold after plotting
end

% Plot Temporal Synergies with different line styles for IMP and TRANSP
for i = 1:4
    for j = 1:9
        subplot(4, 11, (i-1)*11 + j + 1); % Shift by 1 column to accommodate spatial plots
        x = linspace(0, 2.7, length(H_final_IMP{i, j})); % Define x-axis values
        xq = linspace(0, 2.7, 10*length(H_final_IMP{i, j})); % Finer x-axis for smoothing
        y_imp = H_final_IMP{i, j};
        y_transp = H_final_TRANSP{i, j};

        % Interpolate IMP H with dark color
        yq_imp = interp1(x, y_imp, xq, 'spline');
        % Plot IMP H
        area(x, y_imp, 'FaceColor', colors_original(i, :), 'EdgeColor', 'none', 'FaceAlpha', 0.5); % Same color for each row, light opacity
        hold on;
        plot(x, y_imp, 'Color', colors_original(i, :), 'LineWidth', 1.5); % Overlay the line plot with the same color

        % Interpolate TRANSP H with black color
        yq_transp = interp1(x, y_transp, xq, 'spline');
        plot(xq, yq_transp, 'Color', gray_color, 'LineWidth', 1.3); % Black color for TRANSP

        hold off;
        if i == 1
            title(titles{j});
        end
        % Set y-axis limit to ensure flexibility
        ylim([0, 0.5]);
        % Remove x-axis and y-axis numbers
        % set(gca, 'XTick', [], 'YTick', []);
        % Label x-axis for the last row
        if i == 4
            xlabel('Time (s)');
        end
        % Label y-axis for the first column of temporal plots
        if j == 1
            ylabel(['C' num2str(i)], 'FontSize', 12);
        end
    end
end

% Add new plots to the last column
for plot_idx = 1:4
    % Extract the settings for the current plot
    task_list = plot_settings{plot_idx, 1};
    current_length = plot_settings{plot_idx, 2};
    current_color = plot_settings{plot_idx, 3};
    
    % Create a subplot
   
    subplot(4, 11, [plot_idx * 11, plot_idx * 11]); % Enlarge circle plot in the last column
   
    
    hold on;

    % Plot the circle using a slight gray color with dashed-dot pattern
    theta_circle = linspace(0, 2*pi, 100);
    x_circle = length_common * cos(theta_circle);
    y_circle = length_common * sin(theta_circle);
    plot(x_circle, y_circle, 'Color', [0.8, 0.8, 0.8], 'LineStyle', '-'); % light gray dashed circle

    % Overlay with gray diameter lines based on angles_deg
    for k = 1:numel(angles_deg)
        angle_rad = deg2rad(angles_deg(k));
        x_line = length_common * cos(angle_rad);
        y_line = length_common * sin(angle_rad);
        line([-x_line, x_line], [-y_line, y_line], 'Color', [0.8, 0.8, 0.8], 'LineStyle', '-');
    end

    % Set up the axes
    axis equal;
    axis off; % Remove axes
    set(gca, 'Color', 'none'); % Make axes background transparent

    % Plot the arrows for specified tasks and label them
    for i = 1:numel(task_list)
        task_name = task_list{i};

        % Find the task parameters in the cell array
        for j = 1:size(tasks, 1)
            if strcmp(tasks{j, 1}, task_name)
                angle = tasks{j, 2};

                % Calculate the end points of the arrows
                x_end = current_length * cos(angle);
                y_end = current_length * sin(angle);

                % Plot the arrow
                quiver(0, 0, x_end, y_end, 0, 'MaxHeadSize', 0.8, 'LineWidth', 2, 'Color', current_color);
                
                if (strcmp(task_name, 'L3') || strcmp(task_name, 'M3')|| strcmp(task_name, 'T3') )
                    % Calculate position for text label on the perimeter
                    text_x = length_common * cos(angle) + 1;
                    text_y = length_common * sin(angle) + 0.1;
                else
                    text_x = length_common * cos(angle) - 1;
                    text_y = length_common * sin(angle) + 0.1;
                end
                % Add text label on the perimeter of the circle
                text(text_x, text_y, task_name, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 14);

                break;
            end
        end
    end

    % Customize the plot
    
    hold off;
end


% Set the font size for all axes
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 9);

% Adjust figure properties
set(gcf, 'PaperPositionMode', 'auto');

% Save as high-quality PNG
% exportgraphics(gcf, 'merged_plots_enlarged.png', 'Resolution', 600);
