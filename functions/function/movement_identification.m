function [onset, offset] = movement_identification(joints_displ, fs, angular_velocity_threshold, neighborhood)
    time_comm = joints_displ(:,1);
    
    %WE USE JOINT 4 BECAUSE PATTERN IS MORE REGULAR
    % Filtering kinematics data (frequency has been determined by analyzing
    % the periodogram
    [PSD_J4_displ frequencies_J4] = periodogram(joints_displ(:,4),[],[],fs);
    figure;
    semilogy(frequencies_J4, PSD_J4_displ);
    title('Power Spectral Density (PSD) of Kinematics Data');
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency Density');
    grid on;
    [B1, A1] = butter(5, 20/(fs/2), 'low');
    joints_displ_J4_filt = filtfilt(B1, A1, joints_displ(:,4));
    angular_velocity_J4_filt = gradient(joints_displ_J4_filt, 1/fs);

    % Finding peaks related to maximum velocity
    %We find the peaks that are related to the time instants where we have the
    %maximum velocity (we add the control "MinPeakDistance" as the ideal period
    %of the movement in order to find the first peaks related to the
    %antigravity movement 
    [peak, loc_peak] = findpeaks(-angular_velocity_J4_filt, "MinPeakDistance", 5*fs, "MinPeakHeight", 2);
    loc_peak_ok = loc_peak(peak > 0);
    
    % Plotting identified peaks
    figure;
    plot(time_comm, angular_velocity_J4_filt, 'r');
    hold on;
    plot(time_comm(loc_peak_ok), angular_velocity_J4_filt(loc_peak_ok), 'bo');
    legend('Angular Velocity', 'Peaks');
    title('Angular velocity with peaks related to reaching task');
    xlabel('Time (s)','FontSize',12);
    ylabel('Angular Velocity (Degree/Second)','FontSize',12);

    set(gca);
    % Adjust figure properties
    set(gcf, 'PaperPositionMode', 'auto');

    % Save as high-quality PNG
    exportgraphics(gcf, 'high_quality_plot2.png', 'Resolution', 300);

    % Finding time instants with low angular velocity
    idx_still = find(angular_velocity_J4_filt > -angular_velocity_threshold & angular_velocity_J4_filt < angular_velocity_threshold);

    % Initializing onset and offset arrays
    onset = zeros(length(loc_peak_ok), 1);
    offset = zeros(length(loc_peak_ok), 1);

    % Finding onset points
    for i = 1:length(loc_peak_ok)
        candidate_points_onset = idx_still(idx_still < loc_peak_ok(i));
        for j = 1:length(candidate_points_onset)
            time_window_robustness = (candidate_points_onset(end-j) - neighborhood):candidate_points_onset(end-j);
            count = sum(ismember(candidate_points_onset, time_window_robustness));
            if (count > 0.8*neighborhood)
                onset(i) = candidate_points_onset(end-j);
                break;
            end
        end
    end

    % Finding offset points
    for i = 1:length(loc_peak_ok)
        candidate_points_offset = idx_still(idx_still > loc_peak_ok(i));
        for j = 1:length(candidate_points_offset)
            time_window_robustness = (candidate_points_offset(j) - neighborhood):candidate_points_offset(j);
            count = sum(ismember(candidate_points_offset, time_window_robustness));
            if (count > 0.8*neighborhood)
                offset(i) = candidate_points_offset(j);
                break;
            end
        end
    end

    % Plotting onset and offset points
    % figure;
    % plot(time_comm, angular_velocity_J4_filt, 'r');
    % hold on;
    % plot(time_comm(offset), angular_velocity_J4_filt(offset), 'bo');
    % plot(time_comm(onset), angular_velocity_J4_filt(onset), 'ko');
    % legend('Angular Velocity', 'Offset Points', 'Onset Points');
    % title('Angular Velocity with Onset and Offset Points');
    % xlabel('Time (ms)');
    % ylabel('Angular Velocity');
end