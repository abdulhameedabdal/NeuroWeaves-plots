%% Step 1: Initialize Variables
clear; clc; close all;

% Add Intan RHS Functions Path 
 addpath('');<-- add your path

% Parameters
file_duration = 60;                % Duration of each file in seconds
num_channels = 6;                  % Total number of channels (1-6)
stim_threshold = 5;                % Threshold for detecting stimulation
prestim_window = 0.01;             % Pre-stimulus window (10 ms)
poststim_window = 0.1;             % Post-stimulus window (100 ms)
low_cutoff = 1;                    % Bandpass filter lower cutoff (Hz)
high_cutoff = 100;                 % Bandpass filter upper cutoff (Hz)
notch_freq = 60;                   % Notch filter frequency (Hz)
time_start_baseline = 1;           % Start time for baseline (seconds)
time_end_baseline = 5;             % End time for baseline (seconds)

% Initialize Continuous Variables
appended_data = [];
appended_time = [];
appended_stim_waveform = [];
current_file_index = 0;
fs = []; % Sampling rate

% Loop to Append Files 
while true
    % Increment file index
    current_file_index = current_file_index + 1;
    fprintf('Load file %d...\n', current_file_index);

    % Load Intan RHS file
    read_Intan_RHS2000_file;

    % Verify Required Variables
    if ~exist('board_adc_data', 'var') || ~exist('amplifier_data', 'var') || ~exist('t', 'var')
        fprintf('Error: Required data missing. Exiting...\n');
        break;
    end

    % Sampling Rate from First File
    if isempty(fs)
        fs = frequency_parameters.amplifier_sample_rate;
        time_step = 1 / fs; % Time increment per sample
    end

    % Retrieve stimulation waveform and EEG data
    stim_waveform = board_adc_data(1, :); % Assume ADC Channel 1 for stimulation
    eeg_data = amplifier_data;            % EEG Data (6 Channels)

    % Adjust Time for Continuity
    if isempty(appended_time)
        continuous_time = t; % First file: use original time vector
    else
        last_time = appended_time(end);
        continuous_time = last_time + time_step + (0:time_step:(length(stim_waveform)-1)*time_step);
    end

    % Append Data
    appended_stim_waveform = [appended_stim_waveform, stim_waveform];
    appended_data = [appended_data, eeg_data];
    appended_time = [appended_time, continuous_time];

    % Save Progress
    save('appended_rhs_data.mat', 'appended_data', 'appended_time', 'appended_stim_waveform', 'fs');
    fprintf('File %d appended successfully.\n', current_file_index);

    % Ask user to load another file
    choice = questdlg('Load another file?', 'Continue?', 'Yes', 'No', 'Yes');
    if strcmp(choice, 'No')
        break;
    end
end

% Detect Stimulation Onsets
rising_edges = find(diff(appended_stim_waveform > stim_threshold) == 1);
fprintf('Total number of trials detected: %d\n', length(rising_edges));

% Window sizes
samples_prestim = round(prestim_window * fs);
samples_poststim = round(poststim_window * fs);
window_length = samples_prestim + samples_poststim;
time_axis = (-samples_prestim:samples_poststim-1) / fs * 1000; % Time in ms

% Design Filters
bpFilt = designfilt('bandpassiir', 'FilterOrder', 4, 'HalfPowerFrequency1', low_cutoff, 'HalfPowerFrequency2', high_cutoff, 'SampleRate', fs);
notchFilt = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1', notch_freq - 1, 'HalfPowerFrequency2', notch_freq + 1, 'SampleRate', fs);

%% Plot 1: Evoked responses of all channels
figure('Name', 'Evoked Responses - All Channels', 'NumberTitle', 'off');
for ch = 1:num_channels
    response_data = appended_data(ch, :);
    response_data_filtered = filtfilt(notchFilt, response_data);
    response_data_filtered = filtfilt(bpFilt, response_data_filtered);

    % Segment data into trials
    response_trials = zeros(length(rising_edges), window_length);
    for trial = 1:length(rising_edges)
        stim_idx = rising_edges(trial);
        start_idx = max(stim_idx - samples_prestim, 1);
        end_idx = min(stim_idx + samples_poststim - 1, length(response_data_filtered));
        if (end_idx - start_idx + 1) == window_length
            response_trials(trial, :) = response_data_filtered(start_idx:end_idx);
        else
            response_trials(trial, :) = NaN;
        end
    end
    response_trials = response_trials(~any(isnan(response_trials), 2), :);
    average_response = mean(response_trials, 1);

    % Plot
    subplot(3, 2, ch);
    hold on;
    for trial = 1:size(response_trials, 1)
        plot(time_axis, response_trials(trial, :), 'Color', [0.8, 0.8, 0.8]); % Grey individual trials
    end
    plot(time_axis, average_response, 'k', 'LineWidth', 2); % Black trial average
    plot(0, 0, 'ro', 'MarkerSize', 8, 'LineWidth', 2); % Red stimulus marker
    hold off;
    xlabel('Time (ms)'); ylabel('Amplitude (\muV)');
    title(['Channel ', num2str(ch)]);
    xlim([-10, 100]); ylim([-100, 100]); grid on;
end

%% Plot 2: Baseline activity for Channels 5 and 2
baseline_channels = [5, 2];
time_indices = appended_time >= time_start_baseline & appended_time <= time_end_baseline;
time_selected = appended_time(time_indices);
data_selected = appended_data(baseline_channels, time_indices);
raw_channel5 = data_selected(1, :);
raw_channel2 = data_selected(2, :);

% Pad data to mitigate filtering edge effects
padding_length = 1000; % Padding length in samples
padded_raw_channel5 = [flip(raw_channel5(1:padding_length)), raw_channel5, flip(raw_channel5(end-padding_length+1:end))];
padded_raw_channel2 = [flip(raw_channel2(1:padding_length)), raw_channel2, flip(raw_channel2(end-padding_length+1:end))];

% Filter padded data
filtered_channel5 = filtfilt(bpFilt, filtfilt(notchFilt, padded_raw_channel5));
filtered_channel2 = filtfilt(bpFilt, filtfilt(notchFilt, padded_raw_channel2));

% Remove padding
filtered_channel5 = filtered_channel5(padding_length+1:end-padding_length);
filtered_channel2 = filtered_channel2(padding_length+1:end-padding_length);

figure('Name', 'Baseline Activity - Channels 5 and 2', 'NumberTitle', 'off');
offsets = [0, -150, -300, -450];
colors = {'#CB0505', '#00AAFF'}; % Colors for minimal design later

% Plot raw Channel 5
hold on;
plot(time_selected, raw_channel5 + offsets(1), 'k', 'LineWidth', 1.5);
% Plot filtered Channel 5
plot(time_selected, filtered_channel5 + offsets(2), 'k', 'LineWidth', 1.5);
% Plot raw Channel 2
plot(time_selected, raw_channel2 + offsets(3), 'k', 'LineWidth', 1.5);
% Plot filtered Channel 2
plot(time_selected, filtered_channel2 + offsets(4), 'k', 'LineWidth', 1.5);
hold off;

title('Baseline Activity (1-5s)');
xlabel('Time (s)'); ylabel('Amplitude (\muV)');
grid on; legend('Raw Ch5', 'Filtered Ch5', 'Raw Ch2', 'Filtered Ch2');

%% Plot 3: Minimal design for baseline activity
figure('Name', 'Minimal Baseline Activity - Channels 5 and 2', 'NumberTitle', 'off');
hold on;
plot(time_selected, raw_channel5 + offsets(1), 'Color', colors{1}, 'LineWidth', 1.5);
plot(time_selected, raw_channel2 + offsets(3), 'Color', colors{2}, 'LineWidth', 1.5);
hold off;
title('Minimal Baseline Activity');
axis off;

%% Plot 4: Subplots for Evoked Responses of Channels 2,3,5,6
figure('Name', 'Evoked Responses Subplots', 'NumberTitle', 'off');
channels_to_plot = [2, 3, 5, 6];
for idx = 1:length(channels_to_plot)
    subplot(2, 2, idx);
    ch = channels_to_plot(idx);
    response_data = appended_data(ch, :);
    response_data_filtered = filtfilt(notchFilt, response_data);
    response_data_filtered = filtfilt(bpFilt, response_data_filtered);

    % Segment data into trials
    response_trials = zeros(length(rising_edges), window_length);
    for trial = 1:length(rising_edges)
        stim_idx = rising_edges(trial);
        start_idx = max(stim_idx - samples_prestim, 1);
        end_idx = min(stim_idx + samples_poststim - 1, length(response_data_filtered));
        if (end_idx - start_idx + 1) == window_length
            response_trials(trial, :) = response_data_filtered(start_idx:end_idx);
        else
            response_trials(trial, :) = NaN;
        end
    end
    response_trials = response_trials(~any(isnan(response_trials), 2), :);
    average_response = mean(response_trials, 1);

    % Plot
    hold on;
    for trial = 1:size(response_trials, 1)
        plot(time_axis, response_trials(trial, :), 'Color', [0.8, 0.8, 0.8]); % Grey individual trials
    end
    plot(time_axis, average_response, 'k', 'LineWidth', 2); % Black trial average
    plot(0, 0, 'ro', 'MarkerSize', 8, 'LineWidth', 2); % Red stimulus marker
    hold off;
    xlabel('Time (ms)'); ylabel('Amplitude (\muV)');
    title(['Channel ', num2str(ch)]);
    xlim([-10, 100]); ylim([-100, 100]);
end

%% Plot 5: Minimal Design Subplots for Channels 2,3,5,6
figure('Name', 'Minimal Design - Evoked Responses', 'NumberTitle', 'off');
for idx = 1:length(channels_to_plot)
    subplot(2, 2, idx);
    ch = channels_to_plot(idx);
    response_data = appended_data(ch, :);
    response_data_filtered = filtfilt(notchFilt, response_data);
    response_data_filtered = filtfilt(bpFilt, response_data_filtered);

    % Segment data into trials
    response_trials = zeros(length(rising_edges), window_length);
    for trial = 1:length(rising_edges)
        stim_idx = rising_edges(trial);
        start_idx = max(stim_idx - samples_prestim, 1);
        end_idx = min(stim_idx + samples_poststim - 1, length(response_data_filtered));
        if (end_idx - start_idx + 1) == window_length
            response_trials(trial, :) = response_data_filtered(start_idx:end_idx);
        else
            response_trials(trial, :) = NaN;
        end
    end
    response_trials = response_trials(~any(isnan(response_trials), 2), :);
    average_response = mean(response_trials, 1);

    % Plot
    hold on;
    for trial = 1:size(response_trials, 1)
        plot(time_axis, response_trials(trial, :), 'Color', [0.8, 0.8, 0.8]); % Grey individual trials
    end
    plot(time_axis, average_response, 'k', 'LineWidth', 2); % Black trial average
    hold off;
    axis off;
    title(['Channel ', num2str(ch)]);
end
%% Plot 6: Number of trials and plot of stimulation waveform
% Define the timeframe
start_time = 40; % Start of the timeframe in seconds
end_time = 60;   % End of the timeframe in seconds

% Find trials within the timeframe
trials_in_timeframe = rising_edges(appended_time(rising_edges) >= start_time & appended_time(rising_edges) <= end_time);

% Count the number of trials
num_trials = length(trials_in_timeframe);

% Display the result
fprintf('Number of stimulation trials between %.1f s and %.1f s: %d\n', start_time, end_time, num_trials);

% Define the timeframe
start_time = 40; % Start of the timeframe in seconds
end_time = 60;   % End of the timeframe in seconds

% Mask for the selected timeframe
time_mask = appended_time >= start_time & appended_time <= end_time;

% Plot the stimulation waveform
figure('Name', 'Stimulation Waveform (40-60 s)', 'NumberTitle', 'off');
plot(appended_time(time_mask), appended_stim_waveform(time_mask), 'b'); hold on;

% Overlay markers for trials within the timeframe
trials_in_timeframe = rising_edges(appended_time(rising_edges) >= start_time & appended_time(rising_edges) <= end_time);
plot(appended_time(trials_in_timeframe), appended_stim_waveform(trials_in_timeframe), 'ro', 'MarkerSize', 5, 'LineWidth', 1.5);

% Add labels and title
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Stimulation Waveform (40-60 s)');
grid on;
legend('Stimulation Waveform', 'Stimulation Trials');

%% Plot 7: Stimulation waveform
figure('Name', 'Stimulation Waveform', 'NumberTitle', 'off');
plot(appended_time, appended_stim_waveform, 'b'); hold on;
plot(appended_time(rising_edges), appended_stim_waveform(rising_edges), 'ro', 'MarkerSize', 5, 'LineWidth', 1.5);
legend('Stimulation Waveform', ['Trials: ', num2str(length(rising_edges))]);
xlabel('Time (s)'); ylabel('Voltage (V)');
title('Stimulation Waveform and Trials'); grid on;

%% Plot 8: Baseline activity 200 uV offset
figure('Name', 'Minimal Offset Baseline Activity', 'NumberTitle', 'off');
hold on;

% Adjusted offsets
offsets = [0, -200, -400, -600]; % 200 µV step between signals

% Plot Raw and Filtered Channel 5
plot(time_selected, raw_channel5 + offsets(1), 'Color', colors{1}, 'LineWidth', 1.5, 'DisplayName', 'Raw Ch5');
plot(time_selected, filtered_channel5 + offsets(2), 'Color', colors{1}, 'LineWidth', 1.5, 'DisplayName', 'Filtered Ch5');

% Plot Raw and Filtered Channel 2
plot(time_selected, raw_channel2 + offsets(3), 'Color', colors{2}, 'LineWidth', 1.5, 'DisplayName', 'Raw Ch2');
plot(time_selected, filtered_channel2 + offsets(4), 'Color', colors{2}, 'LineWidth', 1.5, 'DisplayName', 'Filtered Ch2');

hold off;

% Add Axes and Labels
xlabel('Time (s)', 'FontSize', 12); % X-axis label in seconds
ylabel('Amplitude (\muV)', 'FontSize', 12); % Y-axis label in microvolts
title('Baseline Activity with 200 µV Offsets', 'FontSize', 14); % Add title
grid on; % Add grid for better visualization
legend show; % Display legend

% Set Axis Limits
xlim([time_selected(1), time_selected(end)]); % X-axis range based on selected time
ylim([min(offsets)-100, max(offsets)+200]); % Y-axis to include offsets

% Save as EPS
file_name = 'Baseline_Offset_Activity_200uV';
set(gcf, 'Renderer', 'painters'); % Use painters renderer for vector output
print(gcf, [file_name, '.eps'], '-depsc', '-painters'); % Save as EPS

% Uncomment the following line to save as SVG if needed:
% print(gcf, [file_name, '.svg'], '-dsvg', '-painters');

fprintf('Baseline Offset Activity plot exported successfully as EPS file.\n');

%% Plot 9: Subplots evoked response TIFF with labels
% Channels to Plot
channels_to_plot = [2, 3, 5, 6];

% Iterate Over Each Channel
for idx = 1:length(channels_to_plot)
    ch = channels_to_plot(idx);
    
    % Extract and Filter Data for the Current Channel
    response_data = appended_data(ch, :);
    response_data_filtered = filtfilt(notchFilt, response_data);
    response_data_filtered = filtfilt(bpFilt, response_data_filtered);

    % Segment Data into Trials
    response_trials = zeros(length(rising_edges), window_length);
    for trial = 1:length(rising_edges)
        stim_idx = rising_edges(trial);
        start_idx = max(stim_idx - samples_prestim, 1);
        end_idx = min(stim_idx + samples_poststim - 1, length(response_data_filtered));
        if (end_idx - start_idx + 1) == window_length
            response_trials(trial, :) = response_data_filtered(start_idx:end_idx);
        else
            response_trials(trial, :) = NaN;
        end
    end
    response_trials = response_trials(~any(isnan(response_trials), 2), :);
    average_response = mean(response_trials, 1);

    % Create a New Figure for Each Channel
    figure('Name', ['Evoked Response - Channel ', num2str(ch)], 'NumberTitle', 'off');
    hold on;
    
    % Plot Individual Trials in Grey
    for trial = 1:size(response_trials, 1)
        plot(time_axis, response_trials(trial, :), 'Color', [0.8, 0.8, 0.8]); % Grey lines
    end
    
    % Plot the Average Response in Black
    plot(time_axis, average_response, 'k', 'LineWidth', 2); % Black line
    
    % Plot the Stimulus Marker
    plot(0, 0, 'ro', 'MarkerSize', 8, 'LineWidth', 2); % Red dot at t=0
    
    hold off;
    
    % Add Axes and Units
    xlabel('Time (ms)', 'FontSize', 12); % X-axis label in milliseconds
    ylabel('Amplitude (\muV)', 'FontSize', 12); % Y-axis label in microvolts
    title(['Evoked Response - Channel ', num2str(ch)], 'FontSize', 14); % Add channel title
    xlim([-10, 100]); % Set x-axis limits to match time axis
    ylim([-100, 100]); % Set y-axis limits to standard amplitude range
    grid on; % Add grid for better readability
    
    % Save Plot as TIFF
    file_base = ['Channel_', num2str(ch), '_Evoked_Response'];
    print(gcf, [file_base, '.tiff'], '-dtiff', '-r600'); % Save as TIFF with 600 DPI
    
    fprintf('Channel %d exported successfully as TIFF file.\n', ch);
end

%% Plot 10: Minimal design Subplots evoked response TIFF 
% Channels to Plot
channels_to_plot = [2, 3, 5, 6];

% Iterate Over Each Channel
for idx = 1:length(channels_to_plot)
    ch = channels_to_plot(idx);
    
    % Extract and Filter Data for the Current Channel
    response_data = appended_data(ch, :);
    response_data_filtered = filtfilt(notchFilt, response_data);
    response_data_filtered = filtfilt(bpFilt, response_data_filtered);

    % Segment Data into Trials
    response_trials = zeros(length(rising_edges), window_length);
    for trial = 1:length(rising_edges)
        stim_idx = rising_edges(trial);
        start_idx = max(stim_idx - samples_prestim, 1);
        end_idx = min(stim_idx + samples_poststim - 1, length(response_data_filtered));
        if (end_idx - start_idx + 1) == window_length
            response_trials(trial, :) = response_data_filtered(start_idx:end_idx);
        else
            response_trials(trial, :) = NaN;
        end
    end
    response_trials = response_trials(~any(isnan(response_trials), 2), :);
    average_response = mean(response_trials, 1);

    % Create a New Figure for Each Channel
    figure('Name', ['Minimal Evoked Response - Channel ', num2str(ch)], 'NumberTitle', 'off');
    hold on;
    
    % Plot Individual Trials in Grey
    for trial = 1:size(response_trials, 1)
        plot(time_axis, response_trials(trial, :), 'Color', [0.8, 0.8, 0.8]); % Grey lines
    end
    
    % Plot the Average Response in Black
    plot(time_axis, average_response, 'k', 'LineWidth', 2); % Black line
    
    hold off;
    
    % Minimal Design
    axis off; % Remove axes for a minimalistic look
    title(['Channel ', num2str(ch)], 'FontSize', 12); % Add minimal title for identification
    
    % Save Plot as TIFF
    file_base = ['Channel_', num2str(ch), '_Minimal_Evoked_Response'];
    print(gcf, [file_base, '.tiff'], '-dtiff', '-r600'); % Save as TIFF with 600 DPI
    
    fprintf('Channel %d exported successfully as minimal TIFF file.\n', ch);
end

%% Plot 11: Stimulation waveform showing from 32 to 35 seconds
figure('Name', 'Stimulation Waveform - 32 to 35 Seconds', 'NumberTitle', 'off');

% Define the time range for visualization (32 to 35 seconds)
time_mask = appended_time >= 32 & appended_time <= 35;

% Plot the stimulation waveform
plot(appended_time(time_mask), appended_stim_waveform(time_mask), 'b'); hold on;

% Overlay markers for trials within the range
trials_within_range = rising_edges(appended_time(rising_edges) >= 32 & appended_time(rising_edges) <= 35);
plot(appended_time(trials_within_range), appended_stim_waveform(trials_within_range), 'ro', 'MarkerSize', 5, 'LineWidth', 1.5);

% Add legend and labels
legend('Stimulation Waveform', ['Trials: ', num2str(length(trials_within_range))], 'Location', 'best');
xlabel('Time (s)', 'FontSize', 12); % X-axis label in seconds
ylabel('Voltage (V)', 'FontSize', 12); % Y-axis label in volts
title('Stimulation Waveform - 32 to 35 Seconds', 'FontSize', 14); % Title
grid on;

% Save Plot as EPS
file_name = 'Stimulation_Waveform_32_to_35_Seconds';
set(gcf, 'Renderer', 'painters'); % Use painters renderer for vector output
print(gcf, [file_name, '.eps'], '-depsc', '-painters'); % Save as EPS

fprintf('Stimulation waveform plot from 32 to 35 seconds exported successfully as EPS file.\n');

%% Plot 12: PSD plot

pre_start = -0.1; % Baseline start time (e.g., -100 ms relative to stimulus)
pre_end = 0;      % Baseline end time (e.g., 0 ms)
post_start = 0;   % Post-stimulus start time (e.g., 0 ms)
post_end = 0.1;   % Post-stimulus end time (e.g., 100 ms)

% Convert timeframes to sample indices
pre_samples = round(pre_start * fs):round(pre_end * fs);
post_samples = round(post_start * fs):round(post_end * fs);

%% Extract Pre- and Post-Stimulus Data
pre_stim_data_all = [];
post_stim_data_all = [];

for trial_idx = 1:length(rising_edges)
    stim_idx = rising_edges(trial_idx);

    % Ensure indices stay within bounds
    if (stim_idx + pre_samples(1) > 0) && (stim_idx + post_samples(end) <= length(appended_time))
        % Extract pre- and post-stimulus data for each trial
        pre_stim_data_all = [pre_stim_data_all, appended_data(:, stim_idx + pre_samples)];
        post_stim_data_all = [post_stim_data_all, appended_data(:, stim_idx + post_samples)];
    end
end

%% Compute PSD for Channels 2 and 5
[psd_pre_ch2, f] = pwelch(pre_stim_data_all(2, :), [], [], [], fs); % Pre-stimulus PSD (Channel 2)
[psd_post_ch2, ~] = pwelch(post_stim_data_all(2, :), [], [], [], fs); % Post-stimulus PSD (Channel 2)

[psd_pre_ch5, ~] = pwelch(pre_stim_data_all(5, :), [], [], [], fs); % Pre-stimulus PSD (Channel 5)
[psd_post_ch5, ~] = pwelch(post_stim_data_all(5, :), [], [], [], fs); % Post-stimulus PSD (Channel 5)

%% Plot PSD (Pre- and Post-Stimulus for Channel 2 and 5)
figure;

% Channel 2
subplot(2, 1, 1);
plot(f, 10*log10(psd_pre_ch2), 'b', 'LineWidth', 2); hold on;
plot(f, 10*log10(psd_post_ch2), 'r', 'LineWidth', 2); hold off;
title('PSD - Channel 2 (Pre- vs Post-Stimulus)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
legend('Pre-Stimulus', 'Post-Stimulus');
xlim([0 100]); grid on;

% Channel 5
subplot(2, 1, 2);
plot(f, 10*log10(psd_pre_ch5), 'b', 'LineWidth', 2); hold on;
plot(f, 10*log10(psd_post_ch5), 'r', 'LineWidth', 2); hold off;
title('PSD - Channel 5 (Pre- vs Post-Stimulus)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
legend('Pre-Stimulus', 'Post-Stimulus');
xlim([0 100]); grid on;%% Compute PSD for Channels 2 and 5
[psd_pre_ch2, f] = pwelch(pre_stim_data_all(2, :), [], [], [], fs); % Pre-stimulus PSD (Channel 2)
[psd_post_ch2, ~] = pwelch(post_stim_data_all(2, :), [], [], [], fs); % Post-stimulus PSD (Channel 2)

[psd_pre_ch5, ~] = pwelch(pre_stim_data_all(5, :), [], [], [], fs); % Pre-stimulus PSD (Channel 5)
[psd_post_ch5, ~] = pwelch(post_stim_data_all(5, :), [], [], [], fs); % Post-stimulus PSD (Channel 5)

%% Plot PSD (Pre- and Post-Stimulus for Channel 2 and 5)
figure;

% Channel 2
subplot(2, 1, 1);
plot(f, 10*log10(psd_pre_ch2), 'b', 'LineWidth', 2); hold on;
plot(f, 10*log10(psd_post_ch2), 'r', 'LineWidth', 2); hold off;
title('PSD - Channel 2 (Pre- vs Post-Stimulus)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
legend('Pre-Stimulus', 'Post-Stimulus');
xlim([0 100]); grid on;

% Channel 5
subplot(2, 1, 2);
plot(f, 10*log10(psd_pre_ch5), 'b', 'LineWidth', 2); hold on;
plot(f, 10*log10(psd_post_ch5), 'r', 'LineWidth', 2); hold off;
title('PSD - Channel 5 (Pre- vs Post-Stimulus)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
legend('Pre-Stimulus', 'Post-Stimulus');
xlim([0 100]); grid on;

%% Save Pre- and Post-Stimulus PSD Results for Later Use
save('pre_post_stim_psd.mat', 'psd_pre_ch2', 'psd_post_ch2', 'psd_pre_ch5', 'psd_post_ch5', 'f');
fprintf('Pre- and post-stimulus PSD data saved.\n');

%% Plot 6: Combined PSD Channels 2 & 3 and Channels 5 & 6
% Extract pre- and post-stimulus data for each group
pre_stim_ch23 = mean(pre_stim_data_all([2, 3], :), 1); % Average Channels 2 & 3 pre-stimulus
post_stim_ch23 = mean(post_stim_data_all([2, 3], :), 1); % Average Channels 2 & 3 post-stimulus

pre_stim_ch56 = mean(pre_stim_data_all([5, 6], :), 1); % Average Channels 5 & 6 pre-stimulus
post_stim_ch56 = mean(post_stim_data_all([5, 6], :), 1); % Average Channels 5 & 6 post-stimulus

%% Compute PSD for Combined Groups
[psd_pre_ch23, f] = pwelch(pre_stim_ch23, [], [], [], fs); % PSD for pre-stimulus Channels 2 & 3
[psd_post_ch23, ~] = pwelch(post_stim_ch23, [], [], [], fs); % PSD for post-stimulus Channels 2 & 3

[psd_pre_ch56, ~] = pwelch(pre_stim_ch56, [], [], [], fs); % PSD for pre-stimulus Channels 5 & 6
[psd_post_ch56, ~] = pwelch(post_stim_ch56, [], [], [], fs); % PSD for post-stimulus Channels 5 & 6

%% Plot PSD for Combined Channels 2 & 3 and 5 & 6
figure;

% Combined Channels 2 & 3
subplot(2, 1, 1);
plot(f, 10*log10(psd_pre_ch23), 'b', 'LineWidth', 2); hold on;
plot(f, 10*log10(psd_post_ch23), 'r', 'LineWidth', 2); hold off;
title('PSD - Combined Channels 2 & 3 (Pre- vs Post-Stimulus)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
legend('Pre-Stimulus', 'Post-Stimulus');
xlim([0 100]); grid on;

% Combined Channels 5 & 6
subplot(2, 1, 2);
plot(f, 10*log10(psd_pre_ch56), 'b', 'LineWidth', 2); hold on;
plot(f, 10*log10(psd_post_ch56), 'r', 'LineWidth', 2); hold off;
title('PSD - Combined Channels 5 & 6 (Pre- vs Post-Stimulus)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
legend('Pre-Stimulus', 'Post-Stimulus');
xlim([0 100]); grid on;

%% Save Results
save('combined_channels_psd.mat', 'psd_pre_ch23', 'psd_post_ch23', 'psd_pre_ch56', 'psd_post_ch56', 'f');
fprintf('Combined PSD for Channels 2 & 3 and 5 & 6 saved.\n');

%% Plot 13: PSD for Combined Post and Pre
% Extract pre- and post-stimulus data for each group
pre_stim_ch23 = mean(pre_stim_data_all([2, 3], :), 1); % Average Channels 2 & 3 pre-stimulus
post_stim_ch23 = mean(post_stim_data_all([2, 3], :), 1); % Average Channels 2 & 3 post-stimulus

pre_stim_ch56 = mean(pre_stim_data_all([5, 6], :), 1); % Average Channels 5 & 6 pre-stimulus
post_stim_ch56 = mean(post_stim_data_all([5, 6], :), 1); % Average Channels 5 & 6 post-stimulus

% Compute PSD for combined groups
[psd_pre_ch23, f] = pwelch(pre_stim_ch23, [], [], [], fs); % PSD for pre-stimulus Channels 2 & 3
[psd_post_ch23, ~] = pwelch(post_stim_ch23, [], [], [], fs); % PSD for post-stimulus Channels 2 & 3

[psd_pre_ch56, ~] = pwelch(pre_stim_ch56, [], [], [], fs); % PSD for pre-stimulus Channels 5 & 6
[psd_post_ch56, ~] = pwelch(post_stim_ch56, [], [], [], fs); % PSD for post-stimulus Channels 5 & 6

%% Plot Combined Pre-Stimulus PSD
figure;
plot(f, 10*log10(psd_pre_ch23), 'b', 'LineWidth', 2); hold on;
plot(f, 10*log10(psd_pre_ch56), 'g', 'LineWidth', 2); hold off;
title('Pre-Stimulus PSD (Combined Channels)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
legend('Channels 2 & 3', 'Channels 5 & 6');
xlim([0 80]); % Frequency range up to 80 Hz
grid on;

% Save Plot as EPS
file_name = 'Pre_Stimulus_PSD_Combined_Channels';
set(gcf, 'Renderer', 'painters'); % Use painters renderer for vector output
print(gcf, [file_name, '.eps'], '-depsc', '-painters'); % Save as EPS
fprintf('Pre-stimulus PSD plot exported successfully as EPS file.\n');

%% Plot Combined Post-Stimulus PSD
figure;
plot(f, 10*log10(psd_post_ch23), 'r', 'LineWidth', 2); hold on;
plot(f, 10*log10(psd_post_ch56), 'm', 'LineWidth', 2); hold off;
title('Post-Stimulus PSD (Combined Channels)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
legend('Channels 2 & 3', 'Channels 5 & 6');
xlim([0 80]); % Frequency range up to 80 Hz
grid on;

% Save Plot as EPS
file_name = 'Post_Stimulus_PSD_Combined_Channels';
set(gcf, 'Renderer', 'painters'); % Use painters renderer for vector output
print(gcf, [file_name, '.eps'], '-depsc', '-painters'); % Save as EPS
fprintf('Post-stimulus PSD plot exported successfully as EPS file.\n');

%% Save Results
save('combined_channels_psd.mat', 'psd_pre_ch23', 'psd_post_ch23', 'psd_pre_ch56', 'psd_post_ch56', 'f');
fprintf('Combined PSD for Channels 2 & 3 and 5 & 6 saved.\n');

% Numerical 1: SNR 
%% Define Parameters
pre_start = -0.1; % Pre-stimulus start time (s)
pre_end = 0;      % Pre-stimulus end time (s)
post_start = 0;   % Post-stimulus start time (s)
post_end = 0.1;   % Post-stimulus end time (s)
fs = frequency_parameters.amplifier_sample_rate; % Sampling rate

% Convert to sample indices
pre_samples = round(pre_start * fs):round(pre_end * fs);
post_samples = round(post_start * fs):round(post_end * fs);

%% Initialize SNR Arrays
snr_ch2_ch3 = []; % SNR for combined Channels 2 & 3
snr_ch5_ch6 = []; % SNR for combined Channels 5 & 6

%% Process Trials
for trial_idx = 1:length(rising_edges)
    stim_idx = rising_edges(trial_idx);
    
    % Ensure indices are within bounds
    if (stim_idx + pre_samples(1) > 0) && (stim_idx + post_samples(end) <= length(appended_time))
        % Combine Channels 2 & 3
        pre_data_2_3 = mean(appended_data([2, 3], stim_idx + pre_samples), 1);
        post_data_2_3 = mean(appended_data([2, 3], stim_idx + post_samples), 1);
        
        % Combine Channels 5 & 6
        pre_data_5_6 = mean(appended_data([5, 6], stim_idx + pre_samples), 1);
        post_data_5_6 = mean(appended_data([5, 6], stim_idx + post_samples), 1);
        
        % Compute PSD
        [psd_pre_2_3, f] = pwelch(pre_data_2_3, [], [], [], fs);
        [psd_post_2_3, ~] = pwelch(post_data_2_3, [], [], [], fs);
        
        [psd_pre_5_6, ~] = pwelch(pre_data_5_6, [], [], [], fs);
        [psd_post_5_6, ~] = pwelch(post_data_5_6, [], [], [], fs);
        
        % Compute AUC (0–100 Hz)
        auc_pre_2_3 = trapz(f(f <= 100), psd_pre_2_3(f <= 100));
        auc_post_2_3 = trapz(f(f <= 100), psd_post_2_3(f <= 100));
        
        auc_pre_5_6 = trapz(f(f <= 100), psd_pre_5_6(f <= 100));
        auc_post_5_6 = trapz(f(f <= 100), psd_post_5_6(f <= 100));
        
        % Compute SNR
        snr_trial_2_3 = 10 * log10(auc_post_2_3 / auc_pre_2_3);
        snr_trial_5_6 = 10 * log10(auc_post_5_6 / auc_pre_5_6);
        
        % Store Results
        snr_ch2_ch3 = [snr_ch2_ch3, snr_trial_2_3];
        snr_ch5_ch6 = [snr_ch5_ch6, snr_trial_5_6];
    end
end

%% Compute Mean and SD of SNR
mean_snr_2_3 = mean(snr_ch2_ch3);
sd_snr_2_3 = std(snr_ch2_ch3);

mean_snr_5_6 = mean(snr_ch5_ch6);
sd_snr_5_6 = std(snr_ch5_ch6);

%% Display Results
fprintf('Channel 2 & 3 - SNR: %.2f ± %.2f dB\n', mean_snr_2_3, sd_snr_2_3);
fprintf('Channel 5 & 6 - SNR: %.2f ± %.2f dB\n', mean_snr_5_6, sd_snr_5_6);

%% Numerical 2: SNR <30 Hz
%% Define Parameters
pre_start = -0.1; % Pre-stimulus start time (s)
pre_end = 0;      % Pre-stimulus end time (s)
post_start = 0;   % Post-stimulus start time (s)
post_end = 0.1;   % Post-stimulus end time (s)
fs = frequency_parameters.amplifier_sample_rate; % Sampling rate

% Convert to sample indices
pre_samples = round(pre_start * fs):round(pre_end * fs);
post_samples = round(post_start * fs):round(post_end * fs);

%% Initialize SNR Arrays
snr_ch2_ch3 = []; % SNR for combined Channels 2 & 3
snr_ch5_ch6 = []; % SNR for combined Channels 5 & 6

%% Process Trials
for trial_idx = 1:length(rising_edges)
    stim_idx = rising_edges(trial_idx);
    
    % Ensure indices are within bounds
    if (stim_idx + pre_samples(1) > 0) && (stim_idx + post_samples(end) <= length(appended_time))
        % Combine Channels 2 & 3
        pre_data_2_3 = mean(appended_data([2, 3], stim_idx + pre_samples), 1);
        post_data_2_3 = mean(appended_data([2, 3], stim_idx + post_samples), 1);
        
        % Combine Channels 5 & 6
        pre_data_5_6 = mean(appended_data([5, 6], stim_idx + pre_samples), 1);
        post_data_5_6 = mean(appended_data([5, 6], stim_idx + post_samples), 1);
        
        % Compute PSD
        [psd_pre_2_3, f] = pwelch(pre_data_2_3, [], [], [], fs);
        [psd_post_2_3, ~] = pwelch(post_data_2_3, [], [], [], fs);
        
        [psd_pre_5_6, ~] = pwelch(pre_data_5_6, [], [], [], fs);
        [psd_post_5_6, ~] = pwelch(post_data_5_6, [], [], [], fs);
        
        % Compute AUC (<30 Hz)
        auc_pre_2_3 = trapz(f(f < 30), psd_pre_2_3(f < 30));
        auc_post_2_3 = trapz(f(f < 30), psd_post_2_3(f < 30));
        
        auc_pre_5_6 = trapz(f(f < 30), psd_pre_5_6(f < 30));
        auc_post_5_6 = trapz(f(f < 30), psd_post_5_6(f < 30));
        
        % Compute SNR
        snr_trial_2_3 = 10 * log10(auc_post_2_3 / auc_pre_2_3);
        snr_trial_5_6 = 10 * log10(auc_post_5_6 / auc_pre_5_6);
        
        % Store Results
        snr_ch2_ch3 = [snr_ch2_ch3, snr_trial_2_3];
        snr_ch5_ch6 = [snr_ch5_ch6, snr_trial_5_6];
    end
end

%% Compute Mean and SD of SNR
mean_snr_2_3 = mean(snr_ch2_ch3);
sd_snr_2_3 = std(snr_ch2_ch3);

mean_snr_5_6 = mean(snr_ch5_ch6);
sd_snr_5_6 = std(snr_ch5_ch6);

%% Display Results
fprintf('Channel 2 & 3 - SNR (<30 Hz): %.2f ± %.2f dB\n', mean_snr_2_3, sd_snr_2_3);
fprintf('Channel 5 & 6 - SNR (<30 Hz): %.2f ± %.2f dB\n', mean_snr_5_6, sd_snr_5_6);


%% 20–80% LATENCY 
channels_lat = [2 3 5 6];
latencies_ms = NaN(size(channels_lat));

for iC = 1:length(channels_lat)
    ch = channels_lat(iC);

    % Filter channel
    resp = appended_data(ch,:);
    resp_f = filtfilt(notchFilt, resp);
    resp_f = filtfilt(bpFilt,   resp_f);

    % Collect trials using your existing window definition
    all_trials = [];

    for tr = 1:length(rising_edges)
        stim_idx = rising_edges(tr);
        s1 = stim_idx - samples_prestim;
        s2 = stim_idx + samples_poststim - 1;

        if s1 >= 1 && s2 <= length(resp_f)
            all_trials(end+1, :) = resp_f(s1:s2);
        end
    end

    if isempty(all_trials)
        warning('No valid trials in channel %d', ch);
        continue;
    end

    % Trial-averaged waveform
    avg_seg = mean(all_trials, 1);

    % Post-stim segment (what your previous code used)
    seg_post = avg_seg(samples_prestim+1:end);
    t_post   = time_axis(samples_prestim+1:end);   % ms

    % Restrict latency search window (validated)
    mask = (t_post >= 5) & (t_post <= 30);
    tE   = t_post(mask);
    yE   = seg_post(mask);

    % Detect early peak
    if ~isempty(tE)
        prom = max(1, 0.15 * range(yE));
        [pk, locs] = findpeaks(yE, tE, ...
            'MinPeakProminence', prom, ...
            'MinPeakDistance', 1.5);

        if ~isempty(pk)
            peak_val = pk(1);
            peak_t   = locs(1);
        else
            [peak_val, idxp] = max(yE);
            peak_t = tE(idxp);
        end
    else
        [peak_val, idxp] = max(seg_post);
        peak_t = t_post(idxp);
    end

    % Peak index relative to full t_post
    [~, pidx] = min(abs(t_post - peak_t));

    % Thresholds
    y20 = 0.2 * peak_val;
    y80 = 0.8 * peak_val;

    % 80% point: last index BEFORE peak where y <= 0.8*peak
    idx80 = find(seg_post(1:pidx) <= y80, 1, 'last');

    % 20%: last index before idx80 where y <= 0.2*peak
    if ~isempty(idx80)
        idx20 = find(seg_post(1:idx80) <= y20, 1, 'last');
    else
        idx20 = [];
    end

    % Compute latency
    if ~isempty(idx20) && ~isempty(idx80) && idx80 > idx20
        xfit = t_post(idx20:idx80);
        yfit = seg_post(idx20:idx80);

        p = polyfit(xfit, yfit, 1);
        latency_est = -p(2)/p(1);
        latencies_ms(iC) = latency_est;
    end

    %% ---- Plot for this channel ----
    figure('Name',sprintf('Latency Channel %d',ch), 'Color','w');
    plot(time_axis, avg_seg,'k','LineWidth',2); hold on;
    xline(0,'r--');

    if ~isnan(latencies_ms(iC))
        plot(t_post(idx20), seg_post(idx20),'bo','MarkerSize',7,'LineWidth',1.5);
        plot(t_post(idx80), seg_post(idx80),'bo','MarkerSize',7,'LineWidth',1.5);
        plot(xfit, polyval(p,xfit),'b--','LineWidth',1.3);
        plot(latency_est, 0,'ro','MarkerSize',8,'LineWidth',1.5);
        title(sprintf('Channel %d Latency = %.2f ms', ch, latency_est));
    else
        title(sprintf('Channel %d Latency = NaN', ch));
    end

    xlabel('Time (ms)');
    ylabel('Amplitude (µV)');
    xlim([-5 40]); grid on;
end


%   Print Mean ± SD for groups

lat_23 = latencies_ms(1:2);
lat_56 = latencies_ms(3:4);

mean_23 = mean(lat_23,'omitnan');
sd_23   = std(lat_23,'omitnan');

mean_56 = mean(lat_56,'omitnan');
sd_56   = std(lat_56,'omitnan');

fprintf('\nChannels 2 & 3: %.2f ± %.2f ms\n', mean_23, sd_23);
fprintf('Channels 5 & 6: %.2f ± %.2f ms\n', mean_56, sd_56);
