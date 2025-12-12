% Load Wireless EEG Data 
folder_name = 'C:'; %% add your path

if ~exist(folder_name, 'dir')
    error('Folder does not exist: %s', folder_name);
end

csv_files = {'ch1.csv', 'ch2.csv', 'ch3.csv', 'ch4.csv'};
channels = cell(1, length(csv_files));
channel_lengths = zeros(1, length(csv_files));

for i = 1:length(csv_files)
    file_path = fullfile(folder_name, csv_files{i});
    if ~exist(file_path, 'file')
        error('File does not exist: %s', file_path);
    end
    data = table2array(readtable(file_path));
    if size(data, 2) < 2
        error('File %s has fewer than 2 columns', file_path);
    end
    channels{i} = data(:, 2) * 1e6;  % Convert to µV
    channel_lengths(i) = length(channels{i});
    fprintf('Loaded %s: %d samples\n', csv_files{i}, channel_lengths(i));
end

min_samples = min(channel_lengths);
fprintf('Minimum number of samples across channels: %d\n', min_samples);
for i = 1:length(channels)
    channels{i} = channels{i}(1:min_samples);
end
[ch1, ch2, ch3, ch4] = deal(channels{:});

% === Initialize Parameters ===
fs = 1000;
combined_amplifier_data = [ch1'; ch2'; ch3'; ch4'];
num_channels = size(combined_amplifier_data, 1);
N = size(combined_amplifier_data, 2);
fprintf('Concatenation complete: %d channels x %d samples\n', num_channels, N);

% === Filter Data (0.5–100 Hz bandpass + notch) ===
N_filt = 4;
d = designfilt('bandpassiir', 'FilterOrder', N_filt, ...
    'HalfPowerFrequency1', 0.5, ...
    'HalfPowerFrequency2', 100, ...
    'SampleRate', fs);

amplifier_data_filtered = zeros(num_channels, N);
for ch = 1:num_channels
    data = combined_amplifier_data(ch, :)';
    for nf = [60, 120, 180, 240, 300, 360, 420, 480]
        Wo = nf / (fs / 2); BW = Wo / 35;
        [b, a] = iirnotch(Wo, BW);
        data = filtfilt(b, a, data);
    end
    amplifier_data_filtered(ch, :) = filtfilt(d, data)';
end

% === Define Segments for Channel 4 ===
ch4_data = amplifier_data_filtered(3, :);
samples_per_minute = fs * 60;

idx_base = (2 * samples_per_minute + 1):(3 * samples_per_minute);  % 2:00–3:00
idx_seiz = (23 * samples_per_minute + 1):(24 * samples_per_minute); % 4:49–4:50

data_base = ch4_data(idx_base);
data_seiz = ch4_data(idx_seiz);
t_base = (0:(length(data_base)-1)) / fs;
t_seiz = (0:(length(data_seiz)-1)) / fs;

%% Plot 1: Baseline (2:00–3:00)
figure('Name', 'Baseline (2:00–3:00) — Channel 3', 'Color', 'w');
plot(t_base, data_base, 'b');
xlabel('Time (s)'); ylabel('Voltage (µV)');
ylim([-2000 2000]); title('Baseline (2:00–3:00) — Channel 4');

% Export the stacked EEG figure as SVG 
print('Zoomed_EEG_baseline', '-dsvg');

% Refine Data for Scalograms (1–300 Hz + notch) 
d_wavelet = designfilt('bandpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency1', 1, ...
    'HalfPowerFrequency2', 300, ...
    'SampleRate', fs);

% Filter baseline
data_base_filt = data_base;
for nf = [60, 120, 180, 240, 300, 360, 420, 480]
    [b, a] = iirnotch(nf / (fs / 2), nf / (fs / 2) / 35);
    data_base_filt = filtfilt(b, a, data_base_filt);
end
data_base_filt = filtfilt(d_wavelet, data_base_filt);

% Filter seizure
data_seiz_filt = data_seiz;
for nf = [60, 120, 180, 240, 300, 360, 420, 480]
    [b, a] = iirnotch(nf / (fs / 2), nf / (fs / 2) / 35);
    data_seiz_filt = filtfilt(b, a, data_seiz_filt);
end
data_seiz_filt = filtfilt(d_wavelet, data_seiz_filt);


%% Plot 2: Full seizure (3:00–5:00), Tight Offset, in mV 
idx_seiz = (18 * samples_per_minute + 1):(20 * samples_per_minute); % 3:00–5:00
t_seiz = (0:(length(idx_seiz)-1)) / fs;

offset = 10;  % Tighter vertical spacing

figure('Name', 'Tightly Stacked EEG — Channels 2–4 (3:00–5:00)', 'Color', 'w'); hold on;

for i = 1:3  % Channels 2, 3, 4
    ch = i + 1;
    data_mV = amplifier_data_filtered(ch, idx_seiz) / 1000;
    plot(t_seiz, data_mV + (i - 1) * offset, 'k');
end
% === Export the stacked EEG figure as SVG ===
print('Stacked_EEG_Ch2to4_3mto5m', '-dsvg');
xlabel('Time (s)');
ylabel('Voltage (mV)');
ylim([-5 (2 * offset + 5)]);
title('Seizure activity');
set(gca, 'Box', 'off', 'TickDir', 'out');
print('Stacked_EEG_Ch2to4_3mto5m', '-dpng', '-r600');


%% Plot 3: Zoomed-In Plot: Channel 4 (20–60 s within 3:00–5:00 window) 
start_sec = 20;  % Relative to 3:00
end_sec = 60;
fs = 1000;

% Get full 3:00–5:00 segment index
idx_seiz = (18 * samples_per_minute + 1):(20 * samples_per_minute); 
data_full = amplifier_data_filtered(4, idx_seiz) / 1000;  % Channel 4, in mV

% Extract zoomed portion
idx_zoom = (start_sec * fs + 1):(end_sec * fs);
t_zoom = (0:(length(idx_zoom)-1)) / fs + start_sec;
data_zoom = data_full(idx_zoom);

% Plot
figure('Name', 'Channel 4 (20–60 s of 3:00–5:00)', 'Color', 'w');
plot(t_zoom, data_zoom, 'k');
xlabel('Time (s)');
ylabel('Voltage (mV)');
ylim([-5 5]);
title('Seizure activity 20-60 s');
set(gca, 'Box', 'off', 'TickDir', 'out');
% Export the stacked EEG figure as SVG 
print('Zoomed_EEG_Ch4_20s_to_60s', '-dsvg');

%% Plot 4: Extract Channel 4 from 3:00–5:00 (in mV) 
idx_seiz = (18 * samples_per_minute + 1):(20 * samples_per_minute);
data_full = amplifier_data_filtered(4, idx_seiz) / 1000;

fs = 1000;

% Baseline (0–10 s)
idx_base = (0 * fs + 1):(10 * fs);
data_baseline = data_full(idx_base);

% Seizure (22–32 s) 
idx_seizure = (22 * fs + 1):(32 * fs);
data_seizure = data_full(idx_seizure);

% Compute Scalogram for Seizure First (to set C-axis) 
[wt_seiz, f_seiz] = cwt(data_seizure, 'morse', fs);
wt_seiz = abs(wt_seiz);
max_c = max(wt_seiz(f_seiz <= 100), [], 'all');  % Max up to 100 Hz

%% Plot 5: Seizure Scalogram (22–32 s as 0–10 s) 
figure('Name', 'Scalogram — Seizure (Ch 4: 22–32 s)', 'Color', 'w');
cwt(data_seizure, 'morse', fs);
ylim([0 100]);
caxis([0 max_c]);
colormap turbo;
title('Seizure');
xticks(0:2:10);
xlim([0 10]);
colorbar;
print('Scalogram_Seizure_Ch4_22s_to_32s', '-dpng', '-r600');

%% Plot 6: Baseline Scalogram (0–10 s) 
figure('Name', 'Scalogram — Baseline (Ch 4: 0–10 s)', 'Color', 'w');
cwt(data_baseline, 'morse', fs);
ylim([0 100]);
caxis([0 max_c]);
colormap turbo;
title('Baseline');
xticks(0:2:10);
xlim([0 10]);
colorbar;
print('Scalogram_Baseline_Ch4_0s_to_10s', '-dpng', '-r600');

%%  Compute High-Resolution Scalogram for Seizure Segment (22–32 s as 0–10 s) 

% Perform CWT to compute max value for caxis normalization
[wt_seiz, f_seiz] = cwt(data_seizure, 'morse', fs, ...
                        'VoicesPerOctave', 48, ...
                        'TimeBandwidth', 60);
wt_seiz = abs(wt_seiz);

% Max value up to 100 Hz
max_c = max(wt_seiz(f_seiz <= 100), [], 'all');

%% Plot 7 High-Resolution Scalogram 
figure('Name', 'Scalogram — Seizure (Ch 4: 22–32 s)', ...
       'Color', 'w', 'Position', [100 100 1200 500]);

cwt(data_seizure, 'morse', fs, ...
    'VoicesPerOctave', 48, ...
    'TimeBandwidth', 60);

ylim([0 100]);           % Frequency axis: 0–100 Hz
caxis([0 max_c]);        % Color axis scaling
colormap(turbo);         % Turbo color scheme
colorbar;

title('Scalogram – Channel 4 Seizure (22–32 s)', 'FontSize', 10);
xlabel('Time (s)', 'FontSize', 8);
ylabel('Frequency (Hz)', 'FontSize', 8);
xticks(0:2:10);          % Relabel time axis
xlim([0 10]);

% Export as SVG 
print(gcf, 'Scalogram_Seizure_Ch4_22s_to_32s', '-dsvg');

% Seizure (22–24 s) 
idx_seizure = (22 * fs + 1):(24 * fs);
data_seizure = data_full(idx_seizure);

% Compute Scalogram for 2 s Seizure Segment
[wt_seiz, f_seiz] = cwt(data_seizure, 'morse', fs);
wt_seiz = abs(wt_seiz);
max_c = max(wt_seiz(f_seiz <= 100), [], 'all');  % Max up to 100 Hz

%% Plot 8: Seizure Scalogram (22–24 s as 0–2 s)
figure('Name', 'Scalogram — Seizure (Ch 4: 22–24 s)', 'Color', 'w', ...
       'Position', [100 100 1200 500]);
cwt(data_seizure, 'morse', fs, ...
    'VoicesPerOctave', 48, ...
    'TimeBandwidth', 60);
ylim([0 100]);
caxis([0 max_c]);
colormap turbo;
title('Scalogram — Channel 4 Seizure (22–24 s)');
xlabel('Time (s)');
xticks(0:0.5:2);
xlim([0 2]);
colorbar;

% === Export as high-res SVG ===
print('Scalogram_Seizure_Ch4_22s_to_24s', '-dsvg');
