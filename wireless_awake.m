% Load CSV data 
folder_name = ''; %% add your path

if ~exist(folder_name, 'dir')
    error('Folder does not exist: %s', folder_name);
end

% All 4 channels
csv_files = {'ch1.csv', 'ch2.csv', 'ch3.csv', 'ch4.csv'};

channels = cell(1, length(csv_files));
channel_lengths = zeros(1, length(csv_files));

%% Read files
for i = 1:length(csv_files)
    file_path = fullfile(folder_name, csv_files{i});
    
    if ~exist(file_path, 'file')
        error('File does not exist: %s', file_path);
    end
    
    data = table2array(readtable(file_path));
    
    if size(data,2) < 2
        error('File %s has fewer than 2 columns', file_path);
    end
    
    channels{i} = data(:,2) * 1e6;   % convert to µV
    channel_lengths(i) = length(channels{i});
    
    fprintf('Loaded %s: %d samples\n', csv_files{i}, channel_lengths(i));
end

%% Ensure all channels same length
min_samples = min(channel_lengths);
fprintf('Minimum number of samples: %d\n', min_samples);

for i = 1:4
    channels{i} = channels{i}(1:min_samples);
end

%% Combine into matrix 4×N
combined_amplifier_data = zeros(4, min_samples);
for i = 1:4
    combined_amplifier_data(i,:) = channels{i}';
end

%% Parameters
fs = 1000;
N = min_samples;
t_amplifier = 0:1/fs:(N-1)/fs;
num_channels = 4;

%% Filtering setup
N_filt = 4;
lf = 0.5;
hf = 100;

d = designfilt('bandpassiir', 'FilterOrder', N_filt, ...
    'HalfPowerFrequency1', lf, ...
    'HalfPowerFrequency2', hf, ...
    'SampleRate', fs);

amplifier_data_filtered = zeros(num_channels, N);

%% Apply notch + bandpass
for ch = 1:num_channels
    fprintf('Filtering channel %d/%d\n', ch, num_channels);
    
    data = combined_amplifier_data(ch,:)';
    
    % Notch filters
    for nf = [60 120 180 240 300 360 420 480]
        Wo = nf/(fs/2);
        BW = Wo/35;
        [b,a] = iirnotch(Wo, BW);
        data = filtfilt(b,a,data);
    end
    
    % Bandpass
    data_filt = filtfilt(d, data);
    amplifier_data_filtered(ch,:) = data_filt';
end


% BASELINE PLOTS


figure('Color','w');

for ch = 1:4
    subplot(4,1,ch);
    plot(t_amplifier, amplifier_data_filtered(ch,:), 'k');
    ylabel(sprintf('Ch%d (\\muV)', ch));
    grid on; box on;
    xlim([0 max(t_amplifier)]);
    
    if ch == 1
        title('Filtered Baseline — Channels 1–4');
    end
    if ch == 4
        xlabel('Time (s)');
    end
end
% Plot 0–5 seconds for Channels 1–4
t0 = 0;
t1 = 5;

idx = (t_amplifier >= t0) & (t_amplifier <= t1);

fig = figure('Color','w');

for ch = 1:4
    subplot(4,1,ch);
    plot(t_amplifier(idx), amplifier_data_filtered(ch, idx), 'k');

    ylim([-310 310]);                  % <<< SAME Y-LIMITS FOR ALL
    xlim([t0 t1]);

    ylabel(sprintf('Ch%d (\\muV)', ch));
    grid on; box on;

    if ch == 1
        title('Filtered Baseline — 0 to 5 seconds (Channels 1–4)');
    end
    if ch == 4
        xlabel('Time (s)');
    end
end

% Export as EPS 
exportgraphics(fig, 'baseline_0to5s_channels1to4.eps', 'ContentType', 'vector');
