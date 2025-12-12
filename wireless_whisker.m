% Define folder name
folder_name = ''; % add your path

% Check if folder exists
if ~exist(folder_name, 'dir')
    error('Folder does not exist: %s', folder_name);
end

% List of CSV files
csv_files = {'ch.1.csv', 'ch.2.csv', 'ch.3.csv', 'ch.4.csv'};

% Initialize cell array to store channel data
channels = cell(1, length(csv_files));
channel_lengths = zeros(1, length(csv_files));

% Read CSV files and extract second column
for i = 1:length(csv_files)
    file_path = fullfile(folder_name, csv_files{i});
    
    % Check if file exists
    if ~exist(file_path, 'file')
        error('File does not exist: %s', file_path);
    end
    
    try
        % Read table and convert to array
        data = table2array(readtable(file_path));
        
        % Check if data has at least 2 columns
        if size(data, 2) < 2
            error('File %s has fewer than 2 columns', file_path);
        end
        
      % Extract second column and convert to microvolts
       channels{i} = data(:, 2) * 1e6;

        
        % Store number of samples
        channel_lengths(i) = length(channels{i});
        
        fprintf('Loaded %s: %d samples\n', csv_files{i}, channel_lengths(i));
    catch e
        error('Failed to read %s: %s', file_path, e.message);
    end
end

% Find the minimum number of samples
min_samples = min(channel_lengths);
fprintf('Minimum number of samples across channels: %d\n', min_samples);

% Truncate all channels to the minimum length
for i = 1:length(channels)
    channels{i} = channels{i}(1:min_samples);
end

% Assign to individual variables (ch1, ch2, ch3, ch4)
ch1 = channels{1};
ch2 = channels{2};
ch3 = channels{3};
ch4 = channels{4};
%% Initialize parameters
fs = 1000; % Sampling frequency (Hz)
segment_duration = 1; % Duration of each segment (seconds)
segment_length = fs * segment_duration; % Number of samples per segment
segment_time = (0:1/fs:(segment_length-1)/fs); % Time vector for 1-second segment

% Concatenate horizontally (along columns/samples)
combined_amplifier_data = [ch1';ch2';ch3';ch4'];

%% Verify and display results
fprintf('Concatenation complete:\n');
fprintf('- Combined amplifier_data: %d channels x %d samples\n', ...
    size(combined_amplifier_data, 1), size(combined_amplifier_data, 2));

t_amplifier = 0:1/fs:(size(combined_amplifier_data,2)-1)/fs;
num_channels = size(combined_amplifier_data,1);
N = size(combined_amplifier_data,2); % Total number of samples

%% Filter the data
% Filter parameters
N_filt = 4; % Butterworth filter order
lf = [0.5]; % Lower bounds in Hz
hf = [100]; % Upper bounds in Hz
indexofhztofilter = 1; % Bandpass range

% Design bandpass filter
d = designfilt('bandpassiir', 'FilterOrder', N_filt, ...
    'HalfPowerFrequency1', lf(indexofhztofilter), ...
    'HalfPowerFrequency2', hf(indexofhztofilter), ...
    'SampleRate', fs);

% Initialize output variables
amplifier_data_filtered = zeros(num_channels, N);
amplifier_data_notch = zeros(num_channels, N);

% Apply bandpass and notch filters to each channel
for ch = 1:num_channels
    % Extract channel data
    data = combined_amplifier_data(ch, :)';
    fprintf('Filtering channel %d/%d\n', ch, num_channels);
    
    % Apply 60 Hz and harmonics notch filters
    for nf = [60, 120, 180, 240, 300, 360, 420, 480]
        Wo = nf / (fs / 2); % Normalized frequency
        BW = Wo / 35; % Bandwidth
        [b, a] = iirnotch(Wo, BW);
        data = filtfilt(b, a, data);
    end
    amplifier_data_notch(ch, :) = data';
    
    % Apply bandpass filter
    data_filtered = filtfilt(d, data);
    amplifier_data_filtered(ch, :) = data_filtered';
end

%% Segment data into 1-second windows
% Calculate number of complete 1-second segments
num_segments = floor(N / segment_length);
fprintf('Number of 1-second segments: %d\n', num_segments);

% Initialize cell array to store segments
segments = cell(num_segments, 1);

% Segment the filtered data
for i = 1:num_segments
    start_idx = (i-1) * segment_length + 1;
    end_idx = i * segment_length;
    segments{i} = amplifier_data_filtered(:, start_idx:end_idx); % Store each segment
end

%% Compute the average of all segments for each channel
average_segments = zeros(num_channels, segment_length); % Initialize average matrix

% Sum all segments
for i = 1:num_segments
    average_segments = average_segments + segments{i};
end

% Divide by number of segments to get the mean
average_segments = average_segments / num_segments;

%% Plot results
for ch = 1:num_channels
    figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.15, 0.4]);
    for k = 1:num_segments
        plot(segment_time, segments{k}(ch, :), 'Color', [0.4 0.4 0.4 0.2]);
        hold on;
    end
    plot(segment_time, average_segments(ch, :), 'Color', [0 0 0], 'LineWidth', 2);
    xlim([0 segment_duration]);
    xlabel('Time (s)');
    ylabel('Voltage (uV)');
    % ylim([-100 100]);
    title(sprintf('Channel %d: 1-Second Average', ch));
end
% %%
% figure();
% cwt(average_segments(3,:), 'morse', fs, 'VoicesPerOctave', 48, 'TimeBandwidth', 50);
% ylim([0 100]);    % Frequency limits
% clim([3 9]);    % Color scale limits

%% waterfall plot
figure();
offset_step = 30;  % Change as needed

for i = 1:size(segments, 1)
    offset = (i - 1) * offset_step;
    plot(segment_time, segments{i}(3, :) + offset);  % Add vertical shift
    hold on;
end

xlabel('Time (s)');
ylabel('Amplitude (µV) + Offset');
ylim('tight')




