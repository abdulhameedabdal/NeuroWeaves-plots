% Load CSV data 
% Define folder name
folder_name = ''; % add your path

% Check if folder exists
if ~exist(folder_name, 'dir')
    error('Folder does not exist: %s', folder_name);
end

% List of CSV files (four channels: ch1–ch4)
csv_files = {'ch1.csv', 'ch2.csv', 'ch3.csv', 'ch4.csv'};

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

% Assign to individual variables (ch1–ch4)
ch1 = channels{1};
ch2 = channels{2};
ch3 = channels{3};
ch4 = channels{4};

%% Initialize parameters
fs = 1000; % Sampling frequency (Hz)
segment_duration = 1; % Duration of each segment (seconds)
segment_length = fs * segment_duration; % Number of samples per segment
segment_time = (0:1/fs:(segment_length-1)/fs); % Time vector for 1-second segment

% Concatenate horizontally (channels x samples)
combined_amplifier_data = [ch1'; ch2'; ch3'; ch4'];

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

%% Plot 10-second traces of Channel 1 across the whole recording
% Uses amplifier_data_filtered(1,:) and t_amplifier

ch_idx          = 1;          % channel 1
x_all           = amplifier_data_filtered(ch_idx, :);
samples_per_seg = fs * 10;    % 10-second segments
total_samples   = length(x_all);
total_segments  = floor(total_samples / samples_per_seg);

segs_per_fig = 10;            % 10 subplots (10 segments) per figure

seg_counter = 0;
for s = 1:total_segments
    % Start a new figure every 10 segments
    if mod(s-1, segs_per_fig) == 0
        figure('Color','w');
        seg_counter = 0;
    end
    
    seg_counter = seg_counter + 1;
    idx_start   = (s-1)*samples_per_seg + 1;
    idx_end     = idx_start + samples_per_seg - 1;
    
    t_seg = t_amplifier(idx_start:idx_end);
    x_seg = x_all(idx_start:idx_end);
    
    subplot(segs_per_fig, 1, seg_counter);
    plot(t_seg, x_seg);
    xlim([t_seg(1) t_seg(end)]);
    ylabel('\muV');
    title(sprintf('Ch1: %.1f–%.1f s', t_seg(1), t_seg(end)));
    grid on;
    
    if seg_counter == segs_per_fig || s == total_segments
        xlabel('Time (s)');
    end
end

%% SPINDLE DETECTION FOR CHANNEL 1 (10–16 Hz) 

x_raw = amplifier_data_filtered(1, :);
t = t_amplifier;
fs = 1000;

% Restrict detection to 0–1300 s
idx_use = t >= 0 & t <= 1300;
x_raw = x_raw(idx_use);
t = t(idx_use);

% 10–16 Hz bandpass (for both detection + plotting)
d_sp = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',10,'HalfPowerFrequency2',16,'SampleRate',fs);
x_filt = filtfilt(d_sp, x_raw);

% Hilbert envelope
amp = abs(hilbert(x_filt));

% Smooth (40 ms)
amp_s = movmean(amp,40);

% Thresholds
mA = mean(amp_s);
TH_low  = 2.0*mA;
TH_high = 3.5*mA;

% Detect spindles
above = amp_s > TH_high;
pks = find(above);
spindles = [];
i = 1;

while i <= length(pks)
    pk = pks(i);

    on = pk;
    while on>1 && amp_s(on)>TH_low
        on = on-1;
    end

    off = pk;
    while off<length(amp_s) && amp_s(off)>TH_low
        off = off+1;
    end

    dur = (off-on)/fs;
    if dur>=0.4 && dur<=2.0
        spindles = [spindles; t(on) t(off) dur];
    end

    i = i + sum(pks>=on & pks<=off);
end

fprintf('Spindles detected (10–16 Hz, 0–1300 s): %d\n', size(spindles,1));

%% PLOT 10-SECOND SEGMENTS OF 10–16 Hz FILTERED TRACE

seg_len = 10;
samples_per_seg = seg_len * fs;
x_all = x_filt;
total_samples = length(x_all);
total_segments = floor(total_samples / samples_per_seg);

segs_per_fig = 10;
subplot_counter = 0;

for s = 1:total_segments

    if mod(s-1, segs_per_fig) == 0
        figure('Color','w');
        subplot_counter = 0;
    end

    subplot_counter = subplot_counter + 1;

    idx_start = (s-1)*samples_per_seg + 1;
    idx_end   = idx_start + samples_per_seg - 1;

    t_seg = t(idx_start:idx_end);
    x_seg = x_all(idx_start:idx_end);

    subplot(segs_per_fig, 1, subplot_counter);
    plot(t_seg, x_seg, 'k');
    hold on;

    if ~isempty(spindles)
        sp = spindles(spindles(:,1) >= t_seg(1) & spindles(:,2) <= t_seg(end), :);
    else
        sp = [];
    end

    for k = 1:size(sp,1)
        xs = sp(k,1);
        xe = sp(k,2);
        patch([xs xe xe xs], [min(x_seg) min(x_seg) max(x_seg) max(x_seg)], ...
            [0.9 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    end

    plot(t_seg, x_seg, 'k');
    hold off;

    title(sprintf('Ch1 (10–16 Hz)   %.1f–%.1f s', t_seg(1), t_seg(end)));
    ylabel('\muV');
    grid on;

    if subplot_counter == segs_per_fig || s == total_segments
        xlabel('Time (s)');
    end
end

%% Plot 170–180 s window for Ch1 (10–16 Hz filtered) 

t_start = 170;
t_end   = 180;

% Extract window
idx = t >= t_start & t <= t_end;
t_seg = t(idx);
x_seg = x_filt(idx);   % 10–16 Hz filtered Ch1

% Long, narrow SI figure (7.73 × 1.21 inches)
fig = figure('Color','w','Units','inches','Position',[1 1 7.73 1.21]);

plot(t_seg, x_seg, 'k', 'LineWidth', 1); hold on;

% Highlight spindles inside the window
sp = spindles(spindles(:,1) >= t_start & spindles(:,2) <= t_end, :);

for k = 1:size(sp,1)
    xs = sp(k,1);
    xe = sp(k,2);
    patch([xs xe xe xs], [-100 -100 100 100], ...
        [0.9 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
end

plot(t_seg, x_seg, 'k', 'LineWidth',1);  % redraw on top

xlabel('Time (s)', 'FontName','Helvetica','FontSize',9);
ylabel('\muV', 'FontName','Helvetica','FontSize',9);

% Fix title
title('Spindle detection filtered (10–16 Hz)', ...
      'FontName','Helvetica','FontSize',10);

ylim([-100 100]);
xlim([t_start t_end]);

set(gca,'FontName','Helvetica','FontSize',9, ...
        'LineWidth',1,'Box','off');
grid on;

% EXPORT AS TIFF
exportgraphics(fig,'SpindleDetection_10-16Hz_170-180s.tiff', ...
               'Resolution',300,'ContentType','image');


%% CROSS-VERIFY SPINDLES ON OTHER CHANNELS (10–16 Hz)

fs = 1000;

fprintf('\n=== CROSS-CHECKING OTHER CHANNELS (10–16 Hz) ===\n');
fprintf('Total spindles on Ch1: %d\n', size(spindles,1));

% rows = spindles, columns = channels 2..num_channels
ch_confirm = zeros(size(spindles,1), num_channels-1);

for ch = 2:num_channels
    x_ch = amplifier_data_filtered(ch, :);
    t_ch = t_amplifier;

    % same restriction 0–1300 s
    idx_use_ch = t_ch >= 0 & t_ch <= 1300;
    x_ch = x_ch(idx_use_ch);
    t_ch = t_ch(idx_use_ch);

    % 10–16 Hz bandpass (same as Ch1)
    d_sp_ch = designfilt('bandpassiir','FilterOrder',4, ...
        'HalfPowerFrequency1',10,'HalfPowerFrequency2',16,'SampleRate',fs);
    x_ch_filt = filtfilt(d_sp_ch, x_ch);

    % Hilbert envelope + smoothing
    amp_ch   = abs(hilbert(x_ch_filt));
    amp_ch_s = movmean(amp_ch,40);

    % Channel-specific thresholds (same relative levels)
    mA_ch      = mean(amp_ch_s);
    TH_low_ch  = 2.0*mA_ch;
    TH_high_ch = 3.5*mA_ch;

    for k = 1:size(spindles,1)
        t_start = spindles(k,1);
        t_end   = spindles(k,2);

        % convert times to indices in this channel vector
        idx_s = find(t_ch >= t_start, 1, 'first');
        idx_e = find(t_ch <= t_end,   1, 'last');

        if isempty(idx_s) || isempty(idx_e)
            continue
        end

        seg_amp = amp_ch_s(idx_s:idx_e);

        % spindle exists on this channel if segment exceeds TH_high_ch
        if any(seg_amp > TH_high_ch)
            ch_confirm(k, ch-1) = 1;
        end
    end

    fprintf('Channel %d also detected spindles for %d of the %d Ch1 events.\n', ...
            ch, sum(ch_confirm(:,ch-1)), size(spindles,1));
end

% Detailed indices for each channel
for ch = 2:num_channels
    fprintf('Spindles confirmed on Channel %d (indices into Ch1 spindle list):\n', ch);
    disp(find(ch_confirm(:,ch-1))');
end

%% PER-SPINDLE SUBPLOTS (BROADBAND, SINGLE COLOR, ±1 s, FIXED Y)
% Single patch color for all spindles (light green)

if isempty(spindles)
    warning('No spindles detected. Skipping broadband per-spindle plots.');
else
    fs = 1000;

    % Restrict broadband to same 0–1300 s window
    x_bb_full = amplifier_data_filtered(1, :);
    t_full    = t_amplifier;

    idx_use_bb = t_full >= 0 & t_full <= 1300;
    x_bb       = x_bb_full(idx_use_bb);
    t_bb       = t_full(idx_use_bb);

    % Window around spindle center (±1 s)
    win_pre  = 1.0;
    win_post = 1.0;
    n_pre    = round(win_pre  * fs);
    n_post   = round(win_post * fs);
    rel_t    = (-n_pre:n_post) / fs;

    seg_cell   = {};
    rel_on_all = [];
    rel_off_all= [];

    for k = 1:size(spindles,1)
        t_on  = spindles(k,1);
        t_off = spindles(k,2);
        t_c   = 0.5*(t_on + t_off);   % center

        % index of center in broadband vector
        [~, idx_c] = min(abs(t_bb - t_c));

        idx_start = idx_c - n_pre;
        idx_end   = idx_c + n_post;

        % skip if window exceeds bounds
        if idx_start < 1 || idx_end > numel(x_bb)
            continue;
        end

        seg_cell{end+1} = x_bb(idx_start:idx_end); %#ok<SAGROW>
        rel_on_all(end+1)  = t_on  - t_c;         %#ok<SAGROW>
        rel_off_all(end+1) = t_off - t_c;         %#ok<SAGROW>
    end

    n_valid = numel(seg_cell);
    fprintf('Broadband spindles (10–16 Hz defined), %d segments (±1 s window).\n', n_valid);

    if n_valid == 0
        warning('No full-window spindles available (broadband).');
    else
        % Subplot layout
        nrows = ceil(sqrt(n_valid));
        ncols = ceil(n_valid / nrows);

        figure('Color','w');
        for i = 1:n_valid
            subplot(nrows, ncols, i);
            seg_i = seg_cell{i};

            % Shade spindle interval (same color for all)
            xs = max(rel_on_all(i),  rel_t(1));
            xe = min(rel_off_all(i), rel_t(end));

            if xs < xe
                patch([xs xe xe xs], [-200 -200 200 200], ...
                      [0.8 1.0 0.8], 'FaceAlpha', 0.35, 'EdgeColor', 'none'); % light green
            end

            hold on;
            plot(rel_t, seg_i, 'k', 'LineWidth', 1);

            ylim([-200 200]);
            xlim([rel_t(1) rel_t(end)]);
            grid on;

            title(sprintf('Spindle %d', i), 'FontSize', 9);

            if i > n_valid - ncols
                xlabel('Time (s)');
            end
            ylabel('\muV');
        end
    end
end

%% SPINDLE-ASSOCIATED SLOW WAVES ONLY (Y-LIM FIXED TO ±150 µV)

fs    = 1000;
x_raw = amplifier_data_filtered(1,:);
t_raw = t_amplifier;

figure('Color','w');

nSp = size(spindles,1);
nrows = ceil(nSp / 5);
ncols = min(5, nSp);

for s = 1:nSp

    t_on  = spindles(s,1);
    t_off = spindles(s,2);

    % window (±0.5 s padding)
    idx   = t_raw >= (t_on - 0.5) & t_raw <= (t_off + 0.5);
    t_win = t_raw(idx);
    x_win = x_raw(idx);

    if numel(t_win) < 100
        continue
    end

    %%  DETECT TRUE SLOW WAVES (0.2–1.0 s NEG HALF-WAVE) 

    zc = find(diff(sign(x_win)) ~= 0);
    SWpeak = [];
    SWseg  = [];

    for i = 1:length(zc)-1

        z1 = zc(i);
        z2 = zc(i+1);

        seg   = x_win(z1:z2);
        seg_t = t_win(z1:z2);

        % must contain a negative deflection
        [minVal,minIdx] = min(seg);
        if minVal >= 0
            continue
        end

        % duration 0.2–1.0 s
        dur = seg_t(end) - seg_t(1);
        if dur < 0.2 || dur > 1.0
            continue
        end

        SWpeak = [SWpeak; seg_t(minIdx)];
        SWseg  = [SWseg; seg_t(1) seg_t(end)];
    end

    %% SPINDLE ASSOCIATION 

    assoc_idx = [];
    for k = 1:length(SWpeak)
        latency = SWpeak(k) - t_on;
        if latency >= -0.8 && latency <= 0.3
            assoc_idx(end+1) = k;
        end
    end

    % KEEP ONLY spindle-associated SW
    SWseg  = SWseg(assoc_idx,:);
    SWpeak = SWpeak(assoc_idx);

    %% PLOT 
    subplot(nrows, ncols, s);

    ymin = -150;
    ymax =  150;

    % spindle window (light red)
    patch([t_on t_off t_off t_on], [ymin ymin ymax ymax], ...
          [1 0.7 0.7], 'FaceAlpha',0.35, 'EdgeColor','none');
    hold on;

    % spindle-associated slow waves ONLY (blue)
    for k = 1:size(SWseg,1)
        xs = SWseg(k,1);
        xe = SWseg(k,2);
        patch([xs xe xe xs], [ymin ymin ymax ymax], ...
              [0.7 0.85 1.0], 'FaceAlpha',0.55, 'EdgeColor','none');
    end

    % broadband trace
    plot(t_win, x_win, 'k', 'LineWidth', 1);

    % formatting
    title(sprintf('Spindle %d', s), 'FontName','Helvetica','FontSize',10);
    xlabel('Time (s)', 'FontName','Helvetica','FontSize',9);
    ylabel('\muV', 'FontName','Helvetica','FontSize',9);

    set(gca, 'FontName','Helvetica', 'FontSize',9, ...
             'LineWidth',1, 'Box','off');

    xlim([t_win(1) t_win(end)]);
    ylim([-150 150]);   % FIXED Y-LIMITS
end

%% CROSS-VERIFY SPINDLE-ASSOCIATED SLOW WAVES ACROSS ALL CHANNELS

fs = 1000;
numCh = size(amplifier_data_filtered,1);

% output: spindle x channel (1 = associated SW exists)
SW_assoc_matrix = zeros(size(spindles,1), numCh);

for ch = 1:numCh

    x_full = amplifier_data_filtered(ch,:);
    t_full = t_amplifier;

    for s = 1:size(spindles,1)

        t_on  = spindles(s,1);
        t_off = spindles(s,2);

        % window (±0.5s)
        idx   = t_full >= (t_on - 0.5) & t_full <= (t_off + 0.5);
        t_win = t_full(idx);
        x_win = x_full(idx);

        if numel(t_win) < 100
            continue
        end

        %% DETECT SLOW WAVES (0.2–1.0s NEG HALF-WAVES) 

        zc = find(diff(sign(x_win)) ~= 0);
        SWpeak = [];

        for i = 1:length(zc)-1

            z1 = zc(i);
            z2 = zc(i+1);

            seg   = x_win(z1:z2);
            seg_t = t_win(z1:z2);

            [minVal, minIdx] = min(seg);
            if minVal >= 0
                continue
            end

            dur = seg_t(end) - seg_t(1);
            if dur < 0.2 || dur > 1.0
                continue
            end

            SWpeak = [SWpeak; seg_t(minIdx)];
        end

        %% SPINDLE ASSOCIATION

        for k = 1:length(SWpeak)
            latency = SWpeak(k) - t_on;
            if latency >= -0.8 && latency <= +0.3
                SW_assoc_matrix(s,ch) = 1;
                break
            end
        end

    end
end

disp('Matrix of spindle-associated slow waves (spindles x channels):');
disp(SW_assoc_matrix);

% Spindle 10: Highlighted EEG traces for Channels 1–4

fs = 1000;
t_on  = spindles(10,1);
t_off = spindles(10,2);

figure('Color','w');

for ch = 1:4

    %%  Window extraction 
    x_full = amplifier_data_filtered(ch,:);
    t_full = t_amplifier;

    idx = t_full >= (t_on - 0.5) & t_full <= (t_off + 0.5);
    t_win = t_full(idx);
    x_win = x_full(idx);

    %% Detect slow waves 
    zc = find(diff(sign(x_win)) ~= 0);
    SWpeak = [];
    SWseg  = [];

    for i = 1:length(zc)-1
        z1 = zc(i);
        z2 = zc(i+1);

        seg   = x_win(z1:z2);
        seg_t = t_win(z1:z2);

        [minVal,minIdx] = min(seg);
        if minVal >= 0, continue; end

        dur = seg_t(end) - seg_t(1);
        if dur < 0.2 || dur > 1.0, continue; end

        SWpeak(end+1,1) = seg_t(minIdx);
        SWseg(end+1,:)  = [seg_t(1) seg_t(end)];
    end

    % Only spindle-associated SW
    assoc_idx = find( (SWpeak - t_on) >= -0.8 & (SWpeak - t_on) <= 0.3 );
    SWseg = SWseg(assoc_idx,:);

    % PLOT 
    subplot(4,1,ch);

    % entire trace in black
    plot(t_win, x_win, 'k', 'LineWidth', 1); hold on;

    % slow-wave segments in blue
    for k = 1:size(SWseg,1)
        idx_sw = t_win >= SWseg(k,1) & t_win <= SWseg(k,2);
        plot(t_win(idx_sw), x_win(idx_sw), 'b', 'LineWidth', 2);
    end

    % spindle in red (draw last → overwrites blue if overlapping)
    idx_sp = t_win >= t_on & t_win <= t_off;
    plot(t_win(idx_sp), x_win(idx_sp), 'r', 'LineWidth', 2);

    % formatting
    title(sprintf('Channel %d – Spindle 10', ch), ...
        'FontName','Helvetica','FontSize',11);
    ylabel('\muV','FontName','Helvetica','FontSize',10);

    set(gca,'FontName','Helvetica','FontSize',10, ...
            'LineWidth',1,'Box','off');

    ylim([-220 200]);           % FIXED Y-LIMITS
    xlim([t_win(1) t_win(end)]);

    if ch == 4
        xlabel('Time (s)','FontName','Helvetica','FontSize',10);
    end
end

%EXPORT AS EPS 
exportgraphics(gcf,'Spindle10_Ch1-4_highlighted.eps','ContentType','vector');

% Spindle 3: Highlighted EEG traces for Channels 1–4 

fs = 1000;
t_on  = spindles(3,1);
t_off = spindles(3,2);

figure('Color','w');

for ch = 1:4

    %% Window extraction
    x_full = amplifier_data_filtered(ch,:);
    t_full = t_amplifier;

    idx = t_full >= (t_on - 0.5) & t_full <= (t_off + 0.5);
    t_win = t_full(idx);
    x_win = x_full(idx);

    %%  Detect slow waves (negative half-wave, 0.2–1.0 s)
    zc = find(diff(sign(x_win)) ~= 0);
    SWpeak = [];
    SWseg  = [];

    for i = 1:length(zc)-1
        z1 = zc(i);
        z2 = zc(i+1);

        seg   = x_win(z1:z2);
        seg_t = t_win(z1:z2);

        [minVal,minIdx] = min(seg);
        if minVal >= 0, continue; end

        dur = seg_t(end) - seg_t(1);
        if dur < 0.2 || dur > 1.0, continue; end

        SWpeak(end+1,1) = seg_t(minIdx);
        SWseg(end+1,:)  = [seg_t(1) seg_t(end)];
    end

    %%  Keep only spindle-associated slow waves 
    assoc_idx = find( (SWpeak - t_on) >= -0.8 & (SWpeak - t_on) <= 0.3 );
    SWseg = SWseg(assoc_idx,:);

    %% Plot 
    subplot(4,1,ch);

    % Full trace in black
    plot(t_win, x_win, 'k', 'LineWidth', 1); hold on;

    % Slow waves in blue
    for k = 1:size(SWseg,1)
        idx_sw = t_win >= SWseg(k,1) & t_win <= SWseg(k,2);
        plot(t_win(idx_sw), x_win(idx_sw), 'b', 'LineWidth', 2);
    end

    % Spindle in red (draw last = overwrites blue)
    idx_sp = t_win >= t_on & t_win <= t_off;
    plot(t_win(idx_sp), x_win(idx_sp), 'r', 'LineWidth', 2);

    % Formatting
    title(sprintf('Channel %d – Spindle 3', ch), ...
        'FontName','Helvetica','FontSize',11);
    ylabel('\muV','FontName','Helvetica','FontSize',10);

    set(gca,'FontName','Helvetica','FontSize',10, ...
            'LineWidth',1,'Box','off');

    ylim([-220 200]);        % <—— NEW LIMITS
    xlim([t_win(1) t_win(end)]);

    if ch == 4
        xlabel('Time (s)','FontName','Helvetica','FontSize',10);
    end
end

%% EXPORT AS EPS 
exportgraphics(gcf,'Spindle3_Ch1-4_highlighted.eps','ContentType','vector');

%% Multitaper spectrogram for Channel 1

ch = 1;
x = double(amplifier_data_filtered(ch,:));
x = x - mean(x,'omitnan');

FREQ_LIMS = [1 20];
WIN_SEC = 2.0;
OL_FRAC = 0.9;
NW = 3.5;
NFFT_FACTOR = 4;

winlen = max(64, round(WIN_SEC * fs));
nover  = max(0, round(OL_FRAC * winlen));
step   = max(1, winlen - nover);
nfft   = 2^nextpow2(max(winlen, round(NFFT_FACTOR * fs)));

nFrames = 1 + floor((numel(x) - winlen) / step);

seg = x(1:winlen);
[Pxx, F] = pmtm(seg, NW, nfft, fs);
fsel = F >= FREQ_LIMS(1) & F <= FREQ_LIMS(2);
Fplt = F(fsel);

Smt  = zeros(nnz(fsel), nFrames);
Tmid = zeros(1, nFrames);

idx0 = 1;
for kf = 1:nFrames
    seg = x(idx0:idx0 + winlen - 1);
    [Pxx,~] = pmtm(seg, NW, nfft, fs);
    Smt(:,kf) = Pxx(fsel);

    mid = idx0 + floor((winlen - 1)/2);
    Tmid(kf) = t_amplifier(mid);

    idx0 = idx0 + step;
end

fig = figure('Color','w','Units','inches','Position',[1 1 7.73 3]);
ax = axes(fig);

imagesc(ax, Tmid, Fplt, 10*log10(Smt + eps));
axis(ax,'xy');
colormap(ax,turbo);
set(ax,'FontName','Helvetica','FontSize',9,'LineWidth',1,'Box','off');
xlabel(ax,'Time (s)');
ylabel(ax,'Frequency (Hz)');
title(ax,'Channel 1 — Multitaper spectrogram (1–20 Hz)');

caxis(ax,[0 25]);
colorbar(ax);


% Zoom around Spindle 10 

t_on  = spindles(10,1);
t_off = spindles(10,2);

lo = t_on - 2;
hi = t_off + 2;

xlim(ax,[lo hi]);
ylim(ax,[1 20]);

exportgraphics(fig,'Ch1_MT_zoom_Spindle10_1-20Hz.eps','ContentType','vector');
