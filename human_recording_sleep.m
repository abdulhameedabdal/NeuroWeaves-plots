%% Multitaper spectrogram (per channel) — customizable file
% Change these:
FILE        = 'full_sleeping';      % tag or exact filename (e.g., 'sleep_2','full_sleeping','sleep_6.csv')
CHS         = {'c4'};               % channels to plot (any of: 'C4','O2','Fp2')
FREQ_LIMS   = [0.5 20];             % Hz band to display
CLIM_DB     = [0 20];               % e.g., [0 60]; [] = auto
WIN_SEC     = 2.0;                  % window length (s)
OL_FRAC     = 0.9;                  % 0–0.95 overlap fraction
NW          = 3.5;                  % time-halfbandwidth (typical 2.5–4)
NFFT_FACTOR = 4;                    % nfft ≈ NFFT_FACTOR*fs (rounded to next pow2)

% --- Load table T (uses S if present, else baseDir or picker) ---
clear T
tag = lower(erase(FILE,'.csv'));
if exist('S','var')==1 && ~isempty(S)
    k = find(contains(lower(string({S.name})), tag), 1, 'first');
    if ~isempty(k)
        T   = S(k).T;
        src = "S:"+S(k).name;
    end
end
if ~exist('T','var') || isempty(T)
    if exist('baseDir','var')==1 && isfolder(baseDir)
        d = dir(fullfile(baseDir, ['*' tag '*']));
        if ~isempty(d)
            fpath = fullfile(d(1).folder, d(1).name);
        else
            [fn,fp] = uigetfile(fullfile(baseDir,'*.csv'),'Pick CSV'); assert(~isequal(fn,0));
            fpath = fullfile(fp,fn);
        end
    else
        [fn,fp] = uigetfile('*.csv','Pick CSV'); assert(~isequal(fn,0));
        fpath = fullfile(fp,fn);
    end
    T   = readtable(fpath,'TextType','string');
    src = "file:"+string(fpath);
end
disp("Using: " + src);

% Timebase & fs 
vars = T.Properties.VariableNames; vlow = lower(vars);
tcol = find(strcmpi(vars,'time') | contains(vlow,'time'),1,'first');
assert(~isempty(tcol),'No time column found.');
t    = T{:, tcol};  t = t(:);
fs   = 1/median(diff(t));
trel = t - t(1);

%  MT spectrogram per channel in CHS 
for c = 1:numel(CHS)
    ci = find(strcmpi(vars, CHS{c}), 1, 'first');
    if isempty(ci)
        warning('Channel %s not found. Skipping.', CHS{c});
        continue;
    end
    x = double(T{:,ci});
    x = x - mean(x,'omitnan');

    winlen = max(64, round(WIN_SEC*fs));
    nover  = max(0, round(OL_FRAC*winlen));
    step   = max(1, winlen - nover);
    nfft   = 2^nextpow2(max(winlen, round(NFFT_FACTOR*fs)));

    nFrames = 1 + floor((numel(x)-winlen)/step);
    if nFrames < 1
        warning('Signal too short for settings.');
        continue;
    end

    % Precompute size using first frame
    seg      = x(1:winlen);
    [Pxx,F]  = pmtm(seg, NW, nfft, fs, 'DropLastTaper', true);
    fsel     = F>=FREQ_LIMS(1) & F<=min(FREQ_LIMS(2),(fs/2)-0.1);
    Fplt     = F(fsel);
    Smt      = zeros(nnz(fsel), nFrames);
    Tmid     = zeros(1, nFrames);

    % Slide & estimate multitaper PSD
    idx0 = 1;
    for kf = 1:nFrames
        seg = x(idx0:idx0+winlen-1);
        [Pxx,~]     = pmtm(seg, NW, nfft, fs, 'DropLastTaper', true);
        Smt(:,kf)   = Pxx(fsel);
        mid         = idx0 + floor((winlen-1)/2);
        Tmid(kf)    = trel(mid);
        idx0        = idx0 + step;
    end

    % Plot
    figure('Color','w','Units','inches','Position',[1 1 12 5]);
    imagesc(Tmid, Fplt, 10*log10(Smt + eps)); axis xy; grid on; box on;
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    title(sprintf('%s — Multitaper spectrogram (NW=%.1f, Win=%.1fs, Overlap=%.0f%%)', ...
          upper(CHS{c}), NW, WIN_SEC, 100*OL_FRAC), 'Interpreter','none');
    colormap(turbo);
    cb = colorbar; cb.Label.String = 'Power (dB)';
    if ~isempty(CLIM_DB)
        caxis(CLIM_DB);
    end

    % EPS export (vector)
    if exist('baseDir','var')==1 && isfolder(baseDir)
        outDir = baseDir;
    else
        outDir = pwd;
    end
    outName = sprintf('%s_mt_spec_%s.eps', tag, lower(CHS{c}));
    outPath = fullfile(outDir, outName);
    exportgraphics(gcf, outPath, 'ContentType','vector');
    fprintf('Saved EPS: %s\n', outPath);
end

%%  Zoom C4 spectrogram to X–X s (freq 0.5–20 Hz) 
lo=445; hi=452; fmin=0.5; fmax=20;
ax=gca; xlim(ax,[lo hi]); ylim(ax,[fmin fmax]); set(ax,'YDir','normal'); grid(ax,'on'); box(ax,'on');
exportgraphics(gcf, fullfile(pwd,'C4_MT_zoom_90_120s_0p5_20Hz.eps'), 'ContentType','vector');

%% Plot raw traces for Fp2, C4, O2 with fixed symmetric Y limits 
% --- settings ---
TIME_LIMS    = [445 452];    % [start end] in seconds
DETREND_MEAN = true;         % subtract mean
FIX_MAG_UV   = 110;          % <<< change this magnitude (e.g., 80, 120). Y-limits = ±FIX_MAG_UV

% load CSV (uses baseDir if available, else picker) 
clear T
if exist('baseDir','var')==1 && isfolder(baseDir)
    [fn,fp] = uigetfile(fullfile(baseDir,'*.csv'),'Pick the CSV'); assert(~isequal(fn,0));
else
    [fn,fp] = uigetfile('*.csv','Pick the CSV'); assert(~isequal(fn,0));
end
T = readtable(fullfile(fp,fn),'TextType','string');
vars = T.Properties.VariableNames; vlow = lower(vars);

%  time (auto-fix ms->s) 
ti = find(strcmpi(vars,'time') | contains(vlow,'time'),1,'first');
t  = T{:,ti}; t = t(:);
fs = 1/median(diff(t),'omitnan'); if fs>2000, t = t/1000; end
trel = t - t(1);

% pull channels in this order 
chs = {'Fp2','C4','O2'};
idx = cellfun(@(c) find(strcmpi(vars,c),1,'first'), chs);
assert(all(~isnan(idx)),'Missing one of Fp2/C4/O2.');
X = double(T{:,idx}); if DETREND_MEAN, X = X - mean(X,1,'omitnan'); end

% select window and plot stacked (fixed ±FIX_MAG_UV for all)
sel = trel >= TIME_LIMS(1) & trel <= TIME_LIMS(2);
figure('Color','w','Units','inches','Position',[1 1 12 6]);
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');
for k = 1:3
    nexttile; plot(trel(sel), X(sel,k), 'LineWidth',1.1);
    ylabel([chs{k} ' (\muV)']); grid on; ylim([-FIX_MAG_UV FIX_MAG_UV]);
    if k<3, set(gca,'XTickLabel',[]); else, xlabel('Time (s)'); end
end
title(sprintf('%s — Fp2, C4, O2 (%.1f–%.1f s) | Y = ±%g µV', fn, TIME_LIMS(1), TIME_LIMS(2), FIX_MAG_UV), 'Interpreter','none');

% export EPS
set(gcf,'Renderer','painters');
exportgraphics(gcf, fullfile(pwd, sprintf('traces_Fp2_C4_O2_%gto%gs_pm%guV.eps', TIME_LIMS(1), TIME_LIMS(2), FIX_MAG_UV)), 'ContentType','vector');

%% check 446-454 for probable spindle 515-535
%% Spindle detector (paper-style) + raw traces in chosen window
% Params (edit)
TIME_LIMS   = [445 452];              % seconds
CHS         = {'Fp2','C4','O2'};
BAND_HZ     = [11 16];                % 11–16 Hz
THRESH_UV   = 3;                      % envelope ≥ 3 µV
DUR_RANGE_S = [0.5 2];                % 0.5–2 s
DETREND_MEAN = true;

% Load once (same pattern as your traces block)
clear T
if exist('baseDir','var')==1 && isfolder(baseDir)
    [fn,fp]=uigetfile(fullfile(baseDir,'*.csv'),'Pick CSV'); assert(~isequal(fn,0));
else
    [fn,fp]=uigetfile('*.csv','Pick CSV'); assert(~isequal(fn,0));
end
T = readtable(fullfile(fp,fn),'TextType','string');
vars = T.Properties.VariableNames; vlow = lower(vars);
ti = find(strcmpi(vars,'time') | contains(vlow,'time'),1,'first');
t  = T{:,ti}; t = t(:);
fs = 1/median(diff(t),'omitnan'); if fs>2000, t = t/1000; end
trel = t - t(1);

% Pull channels
idx = cellfun(@(c) find(strcmpi(vars,c),1,'first'), CHS);
assert(all(~isnan(idx)),'Missing Fp2/C4/O2.');
X = double(T{:,idx}); if DETREND_MEAN, X = X - mean(X,1,'omitnan'); end
sel = trel>=TIME_LIMS(1) & trel<=TIME_LIMS(2);

% Detect spindles (per paper: bandpass 11–16 Hz, Hilbert envelope, ≥3 µV for 0.5–2 s)
wins = cell(1,numel(CHS));   % { [start end; ...] } in seconds
for k = 1:numel(CHS)
    xf = bandpass(X(:,k), BAND_HZ, fs,'Steepness',0.85,'StopbandAttenuation',60);
    env = abs(hilbert(xf));                          % µV
    above = env >= THRESH_UV;
    d = diff([false; above; false]);
    starts = find(d==1); stops = find(d==-1)-1;
    durs = (stops - starts + 1)/fs;
    keep = durs>=DUR_RANGE_S(1) & durs<=DUR_RANGE_S(2);
    wins{k} = [trel(starts(keep)) trel(stops(keep))];
end

% Plot raw traces in window with spindle highlights
figure('Color','w','Units','inches','Position',[1 1 12 6]);
tiledlayout(numel(CHS),1,'Padding','compact','TileSpacing','compact');
for k = 1:numel(CHS)
    tk = trel(sel); xk = X(sel,k);
    mag = 1.05*max(abs(xk));
    nexttile; plot(tk, xk,'LineWidth',1.1); grid on; ylim([-mag mag]);
    ylabel([CHS{k} ' (\muV)']);
    % highlight detections that intersect the window
    W = wins{k};
    for r = 1:size(W,1)
        a = max(W(r,1), TIME_LIMS(1)); b = min(W(r,2), TIME_LIMS(2));
        if a<b, patch([a b b a],[ -mag -mag mag mag],[1 0.6 0.6], ...
                      'FaceAlpha',0.25,'EdgeColor','none'); end
    end
    uistack(findobj(gca,'Type','line'),'top'); % keep trace above patches
    if k<numel(CHS), set(gca,'XTickLabel',[]); else, xlabel('Time (s)'); end
end
title(sprintf('%s — Spindle paper-method (11–16 Hz, ≥%g µV, %g–%g s) | %.1f–%.1f s', ...
      fn, THRESH_UV, DUR_RANGE_S(1), DUR_RANGE_S(2), TIME_LIMS(1), TIME_LIMS(2)), 'Interpreter','none');

%% Spindles - bandpass(11–16 Hz) + Hilbert envelope (per channel)
% Settings
TIME_LIMS   = [446 454];          % seconds
CHS         = {'Fp2','C4','O2'};
BAND_HZ     = [11 16];
THRESH_UV   = 3;                  % envelope ≥ 3 µV
DUR_RANGE_S = [0.5 2];            % 0.5–2 s
DETREND_MEAN = true;

% Load CSV (same pattern as before)
clear T
if exist('baseDir','var')==1 && isfolder(baseDir)
    [fn,fp]=uigetfile(fullfile(baseDir,'*.csv'),'Pick CSV'); assert(~isequal(fn,0));
else
    [fn,fp]=uigetfile('*.csv','Pick CSV'); assert(~isequal(fn,0));
end
T = readtable(fullfile(fp,fn),'TextType','string');
vars = T.Properties.VariableNames; vlow = lower(vars);
ti = find(strcmpi(vars,'time') | contains(vlow,'time'),1,'first');
t  = T{:,ti}; t = t(:); fs = 1/median(diff(t),'omitnan'); if fs>2000, t = t/1000; end
trel = t - t(1);

% Pull channels
idx = cellfun(@(c) find(strcmpi(vars,c),1,'first'), CHS);
assert(all(~isnan(idx)),'Missing Fp2/C4/O2.');
X = double(T{:,idx}); if DETREND_MEAN, X = X - mean(X,1,'omitnan'); end
sel = trel>=TIME_LIMS(1) & trel<=TIME_LIMS(2);
tk  = trel(sel);

% Prep figure (two traces per channel: filtered + envelope)
figure('Color','w','Units','inches','Position',[1 1 12 8]);
tiledlayout(numel(CHS),1,'Padding','compact','TileSpacing','compact');

for k = 1:numel(CHS)
    x = X(:,k);

    % Bandpass 11–16 Hz (zero-phase), Hilbert envelope
    xf = bandpass(x, BAND_HZ, fs,'Steepness',0.85,'StopbandAttenuation',60);
    env = abs(hilbert(xf));              % µV, same scale as x

    % Paper detection (no tweaks): envelope ≥3 µV for 0.5–2 s
    above = env >= THRESH_UV;
    d = diff([false; above; false]);
    starts = find(d==1); stops = find(d==-1)-1;
    durs = (stops - starts + 1)/fs;
    keep = durs>=DUR_RANGE_S(1) & durs<=DUR_RANGE_S(2);
    winS = trel(starts(keep));  winE = trel(stops(keep));

    % Windowed plotting
    xf_w  = xf(sel);
    env_w = env(sel);
    mag   = 1.1*max([abs(xf_w); env_w]);   % symmetric magnitude for visibility

    nexttile; hold on;
    % shade detected spindle windows (that intersect TIME_LIMS)
    for r=1:numel(winS)
        a = max(winS(r), TIME_LIMS(1)); b = min(winE(r), TIME_LIMS(2));
        if a<b, patch([a b b a],[-mag -mag mag mag],[1 .8 .8], ...
                      'EdgeColor','none','FaceAlpha',0.35); end
    end
    % plot filtered 11–16 Hz and its envelope
    plot(tk, xf_w, 'LineWidth',1.0);                 % bandpassed signal (µV)
    plot(tk, env_w, 'LineWidth',1.5);                % Hilbert envelope (µV)
    yline(THRESH_UV,'k--','LineWidth',1);            % 3 µV threshold

    grid on; ylim([-mag mag]);
    ylabel([CHS{k} ' (\muV)']);
    legend({'detected window','11–16 Hz','envelope','3 \muV'},'Location','northeast');
    if k<numel(CHS), set(gca,'XTickLabel',[]); else, xlabel('Time (s)'); end
    title(sprintf('%s — 11–16 Hz bandpass + Hilbert envelope (%.1f–%.1f s)', ...
          CHS{k}, TIME_LIMS(1), TIME_LIMS(2)));
end



%% Slow-wave detector + per-channel plot in chosen window
% settings 
TIME_LIMS     = [445 452];          % [start end] in seconds
CHS           = {'Fp2','C4','O2'};  % channels to analyze
BP_HZ         = [0.2 4];            % bandpass for slow waves
NEG_PK_UV     = -60;                % negative peak threshold (µV)
PP_MIN_UV     = 75;                 % peak-to-peak minimum (µV)
DUR_ZC_RANGE  = [0.2 1.0];          % neg-ZC → pos-ZC duration (s)
DETREND_MEAN  = true;

% --- load CSV (uses baseDir if available, else picker) ---
clear T
if exist('baseDir','var')==1 && isfolder(baseDir)
    [fn,fp]=uigetfile(fullfile(baseDir,'*.csv'),'Pick CSV'); assert(~isequal(fn,0));
else
    [fn,fp]=uigetfile('*.csv','Pick CSV'); assert(~isequal(fn,0));
end
T = readtable(fullfile(fp,fn),'TextType','string');
vars = T.Properties.VariableNames; vlow = lower(vars);

% time (auto-fix ms->s) 
ti = find(strcmpi(vars,'time') | contains(vlow,'time'),1,'first');
t  = T{:,ti}; t = t(:);
fs = 1/median(diff(t),'omitnan'); if fs>2000, t = t/1000; end
trel = t - t(1);

% pull channels 
idx = cellfun(@(c) find(strcmpi(vars,c),1,'first'), CHS);
assert(all(~isnan(idx)),'Missing one of Fp2/C4/O2.');
X = double(T{:,idx}); if DETREND_MEAN, X = X - mean(X,1,'omitnan'); end

% detect slow waves per channel 
SW = cell(1,numel(CHS));  % each: [t_start t_end] seconds
for k = 1:numel(CHS)
    x  = X(:,k);
    xf = bandpass(x, BP_HZ, fs,'Steepness',0.85,'StopbandAttenuation',60);   % 0-phase
    % zero-crossings
    negZ = find(xf(1:end-1)>0 & xf(2:end)<=0);   % positive→negative
    posZ = find(xf(1:end-1)<0 & xf(2:end)>=0);   % negative→positive
    ev   = zeros(0,2);
    pi = 1;
    for i = 1:numel(negZ)
        nz = negZ(i);
        while pi<=numel(posZ) && posZ(pi) <= nz, pi = pi+1; end
        if pi>numel(posZ), break; end
        pz = posZ(pi);
        dur = (pz - nz)/fs;
        if dur < DUR_ZC_RANGE(1) || dur > DUR_ZC_RANGE(2), continue; end
        % negative peak between nz..pz
        [negPk, iMinRel] = min(xf(nz:pz));
        iMin = nz + iMinRel - 1;
        if negPk > NEG_PK_UV, continue; end
        % next positive peak within 1 s after negative peak
        iEnd = min(length(xf), iMin + round(1*fs));
        [posPk, ~] = max(xf(iMin:iEnd));
        if (posPk - negPk) < PP_MIN_UV, continue; end
        ev(end+1,:) = [trel(nz) trel(pz)]; %#ok<AGROW> % mark nz→pz as the slow-wave interval
    end
    SW{k} = ev;
end

% --- plot bandpassed traces with detected slow waves highlighted (TIME_LIMS) ---
sel = trel >= TIME_LIMS(1) & trel <= TIME_LIMS(2);
tk  = trel(sel);
figure('Color','w','Units','inches','Position',[1 1 12 6]);
tiledlayout(numel(CHS),1,'Padding','compact','TileSpacing','compact');
for k = 1:numel(CHS)
    xf = bandpass(X(:,k), BP_HZ, fs,'Steepness',0.85,'StopbandAttenuation',60);
    xw = xf(sel);
    mag = 1.1*max(abs(xw)+eps);
    nexttile; hold on;
    % shade detections intersecting window
    ev = SW{k};
    for r = 1:size(ev,1)
        a = max(ev(r,1), TIME_LIMS(1)); b = min(ev(r,2), TIME_LIMS(2));
        if a<b, patch([a b b a],[-mag -mag mag mag],[0.9 0.95 1], 'EdgeColor','none','FaceAlpha',0.6); end
    end
    plot(tk, xw, 'LineWidth',1.2); grid on; ylim([-mag mag]);
    ylabel([CHS{k} ' (\muV)']);
    if k<numel(CHS), set(gca,'XTickLabel',[]); else, xlabel('Time (s)'); end
end
title(sprintf('%s — Slow waves (%.1f–%.1f s) | BP %.1f–%.1f Hz, neg≤%g µV, P2P≥%g µV, ZC %.1f–%.1f s', ...
      fn, TIME_LIMS(1), TIME_LIMS(2), BP_HZ(1), BP_HZ(2), NEG_PK_UV, PP_MIN_UV, DUR_ZC_RANGE(1), DUR_ZC_RANGE(2)), 'Interpreter','none');



%% Alpha power bar chart (Fp2, C4, O2) over a chosen window
%  settings
TIME_LIMS   = [98 105];      % [start end] in seconds  <-- edit
ALPHA_BAND  = [8 13];         % Hz
DETREND_MEAN= true;           % subtract mean before PSD
WIN_SEC     = 4;              % Welch window (s)
OL_FRAC     = 0.5;            % 50% overlap
NFFT_FACTOR = 8;              % nfft ≈ NFFT_FACTOR*fs
SHOW_DB     = true;           % plot in dB if true
EXPORT_EPS  = false;          % set true to save EPS

% --- load CSV (uses baseDir if available, else picker) ---
clear T
if exist('baseDir','var')==1 && isfolder(baseDir)
    [fn,fp] = uigetfile(fullfile(baseDir,'*.csv'),'Pick the CSV'); assert(~isequal(fn,0));
else
    [fn,fp] = uigetfile('*.csv','Pick the CSV'); assert(~isequal(fn,0));
end
T = readtable(fullfile(fp,fn),'TextType','string');
vars = T.Properties.VariableNames; vlow = lower(vars);

%  time & fs (auto-fix ms->s)
ti = find(strcmpi(vars,'time') | contains(vlow,'time'),1,'first');
t  = T{:,ti}; t = t(:);
fs = 1/median(diff(t),'omitnan'); if fs>2000, t = t/1000; fs = 1/median(diff(t),'omitnan'); end
trel = t - t(1);

%  channels 
chs = {'Fp2','C4','O2'};
cidx = zeros(1,3);
for k=1:3, cidx(k)=find(strcmpi(vars,chs{k}),1,'first'); assert(~isempty(cidx(k)),['Missing ' chs{k}]); end

% select window 
sel = trel>=TIME_LIMS(1) & trel<=TIME_LIMS(2);
assert(nnz(sel)>max(64,round(WIN_SEC*fs)), 'Window too short for Welch.');

% Welch params 
winlen = max(64, round(WIN_SEC*fs));
nover  = round(OL_FRAC*winlen);
nfft   = 2^nextpow2(max(winlen, round(NFFT_FACTOR*fs)));

% alpha power per channel (µV^2) 
Palpha = zeros(1,3);
for k=1:3
    x = double(T{sel,cidx(k)});
    if DETREND_MEAN, x = x - mean(x,'omitnan'); end
    [Pxx,F] = pwelch(x, hamming(winlen), nover, nfft, fs, 'psd');  % µV^2/Hz
    m = F>=ALPHA_BAND(1) & F<=ALPHA_BAND(2);
    Palpha(k) = trapz(F(m), Pxx(m));                               % µV^2
end
vals = Palpha;
ylab = 'Alpha power (\muV^2)';
if SHOW_DB
    vals = 10*log10(Palpha + eps);
    ylab = 'Alpha power (dB \muV^2)';
end

% bar plot
figure('Color','w','Units','inches','Position',[1 1 6 4]);
bh = bar(categorical(chs), vals); grid on; box on;
ylabel(ylab); title(sprintf('Alpha %d–%d Hz | %s | %.1f–%.1f s', ...
       ALPHA_BAND(1), ALPHA_BAND(2), fn, TIME_LIMS(1), TIME_LIMS(2)), 'Interpreter','none');
% value labels
xt = get(gca,'XTick'); yl = ylim; off = 0.02*(yl(2)-yl(1));
for k=1:numel(vals), text(xt(k), vals(k)+off, sprintf('%.2f', vals(k)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom'); end

% --- optional export ---
    set(gcf,'Renderer','painters');
    exportgraphics(gcf, fullfile(fp, sprintf('alpha_bar_%s_%.0f-%.0fs.eps', erase(fn,'.csv'), TIME_LIMS)), 'ContentType','vector');

%% K-complex zoom: raw (unfiltered) trace with inset
% edit these 
CH         = 'Fp2';           % 'Fp2' | 'C4' | 'O2'
TIME_LIMS  = [364 371];       % main window [s]
ZOOM_LIMS  = [369.8 370.5]; % inset window [s] around the K wave
Y_MAG      = 130;             % fixed y-range: [-Y_MAG, +Y_MAG] µV
DETREND_MEAN = false;         % true = subtract mean for display only

% load CSV (picker) 
clear T
[fn,fp] = uigetfile('*.csv','Pick the CSV'); assert(~isequal(fn,0));
T    = readtable(fullfile(fp,fn),'TextType','string');
vars = T.Properties.VariableNames; vlow = lower(vars);

% time (auto-fix ms->s)
ti = find(strcmpi(vars,'time') | contains(vlow,'time'),1,'first');
t  = T{:,ti}; t = t(:);
fs = 1/median(diff(t),'omitnan'); if fs>2000, t = t/1000; end
trel = t - t(1);

% channel 
ci = find(strcmpi(vars,CH),1,'first'); assert(~isempty(ci), 'Channel not found');
x  = double(T{:,ci}); if DETREND_MEAN, x = x - mean(x,'omitnan'); end

% selections 
selMain = trel >= TIME_LIMS(1) & trel <= TIME_LIMS(2);
selZoom = trel >= ZOOM_LIMS(1) & trel <= ZOOM_LIMS(2);

% plot main + inset
figure('Color','w','Units','inches','Position',[1 1 11 5]); 
axes('Position',[0.08 0.15 0.78 0.78]); % main
plot(trel(selMain), x(selMain), 'LineWidth',1.2); grid on; box on;
ylim([-Y_MAG Y_MAG]); xlim(TIME_LIMS);
xlabel('Time (s)'); ylabel([CH ' (\muV)']);
title(sprintf('%s — %s raw trace  (%.1f–%.1f s)', fn, CH, TIME_LIMS(1), TIME_LIMS(2)), 'Interpreter','none');

% rectangle showing inset region
hold on;
yr = ylim; 
patch([ZOOM_LIMS(1) ZOOM_LIMS(2) ZOOM_LIMS(2) ZOOM_LIMS(1)], ...
      [yr(1) yr(1) yr(2) yr(2)], [0.95 0.95 1], ...
      'FaceAlpha',0.15, 'EdgeColor',[0.2 0.4 1], 'LineWidth',1.2);

% inset axes
axInset = axes('Position',[0.70 0.58 0.25 0.30]); %#ok<LAXES>
plot(trel(selZoom), x(selZoom), 'LineWidth',1.3); grid on; box on;
xlim(ZOOM_LIMS); ylim([-Y_MAG Y_MAG]);
title('Zoom', 'FontWeight','normal'); set(axInset, 'XAxisLocation','bottom');

% --- export (EPS & PNG) ---
set(gcf,'Renderer','painters');
outBase = sprintf('Kzoom_%s_%0.1fto%0.1fs', CH, TIME_LIMS(1), TIME_LIMS(2));
exportgraphics(gcf, fullfile(fp, [outBase '.eps']), 'ContentType','vector');
exportgraphics(gcf, fullfile(fp, [outBase '.png']), 'Resolution',300);

%% Spindle segment (C4 shown bandpassed + envelope) and print spindle power (all three)
%  settings
TIME_LIMS    = [449.8 451.3];        % [start end] seconds (edit)
CHS          = {'Fp2','C4','O2'};    % channels
BAND_HZ      = [11 16];              % spindle band
ENV_THR_UV   = 3;                    % envelope threshold (µV)
DUR_RANGE_S  = [0.5 2.0];            % contiguous envelope duration to keep (s)
DETREND_MEAN = true;

% load CSV (uses baseDir if available, else picker) 
clear T
if exist('baseDir','var')==1 && isfolder(baseDir)
    [fn,fp]=uigetfile(fullfile(baseDir,'*.csv'),'Pick CSV'); assert(~isequal(fn,0));
else
    [fn,fp]=uigetfile('*.csv','Pick CSV'); assert(~isequal(fn,0));
end
T = readtable(fullfile(fp,fn),'TextType','string');
vars = T.Properties.VariableNames; vlow = lower(vars);

% time/fs (auto-fix ms->s)
ti = find(strcmpi(vars,'time') | contains(vlow,'time'),1,'first');
t  = T{:,ti}; t = t(:);
fs = 1/median(diff(t),'omitnan');
if fs>2000, t = t/1000; fs = 1/median(diff(t),'omitnan'); end   % timestamps were ms
trel = t - t(1);

% pull channels 
idx = cellfun(@(c) find(strcmpi(vars,c),1,'first'), CHS);
assert(all(~isnan(idx)),'Missing one of Fp2/C4/O2.');
X = double(T{:,idx});
if DETREND_MEAN, X = X - mean(X,1,'omitnan'); end

% compute bandpassed signals and spindle power in the window
sel = trel>=TIME_LIMS(1) & trel<=TIME_LIMS(2);
P_lin = zeros(1,3); P_dB = zeros(1,3);
XF = zeros(nnz(sel),3);

for k=1:3
    x  = X(:,k);
    xf = bandpass(x, BAND_HZ, fs,'Steepness',0.85,'StopbandAttenuation',60);  % zero-phase
    XF(:,k) = xf(sel);
    P_lin(k) = mean(xf(sel).^2,'omitnan');        % µV^2
    P_dB(k)  = 10*log10(P_lin(k)+eps);            % dB(µV^2)
end

% C4 plot: bandpassed trace + Hilbert envelope with shaded detections 
kC4 = find(strcmpi(CHS,'C4'),1,'first');  assert(~isempty(kC4));
xC4_bp = XF(:,kC4);
envC4  = abs(hilbert(xC4_bp));
tk     = trel(sel);

figure('Color','w','Units','inches','Position',[1 1 12 4.2]); hold on;
plot(tk, xC4_bp, 'b','LineWidth',1.0);
plot(tk, envC4,  'r','LineWidth',1.6);
yline(ENV_THR_UV,'k--','3 \muV','LabelHorizontalAlignment','left');

% shade contiguous segments above threshold with duration in range
above = envC4 > ENV_THR_UV;
d = diff([false; above; false]);
starts = find(d==1); ends = find(d==-1)-1;
ym = max([abs(xC4_bp); envC4])+1;
for i=1:numel(starts)
    dur = (ends(i)-starts(i)+1)/fs;
    if dur>=DUR_RANGE_S(1) && dur<=DUR_RANGE_S(2)
        a = tk(starts(i)); b = tk(ends(i));
        patch([a b b a],[-ym -ym ym ym],[1 0.8 0.8], ...
              'EdgeColor','none','FaceAlpha',0.35);
    end
end
grid on; box on;
xlabel('Time (s)');
ylabel('C4 (11–16 Hz) & envelope (\muV)');
title(sprintf('C4 spindle segment — %s (%.1f–%.1f s), %d–%d Hz', ...
      fn, TIME_LIMS(1), TIME_LIMS(2), BAND_HZ(1), BAND_HZ(2)));
legend({'11–16 Hz','Envelope','3 \muV','Detected (0.5–2 s)'},'Location','northeast');

% print spindle-band power (µV^2 and dB) for all channels 
fprintf('\nSpindle-band power %d–%d Hz in %.1f–%.1f s:\n', BAND_HZ(1), BAND_HZ(2), TIME_LIMS(1), TIME_LIMS(2));
for k=1:3
    fprintf('  %-3s: %9.3f µV^2   (%6.2f dB)\n', CHS{k}, P_lin(k), P_dB(k));
end
set(gcf,'Renderer','painters');
print(gcf, fullfile(pwd,'C4_spindle_segment.eps'), '-depsc','-r300');
