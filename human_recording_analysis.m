% C4 is NeuroWeaves here
% Import all extracted CSVs into one MATLAB struct
baseDir = ''; % add your path
files = dir(fullfile(baseDir, '*.csv'));

S = struct('name',{},'T',{},'time',{},'data',{},'vars',{},'fs',{});

for k = 1:numel(files)
    fpath = fullfile(files(k).folder, files(k).name);
    T = readtable(fpath, 'TextType','string');

    vars = T.Properties.VariableNames;
    vlow = lower(vars);

    tidx = find(contains(vlow,'time'), 1, 'first');
    if isempty(tidx), t = (0:height(T)-1)'; else, t = T{:,tidx}; end
    t = t(:);

    want = {'fp2','c4','o2'};
    data = nan(numel(t), numel(want));
    for j = 1:numel(want)
        cidx = find(strcmpi(vars, want{j}), 1, 'first');
        if isempty(cidx), data(:,j) = NaN; else, data(:,j) = T{:,cidx}; end
    end

    if numel(t) > 2
        dt = median(diff(t), 'omitnan');
        fs_est = 1/dt;
    else
        fs_est = NaN;
    end

    S(k).name = erase(files(k).name, '.csv');
    S(k).T = T;
    S(k).time = t;
    S(k).data = data; % [Fp2, C4, O2]
    S(k).vars = vars;
    S(k).fs = fs_est;
end

save(fullfile(baseDir, 'extracted_all.mat'), 'S');
disp("Loaded files:");
disp(string({S.name})');

% Choose baseline file
allNames = string({files.name});
iBase = find(contains(lower(allNames),'baseline'),1,'first');
if isempty(iBase), error('baseline*.csv not found'); end
inFile = fullfile(baseDir, allNames(iBase));

%% Plot 1: Baseline stacked (y=[-100,100] µV), export PNG/EPS
T = readtable(inFile,'TextType','string');
vars = T.Properties.VariableNames;
t = T{:, find(contains(lower(vars),'time'),1,'first')};
ch = {'Fp2','C4','O2'};
idx = cellfun(@(c) find(strcmpi(vars,c),1,'first'), ch);
X = [T{:,idx(1)}, T{:,idx(2)}, T{:,idx(3)}];

f = figure('Color','w','Units','inches','Position',[1 1 10 6]);
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');
for i = 1:3
    nexttile;
    plot(t, X(:,i), 'LineWidth',1.2);
    ylim([-100 100]);
    ylabel(ch{i});
    grid on;
    if i<3, set(gca,'XTickLabel',[]); end
end
xlabel('Time (s)');
title('Baseline | Fp2, C4, O2 (y = [-100, 100] \muV)');
outDir = fileparts(inFile);
exportgraphics(f, fullfile(outDir,'baseline_stacked_fixed.png'), 'Resolution',300);
set(f,'Renderer','painters');
exportgraphics(f, fullfile(outDir,'baseline_stacked_fixed.eps'), 'ContentType','vector');

%% Plot 2: Rotated PSD (C4 vs O2) to 100 Hz, export EPS
xC4 = T{:, find(strcmpi(vars,'C4'),1,'first')};
xO2 = T{:, find(strcmpi(vars,'O2'),1,'first')};
fs = 1/median(diff(t));
L = numel(xC4);
wlen = min(round(4*fs), L);
wlen = max(64, wlen);
nover = floor(0.5*wlen);
nfft = 2^nextpow2(max(wlen, min(L, round(8*fs))));
win = hamming(wlen);
[PC,fpsd] = pwelch(double(xC4), win, nover, nfft, fs);
[PO,~ ] = pwelch(double(xO2), win, nover, nfft, fs);
band = (fpsd>=0.5 & fpsd<=100);

f2 = figure('Color','w','Units','inches','Position',[1 1 6 8]);
plot(10*log10(PC(band)), fpsd(band), 'LineWidth',1.4);
hold on;
plot(10*log10(PO(band)), fpsd(band), 'LineWidth',1.4);
grid on; box off;
ylim([0.5 100]);
xlabel('PSD (dB/\muV^2/Hz)');
ylabel('Frequency (Hz)');
legend({'C4 (bespoke)','O2'},'Location','southeast');
title('Baseline PSD (rotated): C4 vs O2');
set(f2,'Renderer','painters');
exportgraphics(f2, fullfile(outDir,'baseline_psd_rotated_C4_vs_O2.eps'), 'ContentType','vector');

%% Plot 3: PDF histograms (C4 vs O2), shaded, Δ=2 µV, export EPS
xC4 = T{:, find(strcmpi(vars,'C4'),1,'first')};
xO2 = T{:, find(strcmpi(vars,'O2'),1,'first')};
allx = [xC4; xO2];
binw = 2; % Δ = 2 µV
edges = floor(min(allx)/binw)*binw : binw : ceil(max(allx)/binw)*binw;

f3 = figure('Color','w','Units','inches','Position',[1 1 7.5 4.5]); hold on;
histogram(xC4,'BinEdges',edges,'Normalization','pdf','FaceColor',[0 0.45 0.74],'FaceAlpha',0.35,'EdgeColor','none');
histogram(xO2,'BinEdges',edges,'Normalization','pdf','FaceColor',[0.75 0.2 0.2],'FaceAlpha',0.35,'EdgeColor','none');
grid on; box off;
xlabel('\muV'); ylabel('PDF (1/\muV)');
legend({'C4 (bespoke)','O2'},'Location','northeast');
title('Baseline amplitude PDFs (Δ = 2 \muV; prob \approx PDF \times 2 \muV)');
set(f3,'Renderer','painters');
exportgraphics(f3, fullfile(outDir,'baseline_pdf_shaded.eps'), 'ContentType','vector');

% Levene’s test (variance equality) for C4 vs O2
grp = [repmat("C4",numel(xC4),1); repmat("O2",numel(xO2),1)];
X = [xC4; xO2];
p_lev = vartestn(X, grp, 'TestType','LeveneAbsolute','Display','off');
title(sprintf('Baseline PDFs (\\Delta=2 \\muV; prob \\approx PDF\\times2) | Levene p=%.3g', p_lev));

%% Plot 4 chewing (0–2 s window)
% --- locate and load the CHEWING file ---
iChew = find(contains(lower(allNames), 'chewing'), 1, 'first');
assert(~isempty(iChew), 'Could not find a file with "chewing" in its name.');
inFile = fullfile(baseDir, allNames(iChew));
T = readtable(inFile, 'TextType','string');

% --- time + channels ---
vars = T.Properties.VariableNames;
vlow = lower(string(vars));
tidx = find(contains(vlow,"time"), 1, "first");
assert(~isempty(tidx), 'No time-like column found in the chewing file.');

t = T{:,tidx}; t = t(:);
trel = t - t(1);

chIdx   = setdiff(1:numel(vars), tidx);   % every column except time
chNames = vars(chIdx);
nCh     = numel(chIdx);

%  plot all channels stacked with fixed y-limits [-400, 400] 
f = figure('Color','w','Units','inches','Position',[1 1 12 max(3,1.2*nCh)]);
tiledlayout(nCh,1,'Padding','compact','TileSpacing','compact');

for i = 1:nCh
    nexttile;
    plot(trel, T{:, chIdx(i)}, 'LineWidth', 1.1);
    ylim([-400 400]);
    grid on;
    ylabel(chNames{i}, 'Interpreter','none');
    if i < nCh, set(gca,'XTickLabel',[]); else, xlabel('Time (s)'); end
end
title(sprintf('Chewing | %s — all channels (raw)', erase(allNames(iChew), '.csv')));

% --- export as EPS (vector) ---
outName = sprintf('%s_allch_raw_ylim_400.eps', erase(allNames(iChew), '.csv'));
exportgraphics(f, fullfile(baseDir, outName), 'ContentType','vector');

%% Plot 4b chewing (0–0.5 s window)
% locate and load the CHEWING file
iChew = find(contains(lower(allNames), 'chewing'), 1, 'first');
assert(~isempty(iChew), 'Could not find a file with "chewing" in its name.');
inFile = fullfile(baseDir, allNames(iChew));
T = readtable(inFile, 'TextType','string');

%  time + channels 
vars = T.Properties.VariableNames;
vlow = lower(string(vars));
tidx = find(contains(vlow,"time"), 1, "first");
assert(~isempty(tidx), 'No time-like column found in the chewing file.');

t = T{:,tidx}; 
t = t(:);
trel = t - t(1);

% select 0–0.5 s window
sel = trel >= 0 & trel <= 0.5;

chIdx   = setdiff(1:numel(vars), tidx);   % every column except time
chNames = vars(chIdx);
nCh     = numel(chIdx);

%  plot all channels stacked with fixed y-limits [-400, 400] 
f = figure('Color','w','Units','inches','Position',[1 1 12 max(3,1.2*nCh)]);
tiledlayout(nCh,1,'Padding','compact','TileSpacing','compact');

for i = 1:nCh
    nexttile;
    plot(trel(sel), T{sel, chIdx(i)}, 'LineWidth', 1.1);
    ylim([-400 400]);
    grid on;
    ylabel(chNames{i}, 'Interpreter','none');
    if i < nCh
        set(gca,'XTickLabel',[]);
    else
        xlabel('Time (s)');
    end
end

title(sprintf('Chewing | %s — all channels (0–0.5 s, raw)', erase(allNames(iChew), '.csv')));

%  export as EPS 
outName = sprintf('%s_allch_raw_0to0p5s_ylim_400.eps', erase(allNames(iChew), '.csv'));
exportgraphics(f, fullfile(baseDir, outName), 'ContentType','vector');


%% Plot 5 PDR
% Ensure we have the list of CSV names
if ~exist('allNames','var') || isempty(allNames)
    assert(exist('baseDir','var')==1 && isfolder(baseDir), ...
        'baseDir is missing or not a folder.');
    files    = dir(fullfile(baseDir,'*.csv'));
    allNames = string({files.name});
end

% Locate blinking file
iBlink = find(contains(lower(allNames),'blinking'), 1, 'first');
assert(~isempty(iBlink), 'blinking*.csv not found');
inFile = fullfile(baseDir, allNames(iBlink));

% Read table
T    = readtable(inFile,'TextType','string');
vars = T.Properties.VariableNames;
vlow = lower(vars);

% Time (or fallback to sample index)
tidx = find(contains(vlow,'time'), 1, 'first');
if isempty(tidx)
    t   = (0:height(T)-1)'; 
    xlab = 'Samples';
else
    t   = T{:,tidx};
    xlab = 'Time (s)';
end
t    = t(:);
trel = t - t(1);
sel  = trel >= 0 & trel <= 16;   % 0–16 s window

% All numeric channels except time
isNum  = varfun(@(x) isnumeric(x), T, 'OutputFormat','uniform');
chIdx  = find(isNum);
if ~isempty(tidx), chIdx = setdiff(chIdx, tidx); end
chNames = vars(chIdx);
nCh     = numel(chIdx);

% Stacked plot (0–16 s) with fixed y-limits [-100, 100] µV
f = figure('Color','w','Units','inches','Position',[1 1 12 max(3,1.2*nCh)]);
tiledlayout(nCh,1,'Padding','compact','TileSpacing','compact');
for i = 1:nCh
    nexttile;
    plot(trel(sel), T{sel,chIdx(i)}, 'LineWidth', 1.1);
    ylim([-100 100]);                               % <-- fixed y-limits
    ylabel(chNames{i}, 'Interpreter','none'); grid on;
    if i < nCh, set(gca,'XTickLabel',[]); else, xlabel(xlab); end
end

% Title from filename
title(sprintf('Blinking | %s (0–16 s, y=[-100, 100] \\muV)', erase(allNames(iBlink), '.csv')));

% Export EPS
outDir = fileparts(inFile);
set(f,'Renderer','painters');
exportgraphics(f, fullfile(outDir,'blinking_stacked_0to16s_y-100to100.eps'), 'ContentType','vector');

%% Figure 5b (6.5–9.5 s window, all electrodes)

% --- Ensure blinking file is available ---
assert(exist('baseDir','var')==1 && isfolder(baseDir),'baseDir missing.');
if ~exist('allNames','var') || isempty(allNames)
    files = dir(fullfile(baseDir,'*.csv'));
    allNames = string({files.name});
end

iBlink = find(contains(lower(allNames),'blinking'),1,'first');
assert(~isempty(iBlink),'blinking*.csv not found');
inFile = fullfile(baseDir, allNames(iBlink));

%  Read table and time vector 
T    = readtable(inFile,'TextType','string');
vars = T.Properties.VariableNames;
vlow = lower(vars);
tidx = find(contains(vlow,'time'),1,'first');
assert(~isempty(tidx),'No time column found.');
t = T{:,tidx}; t = t(:);
trel = t - t(1);
xlab = 'Time (s)';

%  Select 6.5–9.5 s window 
sel = trel >= 6.5 & trel <= 9.5;

% Channels 
chs = {'Fp2','C4','O2'};
cols = [0.1 0.5 0.9; 0.9 0.4 0.1; 0.2 0.7 0.3];

%  Plot stacked traces 
f = figure('Color','w','Units','inches','Position',[1 1 10 6]);
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

for i = 1:numel(chs)
    cidx = find(strcmpi(vars,chs{i}),1,'first');
    assert(~isempty(cidx), sprintf('%s not found',chs{i}));
    x = double(T{:,cidx});
    nexttile;
    plot(trel(sel), x(sel), 'Color',cols(i,:), 'LineWidth',1.2);
    ylim([-100 100]);
    ylabel([chs{i} ' (\muV)']);
    grid on;
    if i < numel(chs)
        set(gca,'XTickLabel',[]);
    else
        xlabel(xlab);
    end
end

title('Blinking | Fp2, C4, O2 (6.5–9.5 s, y=[-100,100] \muV)');

% --- Export EPS ---
outDir = fileparts(inFile);
set(f,'Renderer','painters');
exportgraphics(f, fullfile(outDir,'blinking_allch_6p5to9p5s_y-100to100.eps'), 'ContentType','vector');



%% Plot 6 Coherence PDR
% Axes:  O2 (x) vs C4 (y) alpha power (dB)
% Color: Fp2 alpha power (dB)  [cool=weak, warm=strong]
% Size:  C4–O2 narrowband correlation^2 (proxy for coherence) in 8–13 Hz
% Arrows: one per-epoch mean (Closed → Open)

% --- time & channels ---
vars = T.Properties.VariableNames; vlow = lower(vars);
t    = T{:, find(contains(vlow,'time'),1,'first')}; t = t(:);
fs   = 1/median(diff(t));
trel = t - t(1);
xO2  = double(T{:, strcmpi(vars,'O2')});
xC4  = double(T{:, strcmpi(vars,'C4')});
xF2  = double(T{:, strcmpi(vars,'Fp2')});

% --- bandpass 8–13 Hz once (zero-phase) ---
xO2b = bandpass(xO2,[8 13],fs,'Steepness',0.85,'StopbandAttenuation',60);
xC4b = bandpass(xC4,[8 13],fs,'Steepness',0.85,'StopbandAttenuation',60);
xF2b = bandpass(xF2,[8 13],fs,'Steepness',0.85,'StopbandAttenuation',60);

% --- 1-s non-overlapping windows: Closed 0–8 s, Open 9–16 s ---
edgesC = 0:1:8;   nC = numel(edgesC)-1;   % 8 windows
edgesO = 9:1:16;  nO = numel(edgesO)-1;   % 7 windows

% helper to compute dB power and corr^2 per window
powdb = @(x,sel) 10*log10(mean(x(sel).^2) + eps);
coh2  = @(a,b,sel) max(0, corr(a(sel), b(sel), 'type','Pearson','rows','complete')^2);

% --- collect per-window metrics ---
O2_C = zeros(nC,1); C4_C = zeros(nC,1); F2_C = zeros(nC,1); Coh_C = zeros(nC,1);
for k=1:nC
    sel = trel>=edgesC(k) & trel<edgesC(k+1);
    O2_C(k)=powdb(xO2b,sel); C4_C(k)=powdb(xC4b,sel); F2_C(k)=powdb(xF2b,sel);
    Coh_C(k)=coh2(xO2b,xC4b,sel);
end
O2_O = zeros(nO,1); C4_O = zeros(nO,1); F2_O = zeros(nO,1); Coh_O = zeros(nO,1);
for k=1:nO
    sel = trel>=edgesO(k) & trel<edgesO(k+1);
    O2_O(k)=powdb(xO2b,sel); C4_O(k)=powdb(xC4b,sel); F2_O(k)=powdb(xF2b,sel);
    Coh_O(k)=coh2(xO2b,xC4b,sel);
end

% --- plot ---
f = figure('Color','w','Units','inches','Position',[1 1 7.2 5.6]); hold on;
% identity line
xx = linspace(min([O2_C;O2_O])-1, max([O2_C;O2_O])+1, 200);
plot(xx, xx, 'k--', 'LineWidth',1);

% sizes: scale coherence^2 into marker area
sC = 40 + 260*Coh_C;   % closed
sO = 40 + 260*Coh_O;   % open

% scatter: Closed (circles), Open (squares). Color = Fp2 power (dB)
sc1 = scatter(O2_C, C4_C, sC, F2_C, 'o', 'filled', 'MarkerFaceAlpha',0.85, 'DisplayName','Closed 0–8 s');
sc2 = scatter(O2_O, C4_O, sO, F2_O, 's', 'filled', 'MarkerFaceAlpha',0.85, 'DisplayName','Open 9–16 s');

% per-epoch mean arrow (Closed → Open)
mc = [mean(O2_C) mean(C4_C)];  mo = [mean(O2_O) mean(C4_O)];
quiver(mc(1), mc(2), mo(1)-mc(1), mo(2)-mc(2), 0, 'k', 'LineWidth',1.5, 'MaxHeadSize',0.4, 'DisplayName','Mean shift');

% axis/labels
xlabel('O2 alpha power (dB)'); ylabel('C4 (bespoke) alpha power (dB)');
title('Alpha concordance: O2 vs C4 | Color = Fp2 dB, Size = (C4–O2) coherence^2');
grid on; box on; axis tight; axis square;

% --- make it TURBO ---
colormap(turbo);                          % <<<<<<<<<<<<<<<<

cb = colorbar; cb.Label.String = 'Fp2 alpha power (dB)';
legend('Location','southeast');

%  export
outDir = fileparts(inFile); if isempty(outDir), outDir = pwd; end
set(f,'Renderer','painters');
exportgraphics(f, fullfile(outDir,'tri_channel_bubble_O2_vs_C4_alpha.eps'), 'ContentType','vector');

%% Plot 7 closed vs open power
% timebase
vars = T.Properties.VariableNames; vlow = lower(string(vars));
t    = T{:, find(contains(vlow,'time'),1,'first')}; t = t(:);
fs   = 1/median(diff(t));
trel = t - t(1);

% channels
O2  = double(T{:, find(strcmpi(vars,'O2'), 1)});
C4  = double(T{:, find(strcmpi(vars,'C4'), 1)});
Fp2 = double(T{:, find(strcmpi(vars,'Fp2'),1)});

% bandpass once (8–13 Hz) for alpha power
O2b  = bandpass(O2, [8 13], fs, 'Steepness',0.85,'StopbandAttenuation',60);
C4b  = bandpass(C4, [8 13], fs, 'Steepness',0.85,'StopbandAttenuation',60);
Fp2b = bandpass(Fp2,[8 13], fs, 'Steepness',0.85,'StopbandAttenuation',60);

% windows
edgesC = 0:1:8;   nC = numel(edgesC)-1;   % closed 8 bins
edgesO = 9:1:16;  nO = numel(edgesO)-1;   % open   7 bins
powdb = @(x,sel) 10*log10(mean(x(sel).^2) + eps);

% per-window alpha power (dB)
O2_C = zeros(nC,1); O2_O = zeros(nO,1);
C4_C = zeros(nC,1); C4_O = zeros(nO,1);
F2_C = zeros(nC,1); F2_O = zeros(nO,1);
for k=1:nC
    sel = trel>=edgesC(k) & trel<edgesC(k+1);
    O2_C(k)=powdb(O2b,sel); C4_C(k)=powdb(C4b,sel); F2_C(k)=powdb(Fp2b,sel);
end
for k=1:nO
    sel = trel>=edgesO(k) & trel<edgesO(k+1);
    O2_O(k)=powdb(O2b,sel); C4_O(k)=powdb(C4b,sel); F2_O(k)=powdb(Fp2b,sel);
end

% labels for stats/ROC
y_closed =  ones(nC,1);   % 1 = closed
y_open   =  zeros(nO,1);  % 0 = open

%% (1) Paired estimation-style plot: Closed vs Open alpha power (per channel)
% Shows distributions, mean difference with bootstrap CI, and permutation p

Bperm = 10000;  Bboot = 10000;
chs = {'O2','C4','Fp2'};
Xc  = {O2_C, C4_C, F2_C};
Xo  = {O2_O, C4_O, F2_O};

f = figure('Color','w','Units','inches','Position',[1 1 9 3.6]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

for ci_idx = 1:3
    nexttile; hold on;
    xc = Xc{ci_idx}; xo = Xo{ci_idx};
    % ensure numeric & finite
    xc = xc(isfinite(xc)); xo = xo(isfinite(xo));

    % --- draw CATEGORICAL boxes FIRST (so axis is categorical) ---
    g = [repmat("open",numel(xo),1); repmat("closed",numel(xc),1)];
    y = [xo; xc];
    boxchart(categorical(g), y, ...
             'BoxFaceColor',[.9 .9 .9], 'WhiskerLineColor','k', ...
             'LineWidth',0.75, 'BoxFaceAlpha',0.25);

    % --- overlay jittered points using categorical x ----
    swarmchart(categorical(repmat("open",  numel(xo),1)), xo, 22, [0.25 0.55 0.90], 'filled', ...
               'MarkerFaceAlpha',0.7);
    swarmchart(categorical(repmat("closed",numel(xc),1)), xc, 22, [0.95 0.60 0.25], 'filled', ...
               'MarkerFaceAlpha',0.7);

    % --- permutation p (difference in means) ---
    labs = [zeros(numel(xo),1); ones(numel(xc),1)];   % 0=open, 1=closed
    obs  = mean(xc) - mean(xo);
    cnt  = 0;
    for b = 1:Bperm
        ls = labs(randperm(numel(labs)));
        cnt = cnt + (abs(mean(y(ls==1)) - mean(y(ls==0))) >= abs(obs));
    end
    p_perm = (cnt+1)/(Bperm+1);

    % --- bootstrap CI for mean difference (closed - open) ---
    diffs = zeros(Bboot,1);
    for b = 1:Bboot
        xb_c = xc(randi(numel(xc), numel(xc), 1));
        xb_o = xo(randi(numel(xo), numel(xo), 1));
        diffs(b) = mean(xb_c) - mean(xb_o);
    end
    ci95 = quantile(diffs,[0.025 0.975]);  % <-- renamed (was clobbering loop var)

    ylabel('Alpha power (dB)');
    title(sprintf('%s: \\Delta=%.2f dB, 95%% CI [%.2f, %.2f]\\np=%.3g (perm)', ...
          chs{ci_idx}, obs, ci95(1), ci95(2), p_perm));
    grid on; box on;
end

%% Plot 8 - closed vs open power for C4
% Assumes you've already run the setup that defines C4_C (closed) and C4_O (open).

Bperm = 10000;   % permutations
Bboot = 10000;   % bootstrap resamples

xc = C4_C;  xo = C4_O;                 % closed, open
xc = xc(isfinite(xc)); xo = xo(isfinite(xo));

% permutation p (difference in means, closed - open)
y    = [xo; xc];
labs = [zeros(numel(xo),1); ones(numel(xc),1)];   % 0=open, 1=closed
obs  = mean(xc) - mean(xo);
cnt  = 0;
for b = 1:Bperm
    ls  = labs(randperm(numel(labs)));
    cnt = cnt + (abs(mean(y(ls==1)) - mean(y(ls==0))) >= abs(obs));
end
p_perm = (cnt+1)/(Bperm+1);

% bootstrap 95% CI for mean difference
diffs = zeros(Bboot,1);
for b = 1:Bboot
    xb_c = xc(randi(numel(xc), numel(xc), 1));
    xb_o = xo(randi(numel(xo), numel(xo), 1));
    diffs(b) = mean(xb_c) - mean(xb_o);
end
ci95 = quantile(diffs, [0.025 0.975]);

% plot C4 only
f = figure('Color','w','Units','inches','Position',[1 1 3.6 5]); hold on;
boxchart(categorical([repmat("open",numel(xo),1); repmat("closed",numel(xc),1)]), ...
         [xo; xc], 'BoxFaceColor',[.9 .9 .9], 'WhiskerLineColor','k', ...
         'LineWidth',0.75, 'BoxFaceAlpha',0.25);

swarmchart(categorical(repmat("open",  numel(xo),1)), xo, 22, [0.25 0.55 0.90], 'filled', 'MarkerFaceAlpha',0.7);
swarmchart(categorical(repmat("closed",numel(xc),1)), xc, 22, [0.95 0.60 0.25], 'filled', 'MarkerFaceAlpha',0.7);

ylabel('Alpha power (dB)');
title(sprintf('C4: \\Delta=%.2f dB, 95%% CI [%.2f, %.2f]\\np=%.4f (perm)', ...
      obs, ci95(1), ci95(2), p_perm));
grid on; box on;

% export EPS
if exist('inFile','var') && ~isempty(inFile)
    outDir = fileparts(inFile);
else
    outDir = pwd;
end
set(f,'Renderer','painters');
exportgraphics(f, fullfile(outDir, 'C4_alpha_closed_vs_open_estimation.eps'), 'ContentType','vector');

%% Plot 9 Photic: C4 traces (first 5 s) + scalograms (first 5 s)
% - Focus channel: C4
% - One figure per photic file: top = first 5 s trace, bottom = scalogram of same segment
% - Adjustable scalogram settings via SC_FREQ_LIMS and SC_CLIM
% - No extra annotations

clearvars -except baseDir; clc;

% ---------- Adjustable ----------
SC_FREQ_LIMS = [1 40];     % Hz range for scalogram
SC_CLIM      = [0 25];     % dB color scale (adjust as needed)
TRACE_DUR_S  = 5.0;        % plot first 5 seconds

% ---------- Folder ----------
if ~exist('baseDir','var') || ~isfolder(baseDir)
    baseDir = uigetdir(pwd,'Select folder with your CSVs');
    if isequal(baseDir,0), error('No folder selected.'); end
end

files = dir(fullfile(baseDir,'*.csv'));
names = string({files.name});
low   = lower(names);

% keep only photic files
keep = contains(low,"photic");
assert(any(keep),'No photic*.csv files found in the folder.');
names = names(keep);

for k = 1:numel(names)
    % ----- Load -----
    T    = readtable(fullfile(baseDir, names(k)), 'TextType','string');
    vars = T.Properties.VariableNames;
    vlow = lower(vars);

    % Time column (fallback to sample index)
    tidx = find(contains(vlow,'time'),1,'first');
    if isempty(tidx)
        t  = (0:height(T)-1)';    % samples
        fs = 1;                   % unknown; only used for labels if no time
        xlab = 'Samples';
    else
        t  = T{:,tidx};
        fs = 1/median(diff(t));
        xlab = 'Time (s)';
    end
    t    = t(:);
    trel = t - t(1);

    % C4 channel
    cidx = find(strcmpi(vars,'C4'),1,'first');
    assert(~isempty(cidx),'C4 column not found in %s', names(k));
    x = double(T{:,cidx}); x = x(:);

    % ----- Select first 5 seconds (or all if shorter) -----
    t0 = trel(1);
    t1 = min(t0 + TRACE_DUR_S, trel(end));
    sel = trel>=t0 & trel<=t1;

    % ----- Compute CWT on the selected window only -----
    % To keep time axes aligned to the selected window, compute on full then slice
    [WT, F] = cwt(x, fs, 'FrequencyLimits', SC_FREQ_LIMS);
    PdB = 10*log10(abs(WT).^2 + eps);
    [~,i0] = min(abs(trel - t0));
    [~,i1] = min(abs(trel - t1));
    twin   = trel(i0:i1);
    PdBw   = PdB(:, i0:i1);

    % ----- Figure: top trace, bottom scalogram (no annotations) -----
    fig = figure('Color','w','Units','inches','Position',[1 1 10 5]); %#ok<NASGU>
    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    % Top: first 5 s raw trace (C4)
    nexttile;
    plot(trel(sel), x(sel), 'LineWidth',1.0);
    box on; grid on;
    xlabel(xlab); ylabel('\muV');
    title(sprintf('%s — C4 (first %.1f s)', names(k), TRACE_DUR_S), 'Interpreter','none');

    % Bottom: scalogram over same window
    nexttile;
    imagesc(twin, F, PdBw); axis xy; box on; grid on;
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    title('Scalogram (CWT power, C4)');
    cb = colorbar; cb.Label.String = 'Power (dB)';
    caxis(SC_CLIM);  % 0–25 dB (adjust via SC_CLIM)

    % Optional export:
    % outpng = fullfile(baseDir, sprintf('%s_C4_first5s_trace_scalo.png', erase(names(k),'.csv')));
    % exportgraphics(gcf, outpng, 'Resolution',300);
end

%% Plot 10 Stacked PSD (10 - 30 Hz)
% Welch PSD (4 s Hamming, 50% overlap), 0.5–40 Hz. No extra annotations.

clearvars -except baseDir; clc;

% ----- settings -----
FMIN = 0.5; FMAX = 40;      % frequency axis (Hz)
WIN_SEC = 4; OL_FRAC = 0.5; % Welch params

% ----- folder -----
if ~exist('baseDir','var') || ~isfolder(baseDir)
    baseDir = uigetdir(pwd,'Select folder with your CSVs'); 
    if isequal(baseDir,0), error('No folder selected.'); end
end
files = dir(fullfile(baseDir,'*.csv')); 
names = string({files.name}); 
low   = lower(names);

% ----- find files for 10/20/30 Hz -----
freqs = [10 20 30];
idx = nan(size(freqs));
for k = 1:numel(freqs)
    ii = find(contains(low,'photic') & contains(low,string(freqs(k))) & contains(low,'hz'), 1, 'first');
    assert(~isempty(ii), 'Missing photic %d Hz file.', freqs(k));
    idx(k) = ii;
end

% ----- compute PSDs and collect y-lims -----
F_all  = cell(1,3); 
C4_all = cell(1,3); 
O2_all = cell(1,3);
ymin = inf; ymax = -inf;

for k = 1:3
    T  = readtable(fullfile(baseDir, names(idx(k))), 'TextType','string');
    vn = T.Properties.VariableNames; vnl = lower(vn);

    ti = find(contains(vnl,'time'),1,'first'); 
    t  = T{:,ti}; t = t(:); 
    fs = 1/median(diff(t));

    iC4 = find(strcmpi(vn,'C4'),1,'first'); 
    iO2 = find(strcmpi(vn,'O2'),1,'first');
    xC4 = double(T{:,iC4}) - mean(T{:,iC4},'omitnan');
    xO2 = double(T{:,iO2}) - mean(T{:,iO2},'omitnan');

    win  = hamming(max(64, round(WIN_SEC*fs)));
    ol   = round(OL_FRAC*numel(win));
    nfft = 2^nextpow2(max(numel(win), round(8*fs)));

    [PC,F] = pwelch(xC4, win, ol, nfft, fs, 'psd');
    [PO,~] = pwelch(xO2, win, ol, nfft, fs, 'psd');

    band = (F>=FMIN & F<=FMAX); F = F(band);
    C4dB = 10*log10(PC(band)+eps);
    O2dB = 10*log10(PO(band)+eps);

    F_all{k}  = F; 
    C4_all{k} = C4dB; 
    O2_all{k} = O2dB;

    ymin = min([ymin; C4dB; O2dB]); 
    ymax = max([ymax; C4dB; O2dB]);
end
YMIN = floor(ymin); YMAX = ceil(ymax);

% ----- plot stacked panels -----
figure('Color','w','Units','inches','Position',[1 1 9 8]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
cols = [0 0.45 0.74; 0.85 0.33 0.10]; % O2 blue, C4 orange

for k = 1:3
    nexttile; hold on;
    plot(F_all{k}, O2_all{k}, 'Color', cols(1,:), 'LineWidth',1.8); % O2
    plot(F_all{k}, C4_all{k}, 'Color', cols(2,:), 'LineWidth',1.8); % C4
    grid on; box on; xlim([FMIN FMAX]); ylim([YMIN YMAX]);
    title(sprintf('Photic %d Hz — PSD (O2 vs C4)', freqs(k)));
    if k < 3
        set(gca,'XTickLabel',[]);
    else
        xlabel('Frequency (Hz)');
    end
    ylabel('dB \muV^2/Hz');
    if k == 1
        legend({'O2','C4'},'Location','northeast');
    end
end

exportgraphics(gcf, fullfile(baseDir,'PSD_stacked_O2_vs_C4_10_20_30Hz.eps'), 'Resolution',300);

%% Plot 11: SNR @ f for photic 10/20/30 Hz — C4 vs O2 
% Reuses Welch PSD approach; averages SNRs when multiple files per frequency exist.

clearvars -except baseDir; clc;

% ---- Settings (tweak if needed) ----
FMIN = 0.5; FMAX = 40;          % PSD band (Hz)
WIN_SEC_DEFAULT = 4;            % target Welch window (s) – auto-shrinks if needed
OL_FRAC = 0.5;                  % 50% overlap
SB_INNER = 1; SB_OUTER = 3;     % sidebands around f: [f-3,f-1] U [f+1,f+3]
LOWER_NOISE_CLIP = 2;           % exclude <2 Hz from noise
TARGET_FREQS = [10 20 30];      % compute only these

% ---- Folder ----
if ~exist('baseDir','var') || ~isfolder(baseDir)
    baseDir = uigetdir(pwd,'Select folder with your CSVs');
    if isequal(baseDir,0), error('No folder selected.'); end
end

dd = dir(fullfile(baseDir,'*.csv'));
names = string({dd.name});
low   = lower(names);

% ---- Collect matching photic files (10/20/30 Hz) ----
file_list = strings(0,1);
f0_list   = [];
for i = 1:numel(names)
    tok = regexp(low(i), 'photic[^0-9]*([0-9]+)\s*hz', 'tokens', 'once');
    if ~isempty(tok)
        f0 = str2double(tok{1});
        if any(TARGET_FREQS == f0)
            file_list(end+1,1) = names(i); %#ok<AGROW>
            f0_list(end+1,1)   = f0;       %#ok<AGROW>
        end
    end
end
assert(~isempty(file_list), 'No photic 10/20/30 Hz files found.');

% ---- Compute SNR per file for C4 and O2 ----
snrC = nan(numel(file_list),1);
snrO = nan(numel(file_list),1);

for i = 1:numel(file_list)
    T  = readtable(fullfile(baseDir, file_list(i)), 'TextType','string');
    vn = T.Properties.VariableNames; vnl = lower(vn);

    % time / fs
    ti = find(contains(vnl,'time'),1,'first'); assert(~isempty(ti),'No time column.');
    t  = T{:,ti}; t = t(:); fs = 1/median(diff(t));

    % channels
    iC4 = find(strcmpi(vn,'C4'),1,'first');
    iO2 = find(strcmpi(vn,'O2'),1,'first');
    assert(~isempty(iC4) && ~isempty(iO2),'Need C4 and O2.');
    xC4 = double(T{:,iC4}) - mean(T{:,iC4},'omitnan');
    xO2 = double(T{:,iO2}) - mean(T{:,iO2},'omitnan');

    % robust Welch params per file
    Lx = numel(xC4);
    winlen = min(round(WIN_SEC_DEFAULT*fs), max(32, floor(Lx/4)));
    if winlen > Lx, winlen = max(16, floor(Lx/2)); end
    win = hamming(winlen);
    ol  = floor(OL_FRAC*winlen); if ol >= winlen, ol = floor(0.5*winlen); end
    nfft = 2^nextpow2(max(winlen, round(4*fs)));

    % PSDs
    try
        [PC,F] = pwelch(xC4, win, ol, nfft, fs, 'psd');
        [PO,~] = pwelch(xO2, win, ol, nfft, fs, 'psd');
    catch
        [PC,F] = periodogram(xC4, hann(min(Lx,1024)), nfft, fs, 'psd');
        [PO,~] = periodogram(xO2, hann(min(Lx,1024)), nfft, fs, 'psd');
    end

    % limit band
    sel = (F>=FMIN & F<=FMAX);
    F  = F(sel); PC = PC(sel); PO = PO(sel);

    % SNR at center frequency f0
    f0 = f0_list(i);
    [~,k0] = min(abs(F - f0));
    PsigC = PC(k0); PsigO = PO(k0);

    lb1 = max(LOWER_NOISE_CLIP, f0 - SB_OUTER);
    ub1 = max(LOWER_NOISE_CLIP, f0 - SB_INNER);
    lb2 = min(FMAX, f0 + SB_INNER);
    ub2 = min(FMAX, f0 + SB_OUTER);
    sb  = (F>=lb1 & F<ub1) | (F>lb2 & F<=ub2);
    if ~any(sb)  % fallback
        sb = (F>=max(LOWER_NOISE_CLIP,f0-5) & F<=min(FMAX,f0-2)) | ...
             (F>=min(FMAX,f0+2) & F<=min(FMAX,f0+5));
    end
    NmedC = median(PC(sb)); NmedO = median(PO(sb));

    snrC(i) = 10*log10(PsigC / max(NmedC, eps));
    snrO(i) = 10*log10(PsigO / max(NmedO, eps));
end

%  Average duplicates by frequency (ensures unique X for bars)
[uf, ~, g] = unique(f0_list);         % unique freqs (numeric)
meanC = accumarray(g, snrC, [], @mean);
meanO = accumarray(g, snrO, [], @mean);

% Bar chart (stacked panels avoided; single grouped bar with numeric xticks) 
figure('Color','w','Units','inches','Position',[1 1 8 4.5]);
bar([meanC meanO], 'grouped'); grid on; box on;
xticks(1:numel(uf)); xticklabels(string(uf));
xlabel('Flash frequency (Hz)'); ylabel('SNR @ f (dB)');
legend({'C4','O2'}, 'Location','northwest');
title('SNR at drive frequency — C4 vs O2 (mean per frequency)');

%  Console table
T = table(uf, meanC, meanO, 'VariableNames',{'Hz','SNR_C4_dB','SNR_O2_dB'});
disp(T);

exportgraphics(gcf, fullfile(baseDir,'SNR.eps'), 'Resolution',300);

%% Plot 12 Scalogram Chewing C4
%% Chewing: C4 scalogram (0–2 s) — extend freq to 100 Hz (turbo colormap)
SC_CLIM = [0 20];   % dB color scale ([] = auto)

% ---- timebase/fs from your existing vars ----
if exist('trel','var') && ~isempty(trel)
    t_rel = trel(:);
elseif exist('t','var') && ~isempty(t)
    t_rel = t(:) - t(1);
else
    error('No time vector found. Run the CHEWING load block first.');
end
if exist('t','var') && ~isempty(t)
    fs = 1/median(diff(t));
else
    fs = 1/median(diff(t_rel));
end

%  get C4 from your existing table T
assert(exist('T','var')==1 && ~isempty(T), 'Table T is missing. Run the CHEWING load block first.');
cidx = find(strcmpi(T.Properties.VariableNames,'C4'), 1, 'first');
assert(~isempty(cidx), 'C4 column not found in the CHEWING file.');
x = double(T{:,cidx}); x = x(:);
x = x - mean(x,'omitnan');

%  0–2 s window
sel = (t_rel >= 0) & (t_rel <= 2);
assert(nnz(sel) > 10, 'Too few samples between 0–2 s.');
twin = t_rel(sel);
xwin = x(sel);

% CWT up to 100 Hz 
fmax = min(100, (fs/2) - 0.1);
fmax = max(1, fmax);
[WT, F] = cwt(xwin, fs, 'FrequencyLimits', [1 fmax]);
PdB = 10*log10(abs(WT).^2 + eps);

% plot 
figure('Color','w','Units','inches','Position',[1 1 8 4.2]);
imagesc(twin, F, PdB); axis xy; grid on; box on;
colormap(turbo);                                 % <<< use turbo colormap
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title(sprintf('Chewing — C4 scalogram (0–2 s, 1–%.1f Hz)', fmax));
cb = colorbar; cb.Label.String = 'Power (dB)';
if ~isempty(SC_CLIM), caxis(SC_CLIM); end

exportgraphics(gcf, fullfile(baseDir,'chewing_C4_scalogram.eps'), 'ContentType','vector');


