%% Test effect of window size on power

% Make some noisey sine waves
len = 50;
t = 0:0.001:len;
nTrial = 1000;
srate = 1/0.001;

ts = zeros(nTrial, size(t,2));
for i = 1:nTrial
    y=sin(2*pi*(4)*t) + rand(size(t));
    data.trial{i} = y;
    data.time{i} = t;
end

data.fsample = srate;
data.label = {'test'};


% define windows
wins = 50:50:len*1000;
window_idx = repmat(wins, 1, nTrial/numel(wins))';
cfg = [];
cfg.begsample = ones(nTrial,1);
cfg.endsample = repmat(wins, 1, nTrial/numel(wins))';
curr = ft_redefinetrial(cfg, data);



cfg = [];
cfg.method = 'mtmfft'; % only single taper for low freqs
cfg.taper = 'hanning';
cfg.ouput = 'pow';
cfg.pad = 'nextpow2';
cfg.foi = 4;
cfg.keeptrials = 'yes';
%cfg.tapsmofrq = 4; %smoothing index - check if same effect is present for others

power = ft_freqanalysis(cfg,data);
pow_by_win = zeros(numel(wins), nTrial/numel(wins));
for i = 1:numel(wins)
    idx = window_idx == wins(i);
    pow_by_win(i,:) = power.powspctrm(idx);
end

% plot
figure(1); clf;
plot(wins, squeeze(mean(pow_by_win,2)),  "color", rgb("steelblue"), 'linewidth', 2); hold on
shade_plot(wins, squeeze(mean(pow_by_win,2))',  squeeze(std(pow_by_win,[],2))', rgb("slategrey"), 0.4);

