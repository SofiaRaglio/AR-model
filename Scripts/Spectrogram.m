% Computing the spectrogram of the LFP for the Monkey B and C.
% This script is used to produce Figure 1a and Suppl. Figure 1a of the paper. 

clearvars; close all;
% clearvars: remove all variables from the workspace (safer than clear all; avoids side effects)
% close all: close any open figures
% clc: clear the command window
%
% Practical motivation:
%   - When working with large files and time-frequency analyses it is easy to leave behind
%     residual variables that change results or consume RAM.

%% Load
load('PreProcData.mat');
load('MEAMUALFP.mat');
% Load two .mat files:
%   - PreProcData.mat: typically metadata, parameters, or preprocessed structures.
%   - MEAMUALFP.mat: signals (MEALFP.values), time vectors (MEALFP.time), and MEAMUA.time.
%
% Performance note:
%   - If MEAMUALFP.mat is ~5 GB, load() can be a bottleneck (I/O + RAM).
%   - matfile is not used here for chunked reading, so the entire file is loaded into memory.

%% Options
Options.SpecgramFreqs        = [0.1 5];
Options.SpecgramFreqSamples  = 100;      % number of frequency samples
Options.SpecgramMovWnd       = 5;        % seconds (moving-average window)
Options.BandWidth            = 0.05;      % relative bandwidth
SD_TIMES                     = 4.0;      % kernel half-width in sigmas
MA_OVERLAP                   = 0.5;      % 50% overlap
useSingle                    = true;     % speed/memory
% Key time-frequency parameters:
%   SpecgramFreqs = [0.1 5] Hz:
%       very slow oscillations (0.1–5 Hz), typical of slow LFP dynamics.
%   SpecgramFreqSamples = 100:
%       number of frequency bins; distributed linearly here (see frequency vector below).
%   SpecgramMovWnd = 5 s:
%       after convolution, apply a moving average over 5 seconds to summarize power in time.
%       This is controlled temporal downsampling (reduces noise and final matrix size).
%   BandWidth = 0.05:
%       relative wavelet bandwidth (time/frequency trade-off).
%       Smaller BandWidth -> longer wavelet -> better frequency resolution, worse time resolution.
%   SD_TIMES = 4:
%       truncate the Gaussian kernel at ±4 sigma; beyond that the Gaussian is negligible.
%       Motivation: finite, computable kernel with minimal truncation error.
%   MA_OVERLAP = 0.5:
%       moving-average windows with 50% overlap, i.e. hop = 0.5 * window.
%       Motivation: smoother temporal continuity without excessive redundancy.
%   useSingle = true:
%       convert data to single precision to save memory and often speed up FFT/convolutions.

%% Time base
Dt = 0.005;              % given
Fs = 1/Dt;
% Dt is the sample interval: 0.005 s -> Fs = 200 Hz.
% Fs is not used directly below (NASGU would suppress warnings), but is conceptually important:
%   - The 0.1–5 Hz band is well below Nyquist (Nyquist = 100 Hz).
%   - At 200 Hz the representation is adequate for slow LFP without the original 30 kHz rate.

%% Time segment (same for both regions)
t0      = MEAMUA.time(1);
t_start = MEAMUA.time(1000);
t_end   = MEAMUA.time(2040000);
% Select a time interval using indices 1000 and 2040000, as in the original script:
% likely to skip initial transients or select a stable portion of the recording.
% Convert absolute times (MEAMUA.time) to indices on MEALFP.time using Dt.

time_s = max(1, round((t_start - t0)/Dt) + 1);
time_e = min(size(MEALFP.values,2), round((t_end - t0)/Dt) + 1);
% Time -> index conversion:
%   - (t_start - t0)/Dt gives the number of samples after t0.
%   - round() maps to the nearest integer index.
%   - +1 because MATLAB indexing starts at 1.
% max/min guard against out-of-range indices.

t_vec = MEALFP.time(time_s:time_e);
t_vec = t_vec(:);
N     = numel(t_vec);
% t_vec is the time vector of the segment (column).
% N is the number of samples in the segment.

%% Frequency vector (log-spaced centers)
f1 = Options.SpecgramFreqs(1);
f2 = Options.SpecgramFreqs(2);
Nf = Options.SpecgramFreqSamples;

% Fedges = exp(linspace(log(f1), log(f2), Nf+1));
Fedges = linspace(f1, f2, Nf+1);
F = sqrt(Fedges(1:end-1) .* Fedges(2:end));
F = F(:);
% Build frequency bins:
%   - Fedges are bin edges (linear spacing here).
%   - F are bin centers, the geometric mean of adjacent edges.
% Motivation:
%   - the geometric mean is natural on a log scale (consistent with log binning);
%   - avoids using edges as effective frequencies (a common mistake), improving interpretation.

%% Build signals: prefrontal + premotor (mean across channels)
seg_pre = MEALFP.values(1:96,   time_s:time_e);
seg_pre = mean(seg_pre, 1).';
seg_pre = seg_pre(:);
% Extract prefrontal channels (1:96) over the time segment.
% mean(seg_pre,1): average across channels -> one representative regional trace.
% Motivation:
%   - reduce dimensionality (96 channels -> 1 signal) and uncorrelated noise.
%   - obtain a population-mean LFP activity per area.
% Trade-off:
%   - spatial detail across channels is lost, but robustness and speed improve.

seg_preM = MEALFP.values(97:192, time_s:time_e);
seg_preM = mean(seg_preM, 1).';
seg_preM = seg_preM(:);
% Same for premotor cortex (channels 97:192).

if useSingle
    seg_pre  = single(seg_pre);
    seg_preM = single(seg_preM);
end
% Convert to single to reduce RAM and speed up FFT operations.
% In spectral analyses, single is often sufficient because:
%   - numerical error is usually much smaller than biological variability,
%   - power is then averaged (moving average), further damping numerical noise.

%% Compute spectrograms (wavelet power + moving-average downsample)
[Spec_pre,  t_out] = wavelet_specgram_fast(seg_pre,  t_vec, F, Dt, Options.SpecgramMovWnd, MA_OVERLAP, Options.BandWidth, SD_TIMES, useSingle);
[Spec_prem, ~    ] = wavelet_specgram_fast(seg_preM, t_vec, F, Dt, Options.SpecgramMovWnd, MA_OVERLAP, Options.BandWidth, SD_TIMES, useSingle);
% Compute two wavelet scalograms:
%   - for each frequency F(k), convolve the signal with a complex wavelet
%   - take |z(t)|^2 as an estimate of instantaneous power at that frequency
%   - then average over temporal windows (5 s, 50% overlap)
%
% Theoretical motivation:
%   - The complex wavelet (Gaussian * complex sinusoid) is a band-pass filter centered at fk.
%   - Power |z|^2 is the local energy of the signal in the band around fk.
%   - Compared with classical STFT: the wavelet has adaptive resolution (longer kernel at low f),
%     which is often desirable for slow oscillations.

%% Plot: 2 panels (2 rows, 1 column), linked x-axes, bone flipped, y linear

Zpre  = log10(Spec_pre);
Zprem = log10(Spec_prem);

fig = figure;
ax1 = subplot(2,1,1);
contourf(t_out, F, Zpre, 10,'edgecolor','none');
ylabel('Frequency (Hz)');
title('Prefrontal (Channels 1:96)');
caxis(log10(exp([-10 -2])))
set(ax1, 'YScale', 'linear');  % requested

ax2 = subplot(2,1,2);
contourf(t_out, F, Zprem, 10,'edgecolor','none');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Premotor (Channels 97:192)');
caxis(log10(exp([-10 -2])))
set(ax2, 'YScale', 'linear');  % requested

linkaxes([ax1 ax2], 'x');

% colormap(fig, flipud(bone));
colormap(flipud(pink))

% ---- Colorbars with label + units ----
cb1 = colorbar(ax1);
cb2 = colorbar(ax2);

% Since you're plotting log10(power), label accordingly.
% If your underlying Spec_* is in "power" units (e.g., uV^2), adapt the unit string.
cb1.Label.String = 'log_{10} Power (a.u.)';
cb2.Label.String = 'log_{10} Power (a.u.)';

%% Plot: 2 panels (2 rows, 1 column), linked x-axes, bone flipped, y linear

% --- Data in log scale (as you already do) ---
Zpre  = log10(Spec_pre);
Zprem = log10(Spec_prem);

% --- Option: upsample via interpolation (NOT smoothing) to improve visual appearance ---
doInterp = true;     % set false to disable
interpT  = 4;        % upsample factor in time (e.g., 2-6)
interpF  = 3;        % upsample factor in freq (e.g., 1-4)

if doInterp
    tq = linspace(t_out(1), t_out(end), interpT*numel(t_out));
    Fq = linspace(F(1),     F(end),     interpF*numel(F));
    [Tg,Fg]   = meshgrid(t_out, F);
    [Tqi,Fqi] = meshgrid(tq,   Fq);

    Zpre  = interp2(Tg, Fg, Zpre,  Tqi, Fqi, 'linear');
    Zprem = interp2(Tg, Fg, Zprem, Tqi, Fqi, 'linear');

    tplot = tq;
    Fplot = Fq;
else
    tplot = t_out;
    Fplot = F;
end

fig = figure;

ax1 = subplot(2,1,1);
imagesc(tplot, Fplot, Zpre);
axis xy
ylabel('Frequency (Hz)');
title('Prefrontal (Channels 1:96)');
caxis(log10(exp([-10 -2])))
set(ax1, 'YScale', 'linear');

ax2 = subplot(2,1,2);
imagesc(tplot, Fplot, Zprem);
axis xy
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Premotor (Channels 97:192)');
caxis(log10(exp([-10 -2])))
set(ax2, 'YScale', 'linear');

linkaxes([ax1 ax2], 'x');

colormap(flipud(bone));   % or alternatively: flipud(pink)

cb1 = colorbar(ax1);
cb2 = colorbar(ax2);
cb1.Label.String = 'log_{10} Power (a.u.)';
cb2.Label.String = 'log_{10} Power (a.u.)';

% (optional) OpenGL renderer often improves image rendering
set(fig, 'Renderer', 'opengl');


%% -------- Local function (R2021b OK in scripts) --------
function [WSP, t_out] = wavelet_specgram_fast(x, t, F, Dt, MovingAverageWindow, MA_OVERLAP, RelBandWidth, SD_TIMES, useSingle)
% Core time-frequency computation.
% INPUT:
%   x: 1D signal (N x 1)
%   t: time vector (N x 1)
%   F: frequency vector (Nf x 1)
%   Dt: sample interval
%   MovingAverageWindow: window length (seconds) for power downsampling/smoothing
%   MA_OVERLAP: overlap between windows (0.5 -> hop = 0.5 * win)
%   RelBandWidth: relative bandwidth (controls wavelet sigma)
%   SD_TIMES: kernel truncation at ±SD_TIMES*sigma
%   useSingle: if true, use single precision
%
% OUTPUT:
%   WSP: Nf x SampleNum matrix of mean power per window
%   t_out: time vector (SampleNum x 1) for each window

    x = x(:);
    t = t(:);
    N = numel(x);
    Nf = numel(F);

    % Moving-average windowing (downsampled time grid)
    WindowSize = round(MovingAverageWindow / Dt);
    HopSize    = max(1, round((1 - MA_OVERLAP) * WindowSize)); % 50% overlap => hop=0.5*win
    starts     = 1:HopSize:(N - WindowSize + 1);
    SampleNum  = numel(starts);
    t_out      = t(starts + floor(WindowSize/2));
    % Define temporal windows:
    %   - WindowSize in samples (5 s -> 1000 samples with Dt=0.005).
    %   - HopSize: step between windows (overlap=0.5 -> hop ~500 samples).
    %   - starts: start index of each window.
    %   - t_out: time at window center.
    %
    % Motivation:
    %   - instead of keeping a power estimate at every sample (huge N),
    %     compress to SampleNum windows: much less RAM and cleaner plots.

    % Preallocate
    if useSingle
        WSP = zeros(Nf, SampleNum, 'single');
        x_ = single(x);
    else
        WSP = zeros(Nf, SampleNum);
        x_ = double(x);
    end
    % Preallocation is crucial in MATLAB:
    %   - avoids repeated reallocations and memory fragmentation inside the frequency loop.

    % FFT sizing based on longest kernel (set by minimum frequency)
    Sigma_minF = 1/(2*pi*min(F)*RelBandWidth);          % IMPORTANT: corrected scaling
    maxHalfLen = ceil(SD_TIMES * Sigma_minF / Dt);
    maxWL      = 2*maxHalfLen + 1;
    NFFT       = 2^nextpow2(N + maxWL - 1);
    % Choose FFT length:
    %   - FFT convolution needs length at least N + M - 1 (M = kernel length).
    %   - The longest kernel is at the minimum frequency (0.1 Hz) because sigma ~ 1/f.
    %   - Round up to a power of 2 for efficiency.
    %
    % Theoretical/practical motivation:
    %   - Direct conv() costs O(N*M) and becomes impractical.
    %   - With FFT: cost ~ O(NFFT log NFFT), much more scalable.
    %
    % Critical note:
    %   Sigma = 1/(2*pi*f*RelBandWidth)
    %   is the corrected choice that avoids huge kernels (the original script had 2*pi in the
    %   numerator, which made sigma enormous and WL very long, especially at 0.1 Hz).

    X = fft(x_, NFFT);
    % FFT the signal once.
    % Major optimization: inside the frequency loop only the kernel FFT is recomputed,
    % then multiplied pointwise with X.

    for k = 1:Nf
        fk = F(k);

        Sigma   = 1/(2*pi*fk*RelBandWidth);
        halfLen = ceil(SD_TIMES * Sigma / Dt);
        xx      = (-halfLen:halfLen) * Dt;
        % Sigma for fk and temporal axis of the kernel.
        % halfLen: samples to left/right of center (±4 sigma).

        WL = (1/sqrt(Sigma*sqrt(pi))) .* exp(-(xx.^2)/(2*Sigma^2)) .* exp(1i*2*pi*fk*xx);
        if useSingle, WL = single(WL); end
        % Build complex wavelet:
        %   - Gaussian envelope: localizes the filter in time
        %   - complex sinusoid: selects frequency fk (band center)
        %   - normalization: keeps comparability (approximately) across scales
        %
        % Theoretical motivation:
        %   - Gabor/Morlet-like wavelet: not the classical Morlet with DC correction,
        %     but often sufficient as a band-pass filter for 0.1–5 Hz LFP.
        %   - RelBandWidth controls relative spectral width (approximately constant Q).

        Y = fft(WL(:), NFFT);
        z = ifft(X .* Y) * Dt;
        p = abs(z).^2;
        % FFT convolution:
        %   - FFT(x) * FFT(WL) -> IFFT -> linear convolution (via implicit zero-padding in NFFT).
        % Power:
        %   - abs(z).^2 is instantaneous power of the band-pass filtered signal around fk.

        % Trim to length N (match your padding removal)
        p = p(halfLen + (1:N));
        % FFT convolution output is longer than N.
        % Remove leading padding samples (related to kernel length)
        % to align with the original N samples.

        % Fast moving-average per window using cumsum
        if useSingle
            cs = cumsum([single(0); single(p)]);
        else
            cs = cumsum([0; p]);
        end
        wsum = cs(starts + WindowSize) - cs(starts);
        WSP(k,:) = (wsum / WindowSize).';
        % Replace a nested loop with a standard technique:
        %   - cumsum gives window sums in O(N) total.
        % Formula:
        %   sum(p(i:i+W-1)) = cs(i+W) - cs(i)
        % Divide by WindowSize for the mean.
        %
        % Motivation:
        %   - speed: removes a loop over j=1:WindowSize and temporary matrices.
        %   - stability: fewer allocations -> less RAM fragmentation.

    end
end
% ---- local function (paste at end of your script) ----
function Zs = gaussSmooth2D_noToolbox(Z, sigma)
    if sigma <= 0
        Zs = Z;
        return;
    end

    % kernel size ~ 6*sigma (odd)
    halfw = ceil(3*sigma);
    x = (-halfw:halfw);
    g = exp(-(x.^2)/(2*sigma^2));
    g = g / sum(g);                 % normalize (sum=1)

    % replicate padding (no padarray)
    Zpad = padReplicate(Z, halfw);

    % separable convolution: first columns then rows
    Ztmp = conv2(Zpad, g(:),  'same');
    Ztmp = conv2(Ztmp, g(:).', 'same');

    % crop back to original size
    Zs = Ztmp(halfw+1:end-halfw, halfw+1:end-halfw);
end

function A = padReplicate(X, p)
    [m,n] = size(X);
    A = zeros(m+2*p, n+2*p, class(X));

    % center
    A(p+1:p+m, p+1:p+n) = X;

    % top/bottom replicate
    A(1:p,     p+1:p+n) = repmat(X(1,:),  p, 1);
    A(p+m+1:end, p+1:p+n) = repmat(X(end,:), p, 1);

    % left/right replicate (including corners)
    A(:, 1:p)       = repmat(A(:, p+1), 1, p);
    A(:, p+n+1:end) = repmat(A(:, p+n), 1, p);
end
