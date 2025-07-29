clear all; clc;
%%
load('C:\Users\sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Preprocessing\Cornelio\PreProcData.mat')
load('C:\Users\sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Data\Cornelio\2Arrays\MEAMUALFP.mat')
%%
Options.SpecgramFreqs = [0.1 5];
Options.SpecgramFreqSamples = 500; % Num. of frequency samples in the log scale.
Options.SpecgramMovWnd = 1;%120;  % Window in seconds for the moving average.
Options.SpecgramRange = 10.^[+1.0 4.0];
Options.SamplingPeriod = 10.;
Options.SpecgramAsBitmap = 0; % If 1 use pcolor to plot spectrograms.
Options.SpecgramLevels = 20; % Num. of contour levels to plot.
Options.LFP.VarRange = [0 150];
Options.BandWidth = 0.1;
%%
Channels=97:192;
t_start = MEAMUA.time(1000);%3000;
t_end = MEAMUA.time(2040000);%6000;
time_s = round(t_start - MEAMUA.time(1))/0.005;
time_e = round(t_end - MEAMUA.time(1))/0.005;

RawSig = MEALFP;
RawSig.value = mean(MEALFP.values(Channels,time_s:time_e));
RawSig.time = MEALFP.time(time_s:time_e);

FEdges = exp(linspace(log(Options.SpecgramFreqs(1)),log(Options.SpecgramFreqs(2)),Options.SpecgramFreqSamples+1));
%%
WF=RawSig;
F =  FEdges
RelBandWidth = Options.BandWidth
MovingAverageWindow = Options.SpecgramMovWnd
SD_TIMES = 4.0;
MA_OVERLAP = 0.5;
Dt = diff(WF.time(1:2));
if exist('MovingAverageWindow','var')
WindowSize = round(MovingAverageWindow/Dt);
StepSize = round(MA_OVERLAP*MovingAverageWindow/Dt);
SampleNum = floor((length(WF.time) - WindowSize)/StepSize) + 1;
t = (0:SampleNum-1) * Dt * StepSize + mean(WF.time(1:WindowSize));
MAData = zeros(WindowSize, SampleNum);
else
t = WF.time;
end
WSP = zeros(length(F), length(t));
for k = 1:length(F)
disp(num2str(F(k),3));
% Makes the convolution kernel (the wavelet)...
Sigma = 2*pi/(F(k)*RelBandWidth);
x = (-ceil(SD_TIMES*Sigma/Dt):floor(SD_TIMES*Sigma/Dt))*Dt;
WL = 1/sqrt(Sigma*sqrt(pi)) * exp(-x.^2/(2*Sigma^2)) .* exp(1i*2*pi*F(k)*x);
% Works out the convolution (using FFT it's faster) and compute the power...
% Method 1: direct convolution, too slow...
z = conv(WF.value,WL) * Dt;
z = abs(z).^2;
%    % Method 2: convolution using FFT, fast but extensive use of memory...
%    X = fft([WF.value' zeros(1,length(WL)-1)]);
%    Y = fft([WL zeros(1,length(WF.value)-1)]);
%    z = abs(ifft(X.*Y) * Dt).^2;
% Method 3: convolution using FFT and dividing vectors in bunches, the
% best...
%    z = ConvPerBunches(WF.value,WL) * Dt;
% %    z = ConvPerBunches(WF.value, WL, 2^19) * Dt;
%    z = abs(z).^2;
% Removes padded elements at the beginning and the end...
z = z(ceil(SD_TIMES*Sigma/Dt)+(1:length(WF.value)));
%    z = z(ceil(length(WL)/2)+(1:length(WF.value)));
% Subsampling by moving average if required...
if exist('MovingAverageWindow','var')
for j = 1:WindowSize
MAdata(j,:) = z(j:StepSize:end-WindowSize+j);
end
WSP(k,:) = mean(MAdata);
else
WSP(k,:) = z;
end
end
Specgram.t = t;
Specgram.F = F;
Specgram.P = WSP;
%%
figure;
contourf(Specgram.t,Specgram.F,Specgram.P,':')

c = prctile(Specgram.P,[5 95]);
caxis([5e-7 5e-3])
xlabel('Time')
ylabel('Frequency')