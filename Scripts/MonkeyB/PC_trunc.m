% Script to compute the PCA truncation at 90% variance for the Monkey B.
% This script is used to produce Figure 4b-top of the paper.

clear all;

%%
ChPF = setdiff(1:96,[66,20,88]);
ChPM = setdiff(97:192,[178,174,113,131,133]);
%%
load('path/to/KernelAndPerfPF_trans1000.mat')
%%
NoT = size(Kernel,1); %time-steps
NoC = size(Kernel,2); %channels-pre
NoCT = size(Kernel,3); % channels-post
dt = 5/1000;

%%
Z = reshape(Kernel,NoT,NoC*NoCT);
[coeff,score,latent] = pca(Z');

%% truncation figure
EV = cumsum(latent);
EV = EV/EV(end);
PCtruncBWaPF=find(EV>0.9,1);
%%
load('path/to/KernelAndPerfPF_trans3000.mat')
%%
NoT = size(Kernel,1); %time-steps
NoC = size(Kernel,2); %channels-pre
NoCT = size(Kernel,3); % channels-post
dt = 5/1000;

%%
Z = reshape(Kernel,NoT,NoC*NoCT);
[coeff,score,latent] = pca(Z');

%% truncation figure
EV = cumsum(latent);
EV = EV/EV(end);
PCtruncBSOPF=find(EV>0.9,1);
%%
load('path/to/KernelAndPerfPF_trans7600.mat')
%%
NoT = size(Kernel,1); %time-steps
NoC = size(Kernel,2); %channels-pre
NoCT = size(Kernel,3); % channels-post
dt = 5/1000;

%%
Z = reshape(Kernel,NoT,NoC*NoCT);
[coeff,score,latent] = pca(Z');

%% truncation figure
EV = cumsum(latent);
EV = EV/EV(end);
PCtruncBAwPF=find(EV>0.9,1);
%%
load('path/to/KernelAndPerfPM_trans100.mat')
%%
NoT = size(Kernel,1); %time-steps
NoC = size(Kernel,2); %channels-pre
NoCT = size(Kernel,3); % channels-post
dt = 5/1000;

%%
Z = reshape(Kernel,NoT,NoC*NoCT);
[coeff,score,latent] = pca(Z');

%% truncation figure
EV = cumsum(latent);
EV = EV/EV(end);
PCtruncBWaPM=find(EV>0.9,1);
%%
load('path/to/KernelAndPerfPM_trans3000.mat')
%%
NoT = size(Kernel,1); %time-steps
NoC = size(Kernel,2); %channels-pre
NoCT = size(Kernel,3); % channels-post
dt = 5/1000;

%%
Z = reshape(Kernel,NoT,NoC*NoCT);
[coeff,score,latent] = pca(Z');

%% truncation figure
EV = cumsum(latent);
EV = EV/EV(end);
PCtruncBSOPM=find(EV>0.9,1);
%%
load('path/to/KernelAndPerfPM_trans7600.mat')

%%
NoT = size(Kernel,1); %time-steps
NoC = size(Kernel,2); %channels-pre
NoCT = size(Kernel,3); % channels-post
dt = 5/1000;

%%
Z = reshape(Kernel,NoT,NoC*NoCT);
[coeff,score,latent] = pca(Z');

%% truncation figure
EV = cumsum(latent);
EV = EV/EV(end);
PCtruncBAwPM=find(EV>0.9,1);


%%
Wa = mean([PCtruncBWaPF PCtruncBWaPM]);% PCtruncBWaPF2PM PCtruncBWaPM2PF]);
errWa = std([PCtruncBWaPF PCtruncBWaPM]);% PCtruncBWaPF2PM PCtruncBWaPM2PF]);
SO = mean([PCtruncBSOPF PCtruncBSOPM]);% PCtruncBSOPF2PM PCtruncBSOPM2PF]);
errSO = std([PCtruncBSOPF PCtruncBSOPM]);% PCtruncBSOPF2PM PCtruncBSOPM2PF]);
Aw = mean([PCtruncBAwPF PCtruncBAwPM]);% PCtruncBAwPF2PM PCtruncBAwPM2PF]);
errAw = std([PCtruncBAwPF PCtruncBAwPM]);% PCtruncBAwPF2PM PCtruncBAwPM2PF]);
%%
figure
errorbar([Wa SO Aw],[errWa errSO errAw],'.-','MarkerSize',15)
xlim([0.5 3.5])
xticks([1:3])
xticklabels({'Wakefulness','Slow Oscillations','Awakening'})
ylim([1 20])