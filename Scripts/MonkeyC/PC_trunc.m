% Script to compute the PCA truncation at 90% variance for the Monkey C.
% This script is used to produce Figure 2b of the paper.
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
PCtruncCWaPF=find(EV>0.9,1);
%%
load('path/to/KernelAndPerfPF_trans4000.mat')

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
PCtruncCSOPF=find(EV>0.9,1);
%%
load('path/to/KernelAndPerfPF_trans6250.mat')
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
PCtruncCAwPF=find(EV>0.9,1);
%%
load('path/to/KernelAndPerfPM_trans1000.mat')
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
PCtruncCWaPM=find(EV>0.9,1);
%%
load('path/to/KernelAndPerfPM_trans4000.mat')
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
PCtruncCSOPM=find(EV>0.9,1);
%%
load('path/to/KernelAndPerfPM_trans6250.mat')
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
PCtruncCAwPM=find(EV>0.9,1);
%%
load('path/to/KernelAndPerfPM2PF_trans1000.mat')
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
PCtruncCWaPM2PF=find(EV>0.9,1);
%%
load('path/to/KernelAndPerfPM2PF_trans4000.mat')

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
PCtruncCSOPM2PF=find(EV>0.9,1);
%%
load('path/to/KernelAndPerfPM2PF_trans6250.mat')
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
PCtruncCAwPM2PF=find(EV>0.9,1);
%%
load('path/to/KernelAndPerfPF2PM_trans1000.mat')
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
PCtruncCWaPF2PM=find(EV>0.9,1);
%%
load('path/to/KernelAndPerfPF2PM_trans4000.mat')
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
PCtruncCSOPF2PM=find(EV>0.9,1);
%%
load('path/to/KernelAndPerfPF2PM_trans6250.mat')
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
PCtruncCAwPF2PM=find(EV>0.9,1);

%%
Wa = mean([PCtruncCWaPF PCtruncCWaPM PCtruncCWaPF2PM PCtruncCWaPM2PF]);
errWa = std([PCtruncCWaPF PCtruncCWaPM PCtruncCWaPF2PM PCtruncCWaPM2PF]);
SO = mean([PCtruncCSOPF PCtruncCSOPM PCtruncCSOPF2PM PCtruncCSOPM2PF]);
errSO = std([PCtruncCSOPF PCtruncCSOPM PCtruncCSOPF2PM PCtruncCSOPM2PF]);
Aw = mean([PCtruncCAwPF PCtruncCAwPM PCtruncCAwPF2PM PCtruncCAwPM2PF]);
errAw = std([PCtruncCAwPF PCtruncCAwPM PCtruncCAwPF2PM PCtruncCAwPM2PF]);
%%
figure
errorbar([Wa SO Aw],[errWa errSO errAw],'.-','MarkerSize',15)
xlim([0.5 3.5])
xticks([1:3])
xticklabels({'Wakefulness','Slow Oscillations','Awakening'})
ylim([1 20])