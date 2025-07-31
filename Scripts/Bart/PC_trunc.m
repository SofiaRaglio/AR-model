clear all;

%%
ChPF = setdiff(1:96,[66,20,88]);
ChPM = setdiff(97:192,[178,174,113,131,133]);
%%
load('C:\Users\Sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Bart\KernelAndPerfPF_trans100.mat')
% load('C:\Users\Sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Cornelio\KernelAndPerfPF_trans1000.mat')
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
PCtruncBartWaPF=find(EV>0.9,1);
%%
load('C:\Users\Sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Bart\KernelAndPerfPF_trans3000.mat')
% load('C:\Users\Sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Cornelio\KernelAndPerfPF_trans4000.mat')
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
PCtruncBartSOPF=find(EV>0.9,1);
%%
load('C:\Users\Sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Bart\KernelAndPerfPF_trans7600.mat')
% load('C:\Users\Sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Cornelio\KernelAndPerfPF_trans6250.mat')
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
PCtruncBartAwPF=find(EV>0.9,1);
%%
load('C:\Users\Sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Bart\KernelAndPerfPM_trans100.mat')
% load('C:\Users\Sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Cornelio\KernelAndPerfPM_trans1000.mat')
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
PCtruncBartWaPM=find(EV>0.9,1);
%%
load('C:\Users\Sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Bart\KernelAndPerfPM_trans3000.mat')
% load('C:\Users\Sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Cornelio\KernelAndPerfPM_trans4000.mat')
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
PCtruncBartSOPM=find(EV>0.9,1);
%%
load('C:\Users\Sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Bart\KernelAndPerfPM_trans7600.mat')
% load('C:\Users\Sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Cornelio\KernelAndPerfPM_trans6250.mat')
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
PCtruncBartAwPM=find(EV>0.9,1);
%%
% load('C:\Users\Sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Bart\KernelAndPerfPM2PF_trans100.mat')
% %%
% NoT = size(Kernel,1); %time-steps
% NoC = size(Kernel,2); %channels-pre
% NoCT = size(Kernel,3); % channels-post
% dt = 5/1000;
% 
% %%
% Z = reshape(Kernel,NoT,NoC*NoCT);
% [coeff,score,latent] = pca(Z');
% 
% %% truncation figure
% EV = cumsum(latent);
% EV = EV/EV(end);
% PCtruncBartWaPM2PF=find(EV>0.9,1);
%%
% load('C:\Users\Sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Bart\KernelAndPerfPM2PF_trans3000.mat')
% 
% %%
% NoT = size(Kernel,1); %time-steps
% NoC = size(Kernel,2); %channels-pre
% NoCT = size(Kernel,3); % channels-post
% dt = 5/1000;
% 
% %%
% Z = reshape(Kernel,NoT,NoC*NoCT);
% [coeff,score,latent] = pca(Z');
% 
% %% truncation figure
% EV = cumsum(latent);
% EV = EV/EV(end);
% PCtruncBartSOPM2PF=find(EV>0.9,1);
% %%
% load('C:\Users\Sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Bart\KernelAndPerfPM2PF_trans7600.mat')
% %%
% NoT = size(Kernel,1); %time-steps
% NoC = size(Kernel,2); %channels-pre
% NoCT = size(Kernel,3); % channels-post
% dt = 5/1000;
% 
% %%
% Z = reshape(Kernel,NoT,NoC*NoCT);
% [coeff,score,latent] = pca(Z');
% 
% %% truncation figure
% EV = cumsum(latent);
% EV = EV/EV(end);
% PCtruncBartAwPM2PF=find(EV>0.9,1);
% %%
% load('C:\Users\Sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Bart\KernelAndPerfPF2PM_trans100.mat')
% %%
% NoT = size(Kernel,1); %time-steps
% NoC = size(Kernel,2); %channels-pre
% NoCT = size(Kernel,3); % channels-post
% dt = 5/1000;
% 
% %%
% Z = reshape(Kernel,NoT,NoC*NoCT);
% [coeff,score,latent] = pca(Z');
% 
% %% truncation figure
% EV = cumsum(latent);
% EV = EV/EV(end);
% PCtruncBartWaPF2PM=find(EV>0.9,1);
% %%
% load('C:\Users\Sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Bart\KernelAndPerfPF2PM_trans3000.mat')
% %%
% NoT = size(Kernel,1); %time-steps
% NoC = size(Kernel,2); %channels-pre
% NoCT = size(Kernel,3); % channels-post
% dt = 5/1000;
% 
% %%
% Z = reshape(Kernel,NoT,NoC*NoCT);
% [coeff,score,latent] = pca(Z');
% 
% %% truncation figure
% EV = cumsum(latent);
% EV = EV/EV(end);
% PCtruncBartSOPF2PM=find(EV>0.9,1);
% %%
% load('C:\Users\Sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Bart\KernelAndPerfPF2PM_trans7600.mat')
% %%
% NoT = size(Kernel,1); %time-steps
% NoC = size(Kernel,2); %channels-pre
% NoCT = size(Kernel,3); % channels-post
% dt = 5/1000;
% 
% %%
% Z = reshape(Kernel,NoT,NoC*NoCT);
% [coeff,score,latent] = pca(Z');
% 
% %% truncation figure
% EV = cumsum(latent);
% EV = EV/EV(end);
% PCtruncBartAwPF2PM=find(EV>0.9,1);

%%
Wa = mean([PCtruncBartWaPF PCtruncBartWaPM]);% PCtruncBartWaPF2PM PCtruncBartWaPM2PF]);
errWa = std([PCtruncBartWaPF PCtruncBartWaPM]);% PCtruncBartWaPF2PM PCtruncBartWaPM2PF]);
SO = mean([PCtruncBartSOPF PCtruncBartSOPM]);% PCtruncBartSOPF2PM PCtruncBartSOPM2PF]);
errSO = std([PCtruncBartSOPF PCtruncBartSOPM]);% PCtruncBartSOPF2PM PCtruncBartSOPM2PF]);
Aw = mean([PCtruncBartAwPF PCtruncBartAwPM]);% PCtruncBartAwPF2PM PCtruncBartAwPM2PF]);
errAw = std([PCtruncBartAwPF PCtruncBartAwPM]);% PCtruncBartAwPF2PM PCtruncBartAwPM2PF]);
%%
figure
errorbar([Wa SO Aw],[errWa errSO errAw],'.-','MarkerSize',15)
xlim([0.5 3.5])
xticks([1:3])
xticklabels({'Wakefulness','Slow Oscillations','Awakening'})
ylim([1 20])