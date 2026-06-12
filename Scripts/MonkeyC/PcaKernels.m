% Script to compute the correlations between the Kernel PCA components for the Monkey C.
% This script is used to produce Suppl. Fig. 3a of the paper.
clear all;
%% load the kernel with transient=100s
load('path/to/KernelAndPerfPM_trans1000.mat');
KernelAw = Kernel;
%% load the kernel with transient=3000s
load('path/to/KernelAndPerfPM_trans4000.mat');
KernelKet = Kernel;
%%
load('path/to/KernelAndPerfPM_trans6250.mat');
KernelAwing = Kernel;
% ChPF = setdiff(1:96,[66,20,88]);
% ChPM = setdiff(97:192,[179,178,174,113,131,133]);
%%
NoT = size(KernelAw,1);
NoC = size(KernelAw,2);
NoCT = size(KernelAw,3);
dt = 5/1000;

%%
Z1 = reshape(KernelAw,NoT,NoC*NoCT);
%%
NoT = size(KernelKet,1);
NoC = size(KernelKet,2);
NoCT = size(KernelKet,3);

Z2 = reshape(KernelKet,NoT,NoC*NoCT);
%%
NoT = size(KernelAwing,1);
NoC = size(KernelAwing,2);
NoCT = size(KernelAwing,3);

Z3 = reshape(KernelAwing,NoT,NoC*NoCT);
%%
[cAw,s,lAw] = pca(Z1');

EV = cumsum(lAw);
EV = EV/EV(end);
nAw=find(EV>0.9,1);
%%
[cSO,s,lSO] = pca(Z2');

EV = cumsum(lSO);
EV = EV/EV(end);
nSO=find(EV>0.9,1);
%%
[cAwing,s,lAwing] = pca(Z3');

EV = cumsum(lAwing);
EV = EV/EV(end);
nAwing=find(EV>0.9,1);
%%
cAwth = cAw(:,1:nAw);
cSOth = cSO(:,1:nSO);
cAwingth = cAwing(:,1:nAwing);
%%
corr=[];
coef = [cAwth cSOth cAwingth];
for i=1:size(coef,2)
    for j=1:size(coef,2)
        corr(i,j) = dot(coef(:,i),coef(:,j));
    end
end
%%
figure
imagesc(abs(corr))
