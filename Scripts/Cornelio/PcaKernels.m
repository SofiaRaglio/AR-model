clear all;
%% load the kernel with transient=100s
load('C:\Users\sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Cornelio\KernelAndPerfPM_trans1000.mat');
KernelAw = Kernel;
%% load the kernel with transient=3000s
load('C:\Users\sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Cornelio\KernelAndPerfPM_trans4000.mat');
KernelKet = Kernel;
%%
load('C:\Users\sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Cornelio\KernelAndPerfPM_trans6250.mat');
KernelAwing = Kernel;
%% select good channels for prefrontal cortex
% Good_Channels = setdiff(1:96,[66,20,88]);
% Good_Channels = setdiff(97:192,[174,113,131,133])-96;
%% reshaping the kernels in the correct way
% for Ch=1:96
%     Kernel_aw(Ch,:) = reshape(KernelAw(:,:,Ch),[1 96*69]);
%     Kernel_ket(Ch,:) = reshape(KernelKet(:,:,Ch),[1 96*69]);
% end
% %% putting the kernel awake and SO together and performing global PCA
% Kernel = [Kernel_aw Kernel_ket];
% [coef, score, latent] = pca(Kernel');
%%
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
%%
Z = [Z2 Z1];
[coeff,score,latent] = pca(Z');
%%
figure
plot3(score(1:9216,1),score(1:9216,2),score(1:9216,3),'.')
hold on
plot3(score(9217:end,1),score(9217:end,2),score(9217:end,3),'.')
set(gca,'TickDir','out')
grid on
xlabel('PC_1')
ylabel('PC_2')
zlabel('PC_3')
%% explained variance for global PCA
Var = cumsum(latent);
ExpVar = Var/Var(end);
%% superimposed kernel awake and SO
% figure
% scatter3(score(1:size(Kernel_aw,2),1),score(1:size(Kernel_aw,2),2),score(1:size(Kernel_aw,2),3))
% hold on
% scatter3(score(size(Kernel_aw,2)+1:end,1),score(size(Kernel_aw,2)+1:end,2),score(size(Kernel_aw,2)+1:end,3))
% legend('awake','SO')
%% reshaping kernel SO
KK3D = reshape(Kernel_ket,96,69,96); 
%% reshaping again kernel SO
KK3D = KK3D(Good_Channels,:,Good_Channels);
KKtime = KK3D(:,:,1);
for k = 2:numel(Good_Channels)
    KKtime = [KKtime; KK3D(:,:,k)];
end
%% performing PCA on kernel SO (riduco il chTarget per vedere come gli altri lo ricostruiscono)
[coef, score, latent] = pca(KKtime);
Var = cumsum(latent);
ExpVar = Var/Var(end);
%% explained variance and scatter of kernel SO
figure
subplot(1,2,1)
plot(ExpVar,'-d')
set(gca,'XScale','log')

subplot(1,2,2)
scatter3(score(:,1),score(:,2),score(:,3),'.')
xlabel('PC_1')
ylabel('PC_2')
zlabel('PC_3')
%% reshaping kernel awake
KA3D = reshape(Kernel_aw,96, 69, 96);%chTarget, time, reconstructing channels
%% reshaping again kernel awake
KA3D = KA3D(Good_Channels,:,Good_Channels);
KAtime = KA3D(:,:,1);
for k = 2:numel(Good_Channels)
    KAtime = [KAtime; KA3D(:,:,k)];
end
%% performing PCA on kernel awake
[coef, score, latent] = pca(KAtime);
Var = cumsum(latent);
ExpVar = Var/Var(end);
%% explained variance and scatter of kernel awake
figure
subplot(1,2,1)
plot(ExpVar,'-d')
set(gca,'XScale','log')

subplot(1,2,2)
scatter3(score(:,1),score(:,2),score(:,3),'.')
xlabel('PC_1')
ylabel('PC_2')
zlabel('PC_3')
%%
Ktime = [KKtime; KAtime];
%%
[coef, score, latent] = pca(Ktime);
Var = cumsum(latent);
ExpVar = Var/Var(end);
%%
figure
subplot(1,2,1)
plot(ExpVar,'-d')
set(gca,'XScale','log')

subplot(1,2,2)
scatter3(score(1:size(KKtime,1),1),score(1:size(KKtime,1),2),score(1:size(KKtime,1),3),'.')
hold on
scatter3(score(size(KKtime,1)+1:end,1),score(size(KKtime,1)+1:end,2),score(size(KKtime,1)+1:end,3),'.')
xlabel('PC_1')
ylabel('PC_2')
zlabel('PC_3')
hold on
for i=1:size(KK3D,1)
    ndx = i+(i-1)*size(KK3D,1);
    scatter3(score(ndx,1),score(ndx,2),score(ndx,3),'ro')
    scatter3(score(ndx+size(KKtime,1),1),score(ndx+size(KKtime,1),2),score(ndx+size(KKtime,1),3),'bo')
end
%% metto dei pallini dove ci sono io che ricostruisco me stesso
figure
scatter(score(1:size(KKtime,1),1),score(1:size(KKtime,1),2),'.')
hold on
scatter(score(size(KKtime,1)+1:end,1),score(size(KKtime,1)+1:end,2),'.')
xlabel('PC_1')
ylabel('PC_2')
% hold on
% for i=1:size(KK3D,1)
%     ndx = i+(i-1)*size(KK3D,1);
%     scatter(score(ndx,1),score(ndx,2),'ro')
%     scatter(score(ndx+size(KKtime,1),1),score(ndx+size(KKtime,1),2),'bo')
% end
%%
figure
scatter(score(1:size(KKtime,1),1),ones(size(KKtime,1),1),'.')
hold on
scatter(score(size(KKtime,1)+1:end,1),ones(size(KKtime,1),1),'.')
xlabel('PC_1')
hold on
ndx = [];
for i=1:size(KK3D,1)
    ndx = i+(i-1)*size(KK3D,1);
    scatter(score(ndx,1),1,'ro')
    scatter(score(ndx+size(KKtime,1),1),1,'bo')
end
    
    