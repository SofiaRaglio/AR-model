% Script to test the percolation analysis on the original vs reconstructed LFP.
% This script reproduces figure 5c of the paper.
clear all;
%%
load('path\to\MEAMUALFP.mat')
load('path\to\PreProcData.mat')
%%
% AllCh=1:192;
% Channels = setdiff(1:192,[20,66,88,174,113,131,133]);
% 
% ChPF = setdiff(1:96,[66,20,88]);
% ChPM = setdiff(97:192,[174,113,131,133]);
%%
AllCh=1:192;
Channels = setdiff(1:192,[50,62]);

ChPF = setdiff(1:96,[50,62]);
ChPM = 97:192;

%%
for c=AllCh
    MUA(c,:) = log(MEAMUA.values(c,:))-DataSet.LogMUAshift(c);
end
LFP = MEALFP.values(AllCh,:);
%%
dt = MEAMUA.dt;
time = MEAMUA.time;
% TransientPeriod= 100;%3000,7600
TransientPeriod= 1000;%4000,6250
LearningPeriod  = 200;
TestPeriod      = 50; 
KernelLength = 69;
Life = TransientPeriod + LearningPeriod + TestPeriod;
t0 = time(1);
TSpan = TransientPeriod + LearningPeriod;
TSpan = time(time>=t0 & time<=t0+TSpan);
t = TSpan;

ndxLP = find(t-TSpan(1)>=+TransientPeriod & t-TSpan(1)<=LearningPeriod+TransientPeriod);
tLP = time(ndxLP);

t1 = TSpan(end);
TSpanTest = TransientPeriod + LearningPeriod + TestPeriod;
TSpanTest = time(time>=t0 & time<=t0+TSpanTest);
tt = TSpanTest;

ndxTP = find(tt-TSpan(1)>=TransientPeriod+LearningPeriod & tt-TSpan(1)<=LearningPeriod+TransientPeriod+TestPeriod);
tTP = tt(ndxTP);

OutOrig_LFP = LFP(Channels,ndxTP+1);
%%
% LFP_tot = LFP(Channels,ndx+1);
Matrix = corr(OutOrig_LFP',OutOrig_LFP');
[perc_threshold,Perc_Matrix_Pre, Sort_value_Pre,n_com_size_Pre,perc_threshold_step] = Percolation(Matrix);
%%
Perc_Matrix_pruned_Pre =[];
for i=1:size(Perc_Matrix_Pre,1)
    for j=1:size(Perc_Matrix_Pre,2)
        if Perc_Matrix_Pre(i,j)<0.1%0.5
           Perc_Matrix_pruned_Pre(i,j)=0;
        else
           Perc_Matrix_pruned_Pre(i,j)=Perc_Matrix_Pre(i,j);
        end
    end
end
%%
d  = graph(Perc_Matrix_pruned_Pre,'upper');
figure;
h = plot(d,'LineStyle','--','MarkerSize',5, 'NodeColor', 'g');
highlight(h, 1:93, 'NodeColor', 'r');
%%
figure;
imagesc(Perc_Matrix_pruned_Pre)
%%
load('path\to\KernelAndPerfPM2PF_trans1000.mat')
OutReconstructTest_PF = OutReconstructTest;

load('path\to\KernelAndPerfPF2PM_trans1000.mat')
OutReconstructTest_PM = OutReconstructTest;
%%
OutTest_LFP = [OutReconstructTest_PF(ChPF,:);OutReconstructTest_PM(ChPM-96,:)];
%%
Matrix_Post = corr(OutTest_LFP',OutTest_LFP');
[perc_threshold,Perc_Matrix_Post, Sort_value_Post,n_com_size_Post,perc_threshold_step] = Percolation(Matrix_Post);
%%
Perc_Matrix_pruned_Post =[];
for i=1:size(Perc_Matrix_Post,1)
    for j=1:size(Perc_Matrix_Post,2)
        if Perc_Matrix_Post(i,j)<0.1%0.5
           Perc_Matrix_pruned_Post(i,j)=0;
        else
           Perc_Matrix_pruned_Post(i,j)=Perc_Matrix_Post(i,j);
        end
    end
end
%%
d  = graph(Perc_Matrix_pruned_Post,'upper');
figure;
h = plot(d,'LineStyle','--','MarkerSize',5, 'NodeColor', 'g');
highlight(h, 1:93, 'NodeColor', 'r');
%%
figure;
imagesc(Perc_Matrix_pruned_Post)