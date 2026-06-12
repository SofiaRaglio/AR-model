%% Boxplot of the performance of the model
% This script is used to produce Figure 5a of the paper. 

clear all;
%%

addpath('path/to/Kernel/B/')
load('KernelAndPerfPM2PF_trans100.mat')
x1 = PerformanceTest;
load('KernelAndPerfPM2PF_trans3000.mat')
x2 = PerformanceTest;
load('KernelAndPerfPM2PF_trans7600.mat')
x3 = PerformanceTest;
%%
addpath('path/to/Kernel/B/')
load('KernelAndPerfPF2PM_trans100.mat')
x4 = PerformanceTest;
load('KernelAndPerfPF2PM_trans3000.mat')
x5 = PerformanceTest;
load('KernelAndPerfPF2PM_trans7600.mat')
x6 = PerformanceTest;

%%
figure;
boxplot([x1',x2',x3'],'Notch','on','Colors','k','Labels',{'Wakefulness','Anesthesia','Awakening'})
hold on
boxplot([x4',x5',x6'],'Notch','on','Colors','r','Labels',{'Wakefulness','Anesthesia','Awakening'})
%%
addpath('path/to/Kernel/C/')
load('KernelAndPerfPM2PF_trans1000.mat')
x1 = PerformanceTest;
load('KernelAndPerfPM2PF_trans4000.mat')
x2 = PerformanceTest;
load('KernelAndPerfPM2PF_trans6250.mat')
x3 = PerformanceTest;
%%
addpath('path/to/Kernel/C/')
load('KernelAndPerfPF2PM_trans1000.mat')
x4 = PerformanceTest;
load('KernelAndPerfPF2PM_trans4000.mat')
x5 = PerformanceTest;
load('KernelAndPerfPF2PM_trans6250.mat')
x6 = PerformanceTest;
%%
figure;
boxplot([x1',x2',x3'],'Notch','on','Colors','k','Labels',{'Wakefulness','Anesthesia','Awakening'})
hold on
boxplot([x4',x5',x6'],'Notch','on','Colors','r','Labels',{'Wakefulness','Anesthesia','Awakening'})
