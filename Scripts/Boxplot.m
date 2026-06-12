%% Boxplot of the performance of the model
% This script is used to produce Figure 2b of the paper. 

clear all;
%%
addpath('path/to/Kernel/B/')
load('KernelAndPerfPF_trans100.mat')
x1 = PerformanceTest;
load('KernelAndPerfPF_trans3000.mat')
x2 = PerformanceTest;
load('KernelAndPerfPF_trans7600.mat')
x3 = PerformanceTest;
%%
addpath('path/to/Kernel/C/')
load('KernelAndPerfPF_trans1000.mat')
x4 = PerformanceTest;
load('KernelAndPerfPF_trans4000.mat')
x5 = PerformanceTest;
load('KernelAndPerfPF_trans6250.mat')
x6 = PerformanceTest;
%%
figure;
boxplot([x1',x2',x3'],'Notch','on','Labels',{'Wakefulness','Anesthesia','Awakening'})
hold on
boxplot([x4',x5',x6'],'Notch','on','Labels',{'Wakefulness','Anesthesia','Awakening'})
%%
addpath('path/to/Kernel/B/')
load('KernelAndPerfPM_trans100.mat')
x1 = PerformanceTest;
load('KernelAndPerfPM_trans3000.mat')
x2 = PerformanceTest;
load('KernelAndPerfPM_trans7600.mat')
x3 = PerformanceTest;
%%
addpath('path/to/Kernel/C/')
load('KernelAndPerfPM_trans1000.mat')
x4 = PerformanceTest;
load('KernelAndPerfPM_trans4000.mat')
x5 = PerformanceTest;
load('KernelAndPerfPM_trans6250.mat')
x6 = PerformanceTest;
%%
figure;
boxplot([x1',x2',x3'],'Notch','on','Labels',{'Wakefulness','Anesthesia','Awakening'})
hold on
boxplot([x4',x5',x6'],'Notch','on','Labels',{'Wakefulness','Anesthesia','Awakening'})
