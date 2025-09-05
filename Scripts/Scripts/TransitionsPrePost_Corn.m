clear all;
%%
% load('C:\Users\Windows\OneDrive - Università di Pavia\Desktop\Roba\Work\ARMAmodel\ARMAModel\Data\Bart\MEAMUALFP.mat')
% load('C:\Users\Windows\OneDrive - Università di Pavia\Desktop\Roba\Work\ARMAmodel\ARMAModel\Data\Bart\PreProcData.mat')

load('MEAMUALFP.mat')
load('PreProcData.mat')
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
% TransientPeriod = 3000;
TransientPeriod = 4000;
LearningPeriod  = 200;
TestPeriod      = 500; 
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

LFP_PF = LFP(ChPF,ndxTP+1);
LFP_PM = LFP(ChPM,ndxTP+1);
MUA_PF = MUA(ChPF,ndxTP+1);
MUA_PM = MUA(ChPM,ndxTP+1);
%%
clear MEAMUA;
clear MEALFP;
%%
SmoothingWindow = 1;
InputTest_PF = MUA(ChPF,ndxTP);
InputTest_PM = MUA(ChPM,ndxTP);
Offsets = SmoothingWindow;
for k = 2:KernelLength
    Offsets = [Offsets Offsets(k-1)+SmoothingWindow];
    InputTest_PF = [InputTest_PF; MUA(ChPF, ndxTP-Offsets(end))];
    InputTest_PM = [InputTest_PM; MUA(ChPM, ndxTP-Offsets(end))];
end
%%
% load('KernelAndPerfPF_trans4000.mat')
load('KernelAndPerfPM2PF_trans4000.mat')
% Out_PF = OutReconstructTest;
for i=1:96
    % Out_PF(i,:) = reshape(Kernel(:,:,i)',1,69*numel(ChPF))*InputTest_PF;
    Out_PF(i,:) = reshape(Kernel(:,:,i)',1,69*numel(ChPM))*InputTest_PM;
end

% load('KernelAndPerfPM_trans4000.mat')
load('KernelAndPerfPF2PM_trans4000.mat')
% Out_PM = OutReconstructTest;
for i=1:96
    % Out_PM(i,:) = reshape(Kernel(:,:,i)',1,69*numel(ChPM))*InputTest_PM;
    Out_PM(i,:) = reshape(Kernel(:,:,i)',1,69*numel(ChPF))*InputTest_PF;
end
%%
MEAMUA.dt = 0.005;
lpCutOffFreq = 12.5;
for nc = 1:numel(ChPF)
   lp_MUA_PF(nc,:) = lowpass(MUA_PF(nc,:),lpCutOffFreq,1/MEAMUA.dt);
   lp_LFP_PF(nc,:) = lowpass(LFP_PF(nc,:),lpCutOffFreq,1/MEAMUA.dt);
   lp_Out_PF(nc,:) = lowpass(Out_PF(nc,:),lpCutOffFreq,1/MEAMUA.dt);
end
for nc = 1:numel(ChPM)
   lp_MUA_PM(nc,:) = lowpass(MUA_PM(nc,:),lpCutOffFreq,1/MEAMUA.dt);
   lp_LFP_PM(nc,:) = lowpass(LFP_PM(nc,:),lpCutOffFreq,1/MEAMUA.dt);
   lp_Out_PM(nc,:) = lowpass(Out_PM(nc,:),lpCutOffFreq,1/MEAMUA.dt);
end

%%
PF.MeanLogMUA = zscore(mean(lp_MUA_PF));%mean(lp_MUA_PF)/std(mean(lp_MUA_PF));
PF.MeanOut = zscore(mean(lp_Out_PF));
PF.MeanLFP = zscore(mean(lp_LFP_PF));%mean(lp_LFP_PF)/std(mean(lp_LFP_PF));

PM.MeanLogMUA = zscore(mean(lp_MUA_PM));%mean(lp_MUA_PM)/std(mean(lp_MUA_PM));
PM.MeanLFP = zscore(mean(lp_LFP_PM));
PM.MeanOut = zscore(mean(lp_Out_PM));%mean(lp_LFP_PM)/std(mean(lp_LFP_PM));
%%
PM.Transitions = zeros(1,numel(ndxTP));
PF.Transitions = zeros(1,numel(ndxTP));
for i=1:numel(ndxTP)
    if PM.MeanLogMUA(i)>0.5
        PM.Transitions(i)= 1;
    else
        PM.Transitions(i)= 0;
    end

    if PF.MeanLogMUA(i)>0.5
        PF.Transitions(i)= 1;
    else
        PF.Transitions(i)= 0;
    end 
end
%%
PMTrans.ndx = find(diff(PM.Transitions) ~= 0)+1;
PMTrans.val = PM.Transitions(PMTrans.ndx); % 1 for upward and 0 for downward transitions.

PFTrans.ndx = find(diff(PF.Transitions) ~= 0)+1;
PFTrans.val = PF.Transitions(PFTrans.ndx); % 1 for upward and 0 for downward transitions.
%%
Threshold = 0.1;
%%
PMTrans.duration = diff(PMTrans.ndx)*0.005;
ndx=find(PMTrans.duration>=Threshold);
Up = find(PMTrans.val(ndx)==1);
Down= find(PMTrans.val(ndx)==0);
PM_Up = PMTrans.ndx(ndx(Up));
PM_Down = PMTrans.ndx(ndx(Down));
%%
PFTrans.duration = diff(PFTrans.ndx)*0.005;
ndx=find(PFTrans.duration>=Threshold);
Up = find(PFTrans.val(ndx)==1);
Down= find(PFTrans.val(ndx)==0);
PF_Up = PFTrans.ndx(ndx(Up));
PF_Down = PFTrans.ndx(ndx(Down));
%%
tWnd=[-0.4,0.4];
dt=0.005;
tVal = [fliplr(-dt:-dt:tWnd(1)) 0:dt:tWnd(2)];
ndxMask=round(tVal/dt);
%%
for i=1:numel(PF_Up)-2
    PF.LogMUA_d2u(i,:)=PF.MeanLogMUA(PF_Up(i+1)+ndxMask);
    PF.LFP_d2u(i,:)=PF.MeanLFP(PF_Up(i+1)+ndxMask);
    PF.Out_d2u(i,:)=PF.MeanOut(PF_Up(i+1)+ndxMask);
end

for i=2:numel(PF_Down)-2
    PF.LogMUA_u2d(i,:)=PF.MeanLogMUA(PF_Down(i+1)+ndxMask);
    PF.LFP_u2d(i,:)=PF.MeanLFP(PF_Down(i+1)+ndxMask);
    PF.Out_u2d(i,:)=PF.MeanOut(PF_Down(i+1)+ndxMask);
end
%%
for i=1:numel(PM_Up)-2
    PM.LogMUA_d2u(i,:)=PM.MeanLogMUA(PM_Up(i+1)+ndxMask);
    PM.LFP_d2u(i,:)=PM.MeanLFP(PM_Up(i+1)+ndxMask);
    PM.Out_d2u(i,:)=PM.MeanOut(PM_Up(i+1)+ndxMask);
end

for i=1:numel(PM_Down)-2
    PM.LogMUA_u2d(i,:)=PM.MeanLogMUA(PM_Down(i+1)+ndxMask);
    PM.LFP_u2d(i,:)=PM.MeanLFP(PM_Down(i+1)+ndxMask);
    PM.Out_u2d(i,:)=PM.MeanOut(PM_Down(i+1)+ndxMask);
end
%%
% mean_LogMUA_PF = mean(PF.LogMUA_u2d(:,41:80),2);
% [val_PF, ndx_PF] = sort(mean_LogMUA_PF);
% figure
% imagesc(tWnd,numel(PF_Down)-2,PF.LogMUA_u2d(ndx_PF,:));
% 
% figure
% imagesc(tWnd,numel(PF_Down)-2,PF.LFP_u2d(ndx_PF,:));
% % cm=gradedColormap([0 0 1],[1 0 0]);
% colormap(hot);
% 
% figure
% imagesc(tWnd,numel(PF_Down)-2,PF.Out_u2d(ndx_PF,:));
% % cm=gradedColormap([0 0 1],[1 0 0]);
% colormap(hot);
% 
% %%
% mean_LogMUA_PM = mean(PM.LogMUA_u2d(:,41:80),2);
% [val_PM, ndx_PM] = sort(mean_LogMUA_PM);
% figure
% imagesc(tWnd,numel(PM_Down)-2,PM.LogMUA_u2d(ndx_PM,:));
% 
% figure
% imagesc(tWnd,numel(PM_Down)-2,PM.LFP_u2d(ndx_PM,:));
% % cm=gradedColormap([0 0 1],[1 0 0]);
% % colormap(cm);
% 
% figure
% imagesc(tWnd,numel(PM_Down)-2,PM.Out_u2d(ndx_PM,:));
% % cm=gradedColormap([0 0 1],[1 0 0]);
% % colormap(cm);
% caxis([-3  3])
% %%
% thr_PF = quantile(val_PF,0.1);
% ind_PF = find(val_PF>=thr_PF,1);
% 
% thr_2_PF = quantile(val_PF,0.9);
% ind_2_PF = find(val_PF>=thr_2_PF,1);
% %%
% thr_PM = quantile(val_PM,0.1);
% ind_PM = find(val_PM>=thr_PM,1);
% 
% thr_2_PM = quantile(val_PM,0.9);
% ind_2_PM = find(val_PM>=thr_2_PM,1);
% %%
% figure;plot(mean(PF.Out_u2d(ndx_PF(1:ind_PF-1),:)))
% hold on
% plot(mean(PF.Out_u2d(ndx_PF(ind_2_PF:end),:)))
% yline(0)
% xline(80)
% ylim([-2 2])
% %%
% figure;plot(mean(PM.Out_u2d(ndx_PM(1:ind_PM-1),:)))
% hold on
% plot(mean(PM.Out_u2d(ndx_PM(ind_2_PM:end),:)))
% yline(0)
% xline(80)
% ylim([-2 2])
% %%
% figure;plot(mean(PF.LFP_u2d(ndx_PF(1:ind_PF-1),:)))
% hold on
% plot(mean(PF.LFP_u2d(ndx_PF(ind_2_PF:end),:)))
% yline(0)
% ylim([-2 2])
% xline(80)
% %%
% figure;plot(mean(PM.LFP_u2d(ndx_PM(1:ind_PM-1),:)))
% hold on
% plot(mean(PM.LFP_u2d(ndx_PM(ind_2_PM:end),:)))
% yline(0)
% ylim([-2 2])
% xline(80)
% %%
% figure;plot(mean(PF.LogMUA_u2d(ndx_PF(1:ind_PF-1),:)))
% hold on
% plot(mean(PF.LogMUA_u2d(ndx_PF(ind_2_PF:end),:)))
% yline(0)
% ylim([-2 2])
% xline(80)
% %%
% figure;plot(mean(PM.LogMUA_u2d(ndx_PM(1:ind_PM-1),:)))
% hold on
% plot(mean(PM.LogMUA_u2d(ndx_PM(ind_2_PM:end),:)))
% yline(0)
% ylim([-2 2])
% xline(80)
% %%
% for i=1:numel(PM_Down)-2
%     PFfromPM.LogMUA_u2d(i,:)=PF.MeanLogMUA(PM_Down(i+1)+ndxMask);
%     PFfromPM.LFP_u2d(i,:)=PF.MeanLFP(PM_Down(i+1)+ndxMask);
%     PFfromPM.Out_u2d(i,:)=PF.MeanOut(PM_Down(i+1)+ndxMask);
% end
% 
% figure
% imagesc(tWnd,numel(PM_Down)-2,PFfromPM.LogMUA_u2d(ndx_PM,:));
% caxis([-1  3])
% 
% figure
% imagesc(tWnd,numel(PM_Down)-2,PFfromPM.LFP_u2d(ndx_PM,:));
% cm=gradedColormap([0 0 1],[1 0 0]);
% colormap(cm);
% caxis([-3  3])
% 
% % figure
% % imagesc(tWnd,numel(PM_Down)-2,PFfromPM.Out_u2d(ndx_PM,:));
% % cm=gradedColormap([0 0 1],[1 0 0]);
% % colormap(cm);
% % caxis([-3  3])
% 
% %%
% for i=1:numel(PF_Down)-2
%     PMfromPF.LogMUA_u2d(i,:)=PM.MeanLogMUA(PF_Down(i+1)+ndxMask);
%     PMfromPF.LFP_u2d(i,:)=PM.MeanLFP(PF_Down(i+1)+ndxMask);
%     PMfromPF.Out_u2d(i,:)=PM.MeanOut(PF_Down(i+1)+ndxMask);
% end
% 
% figure
% imagesc(tWnd,numel(PF_Down)-2,PMfromPF.LogMUA_u2d(ndx_PF,:));
% caxis([-1  3])
% 
% figure
% imagesc(tWnd,numel(PF_Down)-2,PMfromPF.LFP_u2d(ndx_PF,:));
% cm=gradedColormap([0 0 1],[1 0 0]);
% colormap(cm);
% caxis([-3  3])
% 
% % figure
% % imagesc(tWnd,numel(PF_Down)-2,PMfromPF.Out_u2d(ndx_PF,:));
% % cm=gradedColormap([0 0 1],[1 0 0]);
% % colormap(cm);
% % caxis([-3  3])
% 
% %%
% figure;plot(ndxMask*0.005,mean(PFfromPM.LFP_u2d))
% hold on
% plot(ndxMask*0.005,mean(PMfromPF.LFP_u2d))
% plot(ndxMask*0.005,mean(PFfromPM.LFP_u2d)-mean(PMfromPF.LFP_u2d))
% %%
% figure;
% errorbar(mean(PM.LFP_u2d),std(PM.LFP_u2d))
% hold on
% errorbar(mean(PF.LFP_u2d),std(PF.LFP_u2d))
% %%
% meanNorm = mean(PM.LFP_u2d);
% stdNorm = std(PM.LFP_u2d);
% figure
% curve1 = smooth(meanNorm'+ (stdNorm)')';%smooth(meanNorm)'+ stdNorm;
% curve2 = smooth(meanNorm'- (stdNorm)')';%smooth(meanNorm)'- stdNorm;
% x2 = [ndxMask*0.005, fliplr(ndxMask*0.005)];
% inBetween = [curve1, fliplr(curve2)];
% fill(x2, inBetween, 'r','EdgeColor', 'none');
% hold on
% plot(ndxMask*0.005, smooth(meanNorm), '-', 'Color', 'r', 'LineWidth', 1., 'MarkerFaceColor', 'w');
% alpha(.2)
% hold on
% meanNorm = mean(PF.LFP_u2d);
% stdNorm = std(PF.LFP_u2d);
% 
% curve1 = smooth(meanNorm'+ (stdNorm)')';%smooth(meanNorm)'+ stdNorm;
% curve2 = smooth(meanNorm'- (stdNorm)')';%smooth(meanNorm)'- stdNorm;
% x2 = [ndxMask*0.005, fliplr(ndxMask*0.005)];
% inBetween = [curve1, fliplr(curve2)];
% fill(x2, inBetween, 'b','EdgeColor', 'none');
% hold on
% plot(ndxMask*0.005, smooth(meanNorm), '-', 'Color', 'b', 'LineWidth', 1., 'MarkerFaceColor', 'w');
% alpha(.2)



%% New fig
mean_LogMUA_PF = mean(PF.LogMUA_u2d(:,41:80),2);
[val_PF, ndx_PF] = sort(mean_LogMUA_PF);
figure
imagesc(tWnd,numel(PF_Down)-2,PF.LogMUA_u2d(ndx_PF,:));

figure
imagesc(tWnd,numel(PF_Down)-2,PF.LFP_u2d(ndx_PF,:));
% cm=gradedColormap([0 0 1],[1 0 0]);
colormap(hot);

figure
imagesc(tWnd,numel(PF_Down)-2,PF.Out_u2d(ndx_PF,:));
% cm=gradedColormap([0 0 1],[1 0 0]);
colormap(hot);
caxis([-3  3])
%%
a=round(size(PF.LogMUA_u2d,1)/2);
b=size(PF.LogMUA_u2d,1)-a;
figure;plot(mean(PF.LFP_u2d(ndx_PF(1:a),:)))
hold on
plot(mean(PF.LFP_u2d(ndx_PF(a+1:end),:)))
yline(0)
xline(80)
ylim([-2 2])
%%
a=round(size(PF.LogMUA_u2d,1)/2);
b=size(PF.LogMUA_u2d,1)-a;
figure;plot(mean(PF.Out_u2d(ndx_PF(1:a),:)))
hold on
plot(mean(PF.Out_u2d(ndx_PF(a+1:end),:)))
yline(0)
xline(80)
ylim([-2 2])
%%
mean_LogMUA_PM = mean(PM.LogMUA_u2d(:,41:80),2);
[val_PM, ndx_PM] = sort(mean_LogMUA_PM);
figure
imagesc(tWnd,numel(PM_Down)-2,PM.LogMUA_u2d(ndx_PM,:));

figure
imagesc(tWnd,numel(PM_Down)-2,PM.LFP_u2d(ndx_PM,:));
% cm=gradedColormap([0 0 1],[1 0 0]);
colormap(hot);

figure
imagesc(tWnd,numel(PM_Down)-2,PM.Out_u2d(ndx_PM,:));
% cm=gradedColormap([0 0 1],[1 0 0]);
colormap(hot);
caxis([-3  3])
%%
a=round(size(PM.LogMUA_u2d,1)/2);
b=size(PM.LogMUA_u2d,1)-a;
figure;plot(mean(PM.LFP_u2d(ndx_PM(1:a),:)))
hold on
plot(mean(PM.LFP_u2d(ndx_PM(a+1:end),:)))
yline(0)
xline(80)
ylim([-2 2])
%%
a=round(size(PM.LogMUA_u2d,1)/2);
b=size(PM.LogMUA_u2d,1)-a;
figure;plot(mean(PM.Out_u2d(ndx_PM(1:a),:)))
hold on
plot(mean(PM.Out_u2d(ndx_PM(a+1:end),:)))
yline(0)
xline(80)
ylim([-2 2])