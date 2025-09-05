clear all;
%%
load('MEAMUALFP.mat')
%%
TimeRange = [0 250] + 3000;
PFchannels = setdiff(1:96,[20,66,88]);
PMchannels = setdiff(97:192,[174,113,131,133]);

ndx = find(MEAMUA.time>=TimeRange(1) & MEAMUA.time<=TimeRange(2));

t = MEAMUA.time(ndx);
t = t - t(1);

ZPF = log(MEAMUA.values(PFchannels,ndx));
ZPM = log(MEAMUA.values(PMchannels,ndx));

mZPF = mean(ZPF,2);
ZPF = ZPF - repmat(mZPF,1,numel(t));
mZPM = mean(ZPM,2);
ZPM = ZPM - repmat(mZPM,1,numel(t));

[PF.coeff,PF.score,PF.latent] = pca(ZPF');
[PM.coeff,PM.score,PM.latent] = pca(ZPM');

LFP.ZPF = MEALFP.values(PFchannels,ndx);
LFP.ZPM = MEALFP.values(PMchannels,ndx);
%% Low-pass filtering of log(MUA).
lpCutOffFreq = 12.5;
for nc = 1:numel(PFchannels)
   lpZPF(nc,:) = lowpass(ZPF(nc,:),lpCutOffFreq,1/MEAMUA.dt);
   lp_LFP_PF(nc,:) = lowpass(LFP.ZPF(nc,:),lpCutOffFreq,1/MEAMUA.dt);
end
for nc = 1:numel(PMchannels)
   lpZPM(nc,:) = lowpass(ZPM(nc,:),lpCutOffFreq,1/MEAMUA.dt);
   lp_LFP_PM(nc,:) = lowpass(LFP.ZPM(nc,:),lpCutOffFreq,1/MEAMUA.dt);
end

%%
[PF.coeff,PF.score,PF.latent] = pca(lpZPF');
[PM.coeff,PM.score,PM.latent] = pca(lpZPM');

%%
ndx = find(MEALFP.time>=TimeRange(1) & MEALFP.time<=TimeRange(2));

LFP.t = MEALFP.time(ndx);
LFP.t = LFP.t - LFP.t(1);

LFP.ZPF = MEALFP.values(PFchannels,ndx);
LFP.ZPM = MEALFP.values(PMchannels,ndx);

% plot LFP raster
figure, imagesc(LFP.t,1:numel(PFchannels),LFP.ZPF)
CM = gradedColormap([0 0 1],[1 0 0]);
colormap(CM)
colorbar
caxis([-1 1]*0.75)
title('LFP PF cortex')

%% PF
[~, score, ~] = pca([lpZPF; LFP.ZPF]');
% [~, score, ~] = pca([mean(lpZPF); mean(LFP.ZPF)]');
PF.UpDownDect = score(:,1);
PF.MeanLogMUA = mean(lpZPF);

% Normalize.
PF.UpDownDect = PF.UpDownDect/std(PF.UpDownDect);
PF.MeanLogMUA = PF.MeanLogMUA/std(PF.MeanLogMUA);
PF.MeanLFP = mean(LFP.ZPF)/std(mean(LFP.ZPF));

R = corrcoef(PF.UpDownDect,PF.MeanLogMUA);
PF.UpDownDect = PF.UpDownDect*sign(R(1,2));

figure
subplot(2,1,2)
patch([t(1) t t(end)],[0 PF.MeanLogMUA 0],[1 1 1]*0.8)
hold on
plot(t,PF.MeanLogMUA,'k','LineWidth',0.75)
plot(t,PF.MeanLFP,'b','LineWidth',0.75)
plot(t,PF.UpDownDect,'r','LineWidth',0.75)
hcb = colorbar();
set(hcb,'Visible','off')
xlim([85,95])
title('mean MUA and U/D detector, PF cortex')
xlabel('Time, t [s]')

% plot LFP
subplot(2,1,1)
imagesc(LFP.t,1:numel(PMchannels),LFP.ZPM)
CM = gradedColormap([0 0 1],[1 0 0]);
colormap(CM)
colorbar
caxis([-1 1]*0.5)
xlim([85,95])
title('LFP, PF cortex')

%% PM
[~, score, ~] = pca([lpZPM; LFP.ZPM]');
% [~, score, ~] = pca([mean(lpZPM); mean(LFP.ZPM)]');
PM.UpDownDect = score(:,1);
PM.MeanLogMUA = mean(lpZPM);
PM.MeanLFP = mean(LFP.ZPM)/std(mean(LFP.ZPM));

% Normalize.
PM.UpDownDect = PM.UpDownDect/std(PM.UpDownDect);
PM.MeanLogMUA = PM.MeanLogMUA/std(PM.MeanLogMUA);

R = corrcoef(PM.UpDownDect,PM.MeanLogMUA);
PM.UpDownDect = PM.UpDownDect*sign(R(1,2));

figure
subplot(2,1,2)
patch([t(1) t t(end)],[0 PM.MeanLogMUA 0],[1 1 1]*0.8)
hold on
plot(t,PM.MeanLogMUA,'k','LineWidth',0.75)
plot(t,PM.MeanLFP,'b','LineWidth',0.75)
plot(t,PM.UpDownDect,'r','LineWidth',0.75)
xlim([85,95])
title('mean MUA and U/D detector, PM cortex')
xlabel('Time, t [s]')

% plot LFP
subplot(2,1,1)
imagesc(LFP.t,1:numel(PMchannels),LFP.ZPM)
CM = gradedColormap([0 0 1],[1 0 0]);
colormap(CM)
colorbar
caxis([-1 1]*0.75)
xlim([85,95])
title('LFP, PM cortex')
%%
PM.Transitions = zeros(1,numel(t));
PF.Transitions = zeros(1,numel(t));
for i=1:numel(t)
    if PM.UpDownDect(i)>0.4
        PM.Transitions(i)= 1;
    else
        PM.Transitions(i)= 0;
    end

    if PF.UpDownDect(i)>0.4
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
Threshold = 0.3;
%%
PMTrans.duration = diff(PMTrans.ndx)*0.005;
ndx=find(PMTrans.duration>=Threshold);
Up = find(PMTrans.val(ndx)==1);
Down= find(PMTrans.val(ndx)==0);
PM_Corr_Trans = PMTrans.duration(ndx);
PM_UpDuration = PM_Corr_Trans(Up);
PM_DownDuration = PM_Corr_Trans(Down);
%%
PFTrans.duration = diff(PFTrans.ndx)*0.005;
ndx=find(PFTrans.duration>=Threshold);
Up = find(PFTrans.val(ndx)==1);
Down= find(PFTrans.val(ndx)==0);
PF_Corr_Trans = PFTrans.duration(ndx);
%%
PF_UpDuration = PF_Corr_Trans(Up);
PF_DownDuration = PF_Corr_Trans(Down);
%%
PM.MeanUpDur = mean(PM_UpDuration);
PM.StdUpDur = std(PM_UpDuration);
PM.MeanDownDur = mean(PM_DownDuration);
PM.StdDownDur = std(PM_DownDuration);
%%
PF.MeanUpDur = mean(PF_UpDuration);
PF.StdUpDur = std(PF_UpDuration);
PF.MeanDownDur = mean(PF_DownDuration);
PF.StdDownDur = std(PF_DownDuration);
%%
figure
errorbar([PF.MeanUpDur PM.MeanUpDur], [PF.StdUpDur/sqrt(numel(PF_UpDuration)) PM.StdUpDur/sqrt(numel(PM_UpDuration))],'o--')
xlim([0.5 2.5])
ylim([0.2 0.4])
xticks([1:2])
xticklabels(["PF",'PM'])
title("UpStateDuration")

%%
figure
errorbar([PF.MeanDownDur PM.MeanDownDur], [PF.StdDownDur/sqrt(numel(PF_DownDuration)) PM.StdDownDur/sqrt(numel(PM_DownDuration))],'o--')
xlim([0.5 2.5])
ylim([0.4 0.8])
xticks([1:2])
xticklabels(["PF",'PM'])
title("DownStateDuration")
%%
[rPM,tau]= xcorr(mean(lpZPM),'normalize');
[rPF,tau]= xcorr(mean(lpZPF),'normalize');
figure,plot(tau*MEAMUA.dt,rPM,'c',tau*MEAMUA.dt,rPF,'b')
grid on
%%
[rPM,tau]= xcorr(mean(lp_LFP_PM),'normalize');
[rPF,tau]= xcorr(mean(lp_LFP_PF),'normalize');
figure,plot(tau*MEAMUA.dt,rPM,'c',tau*MEAMUA.dt,rPF,'b')
grid on
%%
[rPM,tau]= xcorr(mean(LFP.ZPM),'normalize');
[rPF,tau]= xcorr(mean(LFP.ZPF),'normalize');
figure,plot(tau*MEAMUA.dt,rPM,'c',tau*MEAMUA.dt,rPF,'b')
grid on