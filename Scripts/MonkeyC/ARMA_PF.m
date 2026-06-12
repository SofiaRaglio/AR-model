% Script to train the AR model for the Monkey C PFc.
% Changing the transient period to 1000 (Wa.), 4000 (An.),6250 (Aw.).
clear all;
%%
load('path\to\C\PreProcData.mat')
load('path\to\C\MEAMUALFP.mat')
%%
AllCh=1:96;
Ch = setdiff(AllCh,[50,62]);
for c=1:96
    MUA(c,:) = log(MEAMUA.values(c,:))-DataSet.LogMUAshift(c);
end
LFP = MEALFP.values(AllCh,:);
dt = MEAMUA.dt;
time = MEAMUA.time;
TransientPeriod = 6250;
LearningPeriod  = 200;
TestPeriod      = 50; 
KernelLength = 69;
Life = TransientPeriod + LearningPeriod + TestPeriod;
Kernel = [];
PerformanceTraining = [];
OutReconstructTraining = [];
PvalueTrainig = [];

PerformanceTest = [];
OutReconstructTest = [];
PvalueTest = [];
%%
for ChTarget = Ch%1:length(Ch)
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
    %%
    OutOrig = LFP(ChTarget,ndxLP+1);
    OutOrigTest = LFP(ChTarget,ndxTP+1);
    SmoothingWindow = 1;
    

    InLinComb = MUA(Ch,ndxLP);
    InputTest = MUA(Ch,ndxTP);
    Offsets = SmoothingWindow;
    for k = 2:KernelLength
        Offsets = [Offsets Offsets(k-1)+SmoothingWindow];
        InLinComb = [InLinComb; MUA(Ch,ndxLP-Offsets(end))];
        InputTest = [InputTest; MUA(Ch, ndxTP-Offsets(end))];
    end
    
    Alpha = OutOrig*pinv(InLinComb);
    OutLinComb = Alpha*InLinComb;
    OutTestReconstruct = Alpha*InputTest;
    
    [R, P] = corrcoef(OutOrig,OutLinComb); % Correlation MUA orig-lin.comb.reconst.
    SimuData.RhoLinComb = R(1,2);
    SimuData.Significant = P(1,2);
    
    [R, P] = corrcoef(OutOrigTest,OutTestReconstruct); % Correlation MUA orig-lin.comb.reconst.
    SimuData.RhoTest = R(1,2);
    SimuData.SignTest = P(1,2);

    LC.kernel = reshape(Alpha,numel(Ch),KernelLength)';
    LC.time = Offsets*MEAMUA.dt;
    LC.channels = Ch;
    
    Kernel(:,:,ChTarget) = LC.kernel;
    
    PerformanceTraining(ChTarget) = SimuData.RhoLinComb;
    PvalueTraining(ChTarget) = SimuData.Significant;
    OutReconstructTraining(ChTarget,:) = OutLinComb;
    
    PerformanceTest(ChTarget) = SimuData.RhoTest;
    PvalueTest(ChTarget) = SimuData.SignTest;
    OutReconstructTest(ChTarget,:) = OutTestReconstruct;
 
end
%%    
save('KernelAndPerfPF_trans6250.mat', 'Kernel','PerformanceTraining','OutReconstructTraining','PerformanceTest','PvalueTest','OutReconstructTest','OutOrig','OutOrigTest');
%%
 MEAMap1 = [ 0 17 19 20 21 22 23 24 89 0;
            15 16 18 27 26 25 56 58 57 88;
            13 14 29 28 30 55 54 59 87 90;
            12 3 2 31 51 50 53 60 91 86;
            11 4 32 1 49 47 52 61 92 85;
            10 5 6 48 45 46 63 62 93 84;
            7 41 42 43 44 35 33 64 83 94;
            8 40 39 38 37 36 34 81 82 95;
            9 73 74 75 76 77 78 79 80 96;
            0 72 71 70 69 68 67 66 65 0]; 
       
     
       
       figure
        [~,ndxCh2Map] = sort(MEAMap1(:));
        ndxCh2Map = ndxCh2Map(5:end);
        MeasureMap = zeros(size(MEAMap1));
        MeasureMap(ndxCh2Map(AllCh)) = PerformanceTest;
        imagesc(MeasureMap);
        colorbar();
        caxis([0 0.8])
        hold on
        for channel=Ch
            [r,c] = find(MEAMap1==channel);
            text(c,r,num2str(channel),'HorizontalAlignment','center','VerticalAlignment','middle')
        end
        set(gca,'TickDir','out','Box','on','Layer','top')
        xlabel('X MEA')
        ylabel('Y MEA')