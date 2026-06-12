% Script to train the AR model on PMc for the Monkey C dataset.
% Changing the transient period to 1000(Wa.), 4000 (An.) and 6250 (Aw.).

%%
clear all;
%%
load('path\to\C\PreProcData.mat')
load('path\to\C\MEAMUALFP.mat')
%%
AllCh =97:192;
Ch = AllCh;
% Ch = 1:96;
for c=1:numel(AllCh)
    MUA(c,:) = log(MEAMUA.values(AllCh(c),:))-DataSet.LogMUAshift(AllCh(c));
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
for ChTarget = Ch-96%1:length(Ch)
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
    

    InLinComb = MUA(Ch-96,ndxLP);
    InputTest = MUA(Ch-96,ndxTP);
    Offsets = SmoothingWindow;
    for k = 2:KernelLength
        Offsets = [Offsets Offsets(k-1)+SmoothingWindow];
        InLinComb = [InLinComb; MUA(Ch-96,ndxLP-Offsets(end))];
        InputTest = [InputTest; MUA(Ch-96, ndxTP-Offsets(end))];
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
    %%

end
%%    
save('KernelAndPerfPM_trans6250.mat', 'Kernel','PerformanceTraining','OutReconstructTraining','PerformanceTest','PvalueTest','OutReconstructTest','OutOrig','OutOrigTest');
%%

       
MEAMap1 = [0 113 115 116 117 118 119 120 185 0;
           111 112 114 123 122 121 152 154 153 184;
           109 110 125 124 126 151 150 155 183 186;
           108 99 98 127 147 146 149 156 187 182;
           107 100 128 97 145 143 148 157 188 181; 
           106 101 102 144 141 142 159 158 189 180;
           103 137 138 139 140 131 129 160 179 190;
           104 136 135 134 133 132 130 177 178 191;
           105 169 170 171 172 173 174 175 176 192;
           0 168 167 166 165 164 163 162 161 0]; 

       
       figure
        [~,ndxCh2Map] = sort(MEAMap1(:));
        ndxCh2Map = ndxCh2Map(5:end);
%         [~,ndxM] = max(abs(LC.kernel));%mi dà quale degli istanti di tempo è il più importante per ogni canale
        MeasureMap = zeros(size(MEAMap1));
        MeasureMap(ndxCh2Map(AllCh-96)) = PerformanceTest;
%         MeasureMap(ndxCh2Map(LC.channels)) = PerformanceTest;
        imagesc(MeasureMap);
        colorbar();
%         caxis([0.5 0.85]);
caxis([0 0.8])
%         colormap(BlueRedCM);
        hold on
%         Z = mean(LC.kernel(1:3,:));
%         for k = find(Z > std(LC.kernel(:)))
%             [r,c] = find(MEAMap==LC.channels(k));
%             plot(c,r,'ko','MarkerFaceColor','w','MarkerSize',12*Z(k)/max(Z));
%         end
        for channel=Ch
            [r,c] = find(MEAMap1==channel);
            text(c,r,num2str(channel),'HorizontalAlignment','center','VerticalAlignment','middle')
        end
        set(gca,'TickDir','out','Box','on','Layer','top')
        xlabel('X MEA')
        ylabel('Y MEA')