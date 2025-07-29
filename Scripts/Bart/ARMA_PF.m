clear all;

%%
load('C:\Users\sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Data\Bart\MEAMUALFP.mat')
load('C:\Users\sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\NuoviDatiBart\Dati\PreProcessedData\PreProcData.mat')
% load('C:\Users\sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\NuoviDatiBart\Dati\MEAMUALFP.mat')
%%
% Ch =97:192;
AllCh=1:96;
Ch = setdiff(1:96,[66,20,88]);%1:96;
for c=1:96
    MUA(c,:) = log(MEAMUA.values(c,:))-DataSet.LogMUAshift(c);
end
LFP = MEALFP.values(AllCh,:);
dt = MEAMUA.dt;
time = MEAMUA.time;
TransientPeriod = 7600;
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
    %%

%         MEAMap = [ 0 41 39 37 43 45 47  1  5  0;
%             96 73 95 25 33 24 22  3  7  9;
%             94 75 93 27 35 16 20 10 11 13;
%             92 77 91 29 55 18 14  8  6 15;
%             90 79 89 31 49 57 61 12  4 17;
%             88 81 48 46 51 53 59 71  2 19;
%             86 83 44 42 38 63 65 67 69 21;
%             84 85 50 40 36 34 32 30 28 23;
%             82 87 52 54 74 72 70 62 26 66;
%             0 80 78 76 56 58 60 68 64  0];
%         
%         [~,ndxCh2Map] = sort(MEAMap(:));
%         ndxCh2Map = ndxCh2Map(5:end);
%         [~,ndxM] = max(abs(LC.kernel));%mi dà quale degli istanti di tempo è il più importante per ogni canale
%         MeasureMap = zeros(size(MEAMap));
%         MeasureMap(ndxCh2Map(LC.channels)) = LC.kernel(ndxM+KernelLength*(0:numel(LC.channels)-1));
%         imagesc(MeasureMap);
%         colorbar();
% %         colormap(BlueRedCM);
%         hold on
%         Z = mean(LC.kernel(1:3,:));
% %         for k = find(Z > std(LC.kernel(:)))
% %             [r,c] = find(MEAMap==LC.channels(k));
% %             plot(c,r,'ko','MarkerFaceColor','w','MarkerSize',12*Z(k)/max(Z));
% %         end
%         for channel=Ch
%             [r,c] = find(MEAMap==channel);
%             text(c,r,num2str(channel),'HorizontalAlignment','center','VerticalAlignment','middle')
%         end
%         set(gca,'TickDir','out','Box','on','Layer','top')
%         xlabel('X MEA')
%         ylabel('Y MEA')
end
%%    
save('KernelAndPerfPF_trans7600.mat', 'Kernel','PerformanceTraining','OutReconstructTraining','PerformanceTest','PvalueTest','OutReconstructTest','OutOrig','OutOrigTest');
%%
 MEAMap1 = [ 0 9 8 7 10 11 12 13 15 0;
            72 73 40 41 5 4 3 14 16 17;
            71 74 39 42 6 32 2 29 18 19;
            70 75 38 43 48 1 31 28 27 20;
            69 76 37 44 45 49 51 30 26 21;
            68 77 36 35 46 47 50 55 25 22;
            67 78 34 33 63 52 53 54 56 23;
            66 79 81 64 62 61 60 59 58 24;
            65 80 82 83 93 92 91 87 57 89;
            0 96 95 94 84 85 86 90 88 0]; 
       

% MEAMap1 = [ 0 9 8 7 10 11 12 13 15 0;
%            72 73 40 41 5 4 3 14 16 17;
%            71 74 39 42 6 32 2 29 18 19;
%            70 75 38 43 48 1 31 28 27 0;
%            69 76 37 44 45 49 51 30 26 21;
%            68 77 36 35 46 47 50 55 25 22;
%            67 78 34 33 63 52 53 54 56 23;
%            0 79 81 64 62 61 60 59 58 24;
%            65 80 82 83 93 92 91 87 57 89;
%            0 96 95 94 84 85 86 90 0 0];
       
MEAMap2 = [0 105 104 103 106 107 108 109 0 111;
           168 169 136 137 101 100 99 110 112 113;
           167 170 135 138 102 128 98 125 114 115;
           166 171 134 139 144 97 127 124 123 116;
           165 172 133 140 141 145 147 126 122 117; 
           164 173 132 131 142 143 146 151 121 118;
           163 174 130 129 159 148 149 150 152 119;
           162 175 177 160 158 157 156 155 154 120;
           161 176 178 179 189 188 187 183 153 185;
           0 192 191 190 180 181 182 186 184 0]; 
        
       
       figure
        [~,ndxCh2Map] = sort(MEAMap1(:));
        ndxCh2Map = ndxCh2Map(5:end);
%         ndxCh2Map = ndxCh2Map(8:end);
%         [~,ndxM] = max(abs(LC.kernel));%mi dà quale degli istanti di tempo è il più importante per ogni canale
        MeasureMap = zeros(size(MEAMap1));
%         MeasureMap(ndxCh2Map(LC.channels-96)) = PerformanceTest;
        MeasureMap(ndxCh2Map(AllCh)) = PerformanceTest;
        imagesc(MeasureMap);
        colorbar();
%         caxis([0.5 0.85])
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
        BlueRedCM = gradedColormap([7, 38, 78]/225,[175, 0, 0]/225);
colormap(BlueRedCM)