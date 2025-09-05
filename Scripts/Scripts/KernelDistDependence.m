clear all;
%%
% load('C:\Users\Windows\OneDrive - Università di Pavia\Desktop\Roba\Work\ARMAmodel\ARMAModel\Kernel\Cornelio\KernelAndPerfPF_trans1000.mat')
% KernelWa = Kernel;
% load('C:\Users\Windows\OneDrive - Università di Pavia\Desktop\Roba\Work\ARMAmodel\ARMAModel\Kernel\Cornelio\KernelAndPerfPF_trans4000.mat')
% KernelSO = Kernel;
% load('C:\Users\Windows\OneDrive - Università di Pavia\Desktop\Roba\Work\ARMAmodel\ARMAModel\Kernel\Cornelio\KernelAndPerfPF_trans6250.mat')
% KernelAw = Kernel;
%%
load('C:\Users\Windows\OneDrive - Università di Pavia\Desktop\Roba\Work\ARMAmodel\ARMAModel\Kernel\Cornelio\KernelAndPerfPM_trans1000.mat')
KernelWa = Kernel;
load('C:\Users\Windows\OneDrive - Università di Pavia\Desktop\Roba\Work\ARMAmodel\ARMAModel\Kernel\Cornelio\KernelAndPerfPM_trans4000.mat')
KernelSO = Kernel;
load('C:\Users\Windows\OneDrive - Università di Pavia\Desktop\Roba\Work\ARMAmodel\ARMAModel\Kernel\Cornelio\KernelAndPerfPM_trans6250.mat')
KernelAw = Kernel;
%% Bart
% ChPF = setdiff(1:96,[66,20,88]);
% ChPM = setdiff(97:192,[178,174,113,131,133]);
%% Cornelio
ChPF = setdiff(1:96,[50,62]);
ChPM = 97:192;
%%
NoT = size(Kernel,1); %time-steps
NoC = size(Kernel,2); %channels-pre
NoCT = size(Kernel,3); % channels-post
dt = 5/1000;

%%
% Channels=1:96;
% ChannelSet = setdiff(1:96,[50,62]);
Channels=97:192;
ChannelSet = 97:192;

% MEAMapPF = [ 0 9 8 7 10 11 12 13 15 0;
%             72 73 40 41 5 4 3 14 16 17;
%             71 74 39 42 6 32 2 29 18 19;
%             70 75 38 43 48 1 31 28 27 20;
%             69 76 37 44 45 49 51 30 26 21;
%             68 77 36 35 46 47 50 55 25 22;
%             67 78 34 33 63 52 53 54 56 23;
%             66 79 81 64 62 61 60 59 58 24;
%             65 80 82 83 93 92 91 87 57 89;
%             0 96 95 94 84 85 86 90 88 0];
% 
% MEAMapPM = [0 105 104 103 106 107 108 109 0 111;
%           168 169 136 137 101 100 99 110 112 113;
%           167 170 135 138 102 128 98 125 114 115;
%           166 171 134 139 144 97 127 124 123 116;
%           165 172 133 140 141 145 147 126 122 117;
%           164 173 132 131 142 143 146 151 121 118;
%           163 174 130 129 159 148 149 150 152 119;
%           162 175 177 160 158 157 156 155 154 120;
%           161 176 178 179 189 188 187 183 153 185;
%           0 192 191 190 180 181 182 186 184 0];
% 
% MEAMap= MEAMapPF;

 MEAMapPF = [ 0 17 19 20 21 22 23 24 89 0;
            15 16 18 27 26 25 56 58 57 88;
            13 14 29 28 30 55 54 59 87 90;
            12 3 2 31 51 50 53 60 91 86;
            11 4 32 1 49 47 52 61 92 85;
            10 5 6 48 45 46 63 62 93 84;
            7 41 42 43 44 35 33 64 83 94;
            8 40 39 38 37 36 34 81 82 95;
            9 73 74 75 76 77 78 79 80 96;
            0 72 71 70 69 68 67 66 65 0]; 
       

       
MEAMapPM = [0 113 115 116 117 118 119 120 185 0;
           111 112 114 123 122 121 152 154 153 184;
           109 110 125 124 126 151 150 155 183 186;
           108 99 98 127 147 146 149 156 187 182;
           107 100 128 97 145 143 148 157 188 181; 
           106 101 102 144 141 142 159 158 189 180;
           103 137 138 139 140 131 129 160 179 190;
           104 136 135 134 133 132 130 177 178 191;
           105 169 170 171 172 173 174 175 176 192;
           0 168 167 166 165 164 163 162 161 0]; 

MEAMap= MEAMapPM;

%% quanto ciascun pre contribuisce in media a ricostruire tutti i post

% figure
% [~,ndxCh2Map] = sort(MEAMap(:));
% ndxCh2Map = ndxCh2Map(5:end);
% MeasureMap = zeros(size(MEAMap));
% if MEAMap(2,2)>96
%     MeasureMap(ndxCh2Map(ChannelSet-96)) = ChNorm;
% else
%     MeasureMap(ndxCh2Map(ChannelSet)) = ChNorm;
% end
% imagesc(MeasureMap);
% colorbar();
% hold on
% for channel=ChannelSet
%    [r,c] = find(MEAMap==channel);
%    text(c,r,num2str(channel),'HorizontalAlignment','center','VerticalAlignment','middle','Color','k')
% end
% % [r,c] = find(MEAMap==pre);
% % text(c,r,num2str(pre),'HorizontalAlignment','center','VerticalAlignment','middle','Color','r')
% 
% % caxis(ChNormRange)
% set(gca,'TickDir','out','Box','on','Layer','top')
% xlabel('X MEA')
% ylabel('Y MEA')
%% quanto un pre contribuisce a ricostruire tutti i post
% pre = 1;
% scorepre = (coeff(:,1:PCtrunc)'*squeeze(Kernel(:,pre,:)))';
% ChNorm = sqrt(sum(scorepre.^2,2));
% 
% figure
% [~,ndxCh2Map] = sort(MEAMap(:));
% ndxCh2Map = ndxCh2Map(5:end);
% MeasureMap = zeros(size(MEAMap));
% if MEAMap(2,2)>96
%     MeasureMap(ndxCh2Map(Channels-96)) = ChNorm;
% else
%     MeasureMap(ndxCh2Map(Channels)) = ChNorm;
% end
% imagesc(MeasureMap);
% colorbar();
% hold on
% for channel=setdiff(ChannelSet,pre)
%    [r,c] = find(MEAMap==channel);
%    text(c,r,num2str(channel),'HorizontalAlignment','center','VerticalAlignment','middle','Color','k')
% end
% [r,c] = find(MEAMap==pre);
% text(c,r,num2str(pre),'HorizontalAlignment','center','VerticalAlignment','middle','Color','r')
% 
% % caxis(ChNormRange)
% set(gca,'TickDir','out','Box','on','Layer','top')
% xlabel('X MEA')
% ylabel('Y MEA')
%%
distance=[];
distance_all=[];
n=0;
for ch=ChannelSet
    n=n+1;
    m=0;
    [y,x]=find(MEAMap==ch);
    for channels=ChannelSet
        m=m+1;
        [y1,x1]=find(MEAMap==channels);
        distance(m)=sqrt((y1-y)^2+(x1-x)^2);
    end
    distance_all(n,:) = distance;
end
%%
Z = reshape(KernelWa,NoT,NoC*NoCT);
[coeff,score,latent] = pca(Z');
%%
EV = cumsum(latent);
EV = EV/EV(end);
PCtrunc=find(EV>0.9,1);
%%
ChNormWa=[];
for ch=1:NoC
    % scorepre = (coeff(:,1:PCtrunc)'*squeeze(KernelWa(:,ch,ChannelSet)))';
    scorepre = (coeff(:,1:PCtrunc)'*squeeze(KernelWa(:,ch,ChannelSet-96)))';
    ChNormWa(:,ch) = sqrt(sum(scorepre.^2,2));
end
%%
meanNorm=[];
dist=unique(distance);

for i=1:numel(dist)
    ind = find(distance_all==dist(i));
    meanNormWa(i) = mean(ChNormWa(ind));
    stdNormWa(i) = std(ChNormWa(ind))/sqrt(numel(ind));
end

%%
Z = reshape(KernelSO,NoT,NoC*NoCT);
[coeff,score,latent] = pca(Z');
%%
EV = cumsum(latent);
EV = EV/EV(end);
PCtrunc=find(EV>0.9,1);
%%
ChNormSO=[];
for ch=1:NoC
    % scorepre = (coeff(:,1:PCtrunc)'*squeeze(KernelSO(:,ch,ChannelSet)))';
    scorepre = (coeff(:,1:PCtrunc)'*squeeze(KernelSO(:,ch,ChannelSet-96)))';
    ChNormSO(:,ch) = sqrt(sum(scorepre.^2,2));
end
%%
meanNorm=[];
dist=unique(distance);

for i=1:numel(dist)
    ind = find(distance_all==dist(i));
    meanNormSO(i) = mean(ChNormSO(ind));
    stdNormSO(i) = std(ChNormSO(ind))/sqrt(numel(ind));
end

%%
Z = reshape(KernelAw,NoT,NoC*NoCT);
[coeff,score,latent] = pca(Z');
%%
EV = cumsum(latent);
EV = EV/EV(end);
PCtrunc=find(EV>0.9,1);
%%
ChNormAw=[];
for ch=1:NoC
    % scorepre = (coeff(:,1:PCtrunc)'*squeeze(KernelAw(:,ch,ChannelSet)))';
    scorepre = (coeff(:,1:PCtrunc)'*squeeze(KernelAw(:,ch,ChannelSet-96)))';
    ChNormAw(:,ch) = sqrt(sum(scorepre.^2,2));
end
%%
meanNorm=[];
dist=unique(distance);

for i=1:numel(dist)
    ind = find(distance_all==dist(i));
    meanNormAw(i) = mean(ChNormAw(ind));
    stdNormAw(i) = std(ChNormAw(ind))/sqrt(numel(ind));
end
%% For Camille
% figure
% curve1 = smooth(meanNormWa'+ (stdNormWa)')';
% curve2 = smooth(meanNormWa'- (stdNormWa)')';
% x2 = [dist, fliplr(dist)];
% inBetween = [curve1, fliplr(curve2)];
% fill(x2, inBetween, 'r','EdgeColor', 'none');
% hold on
% plot(dist, smooth(meanNormWa), '-', 'Color', 'r', 'LineWidth', 1., 'MarkerFaceColor', 'w');
% alpha(.2)
% 
% curve3 = smooth(meanNormSO'+ (stdNormSO)')';
% curve4 = smooth(meanNormSO'- (stdNormSO)')';
% x2 = [dist, fliplr(dist)];
% inBetween = [curve3, fliplr(curve4)];
% fill(x2, inBetween, 'b','EdgeColor', 'none');
% hold on
% plot(dist, smooth(meanNormSO), '-', 'Color', 'r', 'LineWidth', 1., 'MarkerFaceColor', 'w');
% alpha(.2)
% 
% curve5 = smooth(meanNormAw'+ (stdNormAw)')';
% curve6 = smooth(meanNormAw'- (stdNormAw)')';
% x2 = [dist, fliplr(dist)];
% inBetween = [curve5, fliplr(curve6)];
% fill(x2, inBetween, 'r','EdgeColor', 'none');
% hold on
% plot(dist, smooth(meanNormAw), '-', 'Color', 'g', 'LineWidth', 1., 'MarkerFaceColor', 'w');
% alpha(.2)
%% For Sofia (without toolboxes)
figure
curve1 = (meanNormWa'+ (stdNormWa)')';
curve2 = (meanNormWa'- (stdNormWa)')';
x2 = [dist, fliplr(dist)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'r','EdgeColor', 'none');
hold on
plot(dist, (meanNormWa), '-', 'Color', 'r', 'LineWidth', 1., 'MarkerFaceColor', 'w');
alpha(.2)

curve3 = (meanNormSO'+ (stdNormSO)')';
curve4 = (meanNormSO'- (stdNormSO)')';
x2 = [dist, fliplr(dist)];
inBetween = [curve3, fliplr(curve4)];
fill(x2, inBetween, 'b','EdgeColor', 'none');
hold on
plot(dist, (meanNormSO), '-', 'Color', 'b', 'LineWidth', 1., 'MarkerFaceColor', 'w');
alpha(.2)

curve5 = (meanNormAw'+ (stdNormAw)')';
curve6 = (meanNormAw'- (stdNormAw)')';
x2 = [dist, fliplr(dist)];
inBetween = [curve5, fliplr(curve6)];
fill(x2, inBetween, 'g','EdgeColor', 'none');
hold on
plot(dist, (meanNormAw), '-', 'Color', 'g', 'LineWidth', 1., 'MarkerFaceColor', 'w');
alpha(.2)