clear all;
%%
load('C:\Users\Sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Cornelio\KernelAndPerfPM_trans1000.mat')
%%
ChPF = setdiff(1:96,[66,20,88]);
ChPM = setdiff(97:192,[178,174,113,131,133]);
%%
NoT = size(Kernel,1); %time-steps
NoC = size(Kernel,2); %channels-pre
NoCT = size(Kernel,3); % channels-post
dt = 5/1000;

%%
Z = reshape(Kernel,NoT,NoC*NoCT);
[coeff,score,latent] = pca(Z');
%%
EV = cumsum(latent);
EV = EV/EV(end);
PCtrunc=find(EV>0.9,1);
%%
ChNorm = zeros(NoC,1);
for pre = 1:NoC
   scorepre = (coeff(:,1:PCtrunc)'*squeeze(Kernel(:,pre,:)))';
   ChNorm(pre) = mean(sqrt(sum(scorepre.^2,2)));
end
figure
stem(ChNorm)
%%
% Channels=1:96;
% ChannelSet = setdiff(1:96,[66,20,88]);
Channels=97:192;
ChannelSet = Channels;

MEAMapPF = [ 0 9 8 7 10 11 12 13 15 0;
            72 73 40 41 5 4 3 14 16 17;
            71 74 39 42 6 32 2 29 18 19;
            70 75 38 43 48 1 31 28 27 20;
            69 76 37 44 45 49 51 30 26 21;
            68 77 36 35 46 47 50 55 25 22;
            67 78 34 33 63 52 53 54 56 23;
            66 79 81 64 62 61 60 59 58 24;
            65 80 82 83 93 92 91 87 57 89;
            0 96 95 94 84 85 86 90 88 0];

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
% ChNormRange = [0 0.275];
% pre = 1;
% scorepre = (coeff(:,1:3)'*squeeze(Kernel(:,pre,:)))';
% ChNorm = sqrt(sum(scorepre.^2,2));

figure
[~,ndxCh2Map] = sort(MEAMap(:));
ndxCh2Map = ndxCh2Map(5:end);
MeasureMap = zeros(size(MEAMap));
if MEAMap(2,2)>96
    MeasureMap(ndxCh2Map(ChannelSet-96)) = ChNorm;
else
    MeasureMap(ndxCh2Map(ChannelSet)) = ChNorm;
end
imagesc(MeasureMap);
colorbar();
hold on
for channel=ChannelSet
   [r,c] = find(MEAMap==channel);
   text(c,r,num2str(channel),'HorizontalAlignment','center','VerticalAlignment','middle','Color','k')
end
% [r,c] = find(MEAMap==pre);
% text(c,r,num2str(pre),'HorizontalAlignment','center','VerticalAlignment','middle','Color','r')

% caxis(ChNormRange)
set(gca,'TickDir','out','Box','on','Layer','top')
xlabel('X MEA')
ylabel('Y MEA')
%% quanto un pre contribuisce a ricostruire tutti i post
pre = 91;
scorepre = (coeff(:,1:PCtrunc)'*squeeze(Kernel(:,pre,:)))';
ChNorm = sqrt(sum(scorepre.^2,2));

figure
[~,ndxCh2Map] = sort(MEAMap(:));
ndxCh2Map = ndxCh2Map(5:end);
MeasureMap = zeros(size(MEAMap));
if MEAMap(2,2)>96
    MeasureMap(ndxCh2Map(Channels-96)) = ChNorm;
else
    MeasureMap(ndxCh2Map(Channels)) = ChNorm;
end
imagesc(MeasureMap);
colorbar();
hold on
for channel=setdiff(ChannelSet,pre)
   [r,c] = find(MEAMap==channel);
   text(c,r,num2str(channel),'HorizontalAlignment','center','VerticalAlignment','middle','Color','k')
end
[r,c] = find(MEAMap==pre);
text(c,r,num2str(pre),'HorizontalAlignment','center','VerticalAlignment','middle','Color','r')

% caxis(ChNormRange)
set(gca,'TickDir','out','Box','on','Layer','top')
xlabel('X MEA')
ylabel('Y MEA')
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
ChNorm=[];
for ch=1:NoC
    scorepre = (coeff(:,1:PCtrunc)'*squeeze(Kernel(:,ch,ChannelSet-96)))';
    ChNorm(:,ch) = sqrt(sum(scorepre.^2,2));
end
%%
meanNorm=[];
dist=unique(distance);

for i=1:numel(dist)
    ind = find(distance_all==dist(i));
    meanNorm(i) = mean(ChNorm(ind));
end
%%
figure
% scatter(distance_all,ChNorm,'c.')
% hold on
plot(dist,smooth(meanNorm),'r.-')