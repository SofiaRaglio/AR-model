clear all;
%%
load('C:\Users\Sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Bart\KernelAndPerfPM_trans7600.mat')
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
ChannelSet = setdiff(97:192, [174, 131, 133, 113]);

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

MEAMapPM = [0 105 104 103 106 107 108 109 0 111;
          168 169 136 137 101 100 99 110 112 113;
          167 170 135 138 102 128 98 125 114 115;
          166 171 134 139 144 97 127 124 123 116;
          165 172 133 140 141 145 147 126 122 117;
          164 173 132 131 142 143 146 151 121 118;
          163 174 130 129 159 148 149 150 152 119;
          162 175 177 160 158 157 156 155 154 120;
          161 176 178 179 189 188 187 183 153 185;
          0 192 191 190 180 181 182 186 184 0];

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