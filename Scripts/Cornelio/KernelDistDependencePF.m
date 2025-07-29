clear all;
%%
load('C:\Users\Sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Cornelio\KernelAndPerfPF_trans6250.mat')
%%
% ChPF = setdiff(1:96,[66,20,88]);
% ChPM = setdiff(97:192,[178,174,113,131,133]);
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
Channels=1:96;
ChannelSet = setdiff(Channels,[50,62]);
% Channels=97:192;
% ChannelSet = setdiff(97:192, [174, 131, 133, 113]);

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



MEAMap= MEAMapPF;

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
    MeasureMap(ndxCh2Map(ChannelSet-96)) = ChNorm;
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
    scorepre = (coeff(:,1:PCtrunc)'*squeeze(Kernel(:,ch,ChannelSet)))';
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