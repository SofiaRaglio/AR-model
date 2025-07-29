clear all;
%%
load('C:\Users\Sofia\OneDrive - Istituto Superiore di Sanità\ARMAmodel\ARMAModel\Kernel\Bart\KernelAndPerfPM2PF_trans3000.mat')
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

%% truncation figure
EV = cumsum(latent);
EV = EV/EV(end);
figure
stem(EV*100)
set(gca,'XScale','log')
xlabel('PCs')
ylabel('Expl. var (%)')
thr=find(EV>0.9,1);
title(sprintf('PCs truncation=%s'),thr)
%% pcs shape -- come ogni tempo nel passato contribuisce a ciascuna PC
NoPC = 3;
figure
hold all
hlgn = zeros(1,NoPC);
for k = 1:NoPC
   hlgn(k) = plot(([1:69]*0.005),coeff(:,k));
   slgn{k} = ['PC_' num2str(k)];
end
legend(hlgn,slgn)
plot(([1:69]*0.005),[0 0],'k--');
set(gca,'TickDir','out')
% grid on

%% star
figure
plot3(score(:,1),score(:,2),score(:,3),'.')
set(gca,'TickDir','out')
grid on
xlabel('PC_1')
ylabel('PC_2')
zlabel('PC_3')

hold on
for i=1
    plot3(score(1+(i-1)*93+(i-1),1),score(1+(i-1)*93+(i-1),2),score(1+(i-1)*93+(i-1),3),'ro')
end

%%
figure
plot(score(:,1),score(:,2),'b.')
set(gca,'TickDir','out')
% grid on
xlabel('PC_1')
ylabel('PC_2')

hold on
% for i=1
%     plot(score(1+(i-1)*93+(i-1),1),score(1+(i-1)*93+(i-1),2),'ro')
% end
pre=1;
for i=1:96
    plot(score(pre+(i-1)*93,1),score(pre+(i-1)*93,2),'Marker','.','Color','m', 'MarkerSize',10)
end
plot(score(pre+((pre-1)*93),1),score(pre+((pre-1)*93),2),'ko')

pre=8;
for i=1:96
    plot(score(pre+(i-1)*93,1),score(pre+(i-1)*93,2),'Marker','.','Color','m', 'MarkerSize',10)
end
% plot(score(pre+((pre-1)*93),1),score(pre+((pre-1)*93),2),'ko')

pre=67;
for i=1:96
    plot(score(pre+(i-1)*93,1),score(pre+(i-1)*93,2),'Marker','.','Color','g', 'MarkerSize',10)
end
% plot(score(pre+((pre-1)*93),1),score(pre+((pre-1)*93),2),'ko')

% pre=88;
% for i=1:96
%     plot(score(pre+(i-1)*93,1),score(pre+(i-1)*93,2),'Marker','.','Color','c', 'MarkerSize',10)
% end
% plot(score(pre+((pre-1)*93),1),score(pre+((pre-1)*93),2),'ko')

pre=75;
for i=1:96
    plot(score(pre+(i-1)*93,1),score(pre+(i-1)*93,2),'Marker','.','Color','y', 'MarkerSize',10)
end
% plot(score(pre+((pre-1)*93),1),score(pre+((pre-1)*93),2),'ko')

pre=35;
for i=1:96
    plot(score(pre+(i-1)*93,1),score(pre+(i-1)*93,2),'Marker','.','Color','m', 'MarkerSize',10)
end
% plot(score(pre+((pre-1)*93),1),score(pre+((pre-1)*93),2),'ko')

% pre=36;
% for i=1:96
%     plot(score(pre+(i-1)*93,1),score(pre+(i-1)*93,2),'Marker','.','Color','c', 'MarkerSize',10)
% end
% plot(score(pre+((pre-1)*93),1),score(pre+((pre-1)*93),2),'ko')

pre=2;
for i=1:96
    plot(score(pre+(i-1)*93,1),score(pre+(i-1)*93,2),'Marker','.','Color','c', 'MarkerSize',10)
end
% plot(score(pre+((pre-1)*93),1),score(pre+((pre-1)*93),2),'ko')

pre=41;
for i=1:96
    plot(score(pre+(i-1)*93,1),score(pre+(i-1)*93,2),'Marker','.','Color','r', 'MarkerSize',10)
end
% plot(score(pre+((pre-1)*93),1),score(pre+((pre-1)*93),2),'ko')

pre=49;
for i=1:96
    plot(score(pre+(i-1)*93,1),score(pre+(i-1)*93,2),'Marker','.','Color','r', 'MarkerSize',10)
end
% plot(score(pre+((pre-1)*93),1),score(pre+((pre-1)*93),2),'ko')
% pre = 24;
% scorepre = (coeff(:,1:3)'*squeeze(Kernel(:,pre,:)))';
% hold on
% plot3(scorepre(:,1),scorepre(:,2),scorepre(:,3),'c.');
% 
% pre = 25;
% scorepre = (coeff(:,1:3)'*squeeze(Kernel(:,pre,:)))';
% hold on
% plot3(scorepre(:,1),scorepre(:,2),scorepre(:,3),'r.');
% 

% pre = 1;
% scorepre = (coeff(:,1:3)'*squeeze(Kernel(:,pre,:)))';
% hold on
% plot3(scorepre(:,1),scorepre(:,2),scorepre(:,3),'m.');
plot(score(1,1),score(1,2),'ko')
% plot(score(pre+(pre*93),1),score(pre+(pre*93),2),'ko')
%%
ChNorm = zeros(NoC,1);
for pre = 1:NoC
   scorepre = (coeff(:,1:3)'*squeeze(Kernel(:,pre,:)))';
   ChNorm(pre) = mean(sqrt(sum(scorepre.^2,2)));
end
figure
stem(ChNorm)

%%
Channels=1:96;
ChannelSet = setdiff(1:96,[66,20,88]);
% Channels=97:192;
% ChannelSet = setdiff(97:192, [174, 131, 133, 113]);

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

MEAMap= MEAMapPF;

%%
% ChNormRange = [0 0.275];
pre = 89;
scorepre = (coeff(:,1:3)'*squeeze(Kernel(:,pre,:)))';
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

ChNorm=[];
for ch=1:NoC
    scorepre = (coeff(:,1:3)'*squeeze(Kernel(:,ch,:)))';
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
scatter(distance_all,ChNorm,'c.')
% hold on
% plot(dist,smooth(meanNorm),'r.-')
%%
ChPF = setdiff(1:96,[66,20,88]);
group1PF = [96 8 44 46 36 41 50 80 73 81 93 65 67 95 78 84 71 79 94 68 82 83 77 47 76 69 75 74 70 48 72 45];%bimodali
group2PF = setdiff(ChPF,group1PF);
for i=1:numel(group1PF)
    if group1PF(i)>20 && group1PF(i)<60
        group1PF(i) = group1PF(i)-1;
    elseif group1PF(i)>60 && group1PF(i)<88
        group1PF(i) = group1PF(i)-2;
    elseif group1PF(i)>88
        group1PF(i) = group1PF(i)-3;
    end
end

for i=1:numel(group2PF)
    if group2PF(i)>20 && group2PF(i)<60
        group2PF(i) = group2PF(i)-1;
    elseif group2PF(i)>60 && group2PF(i)<88
        group2PF(i) = group2PF(i)-2;
    elseif group2PF(i)>88
        group2PF(i) = group2PF(i)-3;
    end
end
%%
ChNorm1=[];
ChNorm2=[];
n=0;
for ch=group1PF
    n=n+1;
    scorepre = (coeff(:,1:3)'*squeeze(Kernel(:,ch,ChannelSet)))';
    ChNorm1(:,n) = sqrt(sum(scorepre.^2,2));
end

for i=1:numel(dist)
    ind = find(distance_all(group1PF,:)==dist(i));
    meanNorm1(i) = mean(ChNorm1(ind));
end

n=0;
for ch=group2PF
    n=n+1;
    scorepre = (coeff(:,1:3)'*squeeze(Kernel(:,ch,ChannelSet)))';
    ChNorm2(:,n) = sqrt(sum(scorepre.^2,2));
end

for i=1:numel(dist)
    ind = find(distance_all(group2PF,:)==dist(i));
    meanNorm2(i) = mean(ChNorm2(ind));
end


figure
plot(dist,meanNorm1,'r.-')
hold on
plot(dist,meanNorm2,'b.-')
%%
% figure
% plot(dist,mean([meanNorm1; meanNorm2]),'r.-')
%%
% % ChNormRange = [0 0.275];
% post = 25;
% scorepost = (coeff(:,1:3)'*squeeze(Kernel(:,:,post)))';
% ChNorm = sqrt(sum(scorepost.^2,2));
% 
% figure
% [~,ndxCh2Map] = sort(MEAMap(:));
% ndxCh2Map = ndxCh2Map(5:end);
% MeasureMap = zeros(size(MEAMap));
% if MEAMap(2,2)>NoC
%     MeasureMap(ndxCh2Map(Channels-NoC)) = ChNorm;
% else
%     MeasureMap(ndxCh2Map(ChPF)) = ChNorm;
% end
% imagesc(MeasureMap);
% colorbar();
% hold on
% for channel=setdiff(ChannelSet,post)
%    [r,c] = find(MEAMap==channel);
%    text(c,r,num2str(channel),'HorizontalAlignment','center','VerticalAlignment','middle','Color','k')
% end
% [r,c] = find(MEAMap==post);
% text(c,r,num2str(post),'HorizontalAlignment','center','VerticalAlignment','middle','Color','r')
% 
% % caxis(ChNormRange)
% set(gca,'TickDir','out','Box','on','Layer','top')
% xlabel('X MEA')
% ylabel('Y MEA')

%%
% % ChNormRange = [0 1];
% pre = 1;
% scorepre = (coeff(:,1:3)'*squeeze(Kernel(:,pre,:)))';
% 
% figure
% [~,ndxCh2Map] = sort(MEAMap(:));
% ndxCh2Map = ndxCh2Map(5:end);
% MeasureMap = zeros(size(MEAMap));
% if MEAMap(2,2)>NoC
%     MeasureMap(ndxCh2Map(Channels-NoC)) = PerformanceTest;
% else
%     MeasureMap(ndxCh2Map(Channels)) = PerformanceTest;
% end
% imagesc(MeasureMap);
% colorbar();
% hold on
% for channel = ChannelSet
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
