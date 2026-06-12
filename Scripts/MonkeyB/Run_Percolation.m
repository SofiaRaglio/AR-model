clear all;%close all; clc;
%%
load('Wakefulness.mat')
%%
Ch = setdiff(1:192,[66,20,88,174,113,131,133]);
% Ch = setdiff(1:192,[66,20,88,174,113,131,133]);
ind = round(size(LFP_W1,2)/16)-1;
%%
for i=1:16
    ndxT = 1+ind*(i-1):ind+ind*(i-1);
    Matrix(:,:,i) = corr(MUA_W1(Ch,ndxT)',MUA_W1(Ch,ndxT)');
    [perc_threshold,Veglia.Perc_Matrix(:,:,i), Veglia.Sort_value{i},Veglia.n_com_size{i},perc_threshold_step] = Percolation(Matrix(:,:,i));
    Veglia.Perc_Matrix(:,:,i) = tanh(Veglia.Perc_Matrix(:,:,i));
end
Veglia.Perc_Matrix = atanh(mean(Veglia.Perc_Matrix,3));
%%
Veglia.Perc_Matrix_pruned =[];
for i=1:size(Veglia.Perc_Matrix,1)
    for j=1:size(Veglia.Perc_Matrix,2)
        if Veglia.Perc_Matrix(i,j)<0.5
            Veglia.Perc_Matrix_pruned(i,j)=0;
        else
            Veglia.Perc_Matrix_pruned(i,j)=Veglia.Perc_Matrix(i,j);
        end
    end
end
%%
d  = graph(Veglia.Perc_Matrix_pruned,'upper');
figure;
h = plot(d,'LineStyle','--','MarkerSize',5, 'NodeColor', 'g');
highlight(h, 1:93, 'NodeColor', 'r');
%%
figure;
imagesc(Veglia.Perc_Matrix_pruned)
%%
load('Anaesthesia.mat')
%%
Ch = setdiff(1:192,[66,20,88,174,113,131,133]);
ind = round(size(LFP_SO,2)/16)-1;
%%
for i=1:16
    ndxT = 1+ind*(i-1):ind+ind*(i-1);
    Matrix(:,:,i) = corr(MUA_SO(Ch,ndxT)',MUA_SO(Ch,ndxT)');
    [perc_threshold,Anestesia.Perc_Matrix(:,:,i), Anestesia.Sort_value{i},Anestesia.n_com_size{i},perc_threshold_step] = Percolation(Matrix(:,:,i));
    Anestesia.Perc_Matrix(:,:,i) = tanh(Anestesia.Perc_Matrix(:,:,i));
end
Anestesia.Perc_Matrix = atanh(mean(Anestesia.Perc_Matrix,3));
%%
Anestesia.Perc_Matrix_pruned =[];
for i=1:size(Anestesia.Perc_Matrix,1)
    for j=1:size(Anestesia.Perc_Matrix,2)
        if Anestesia.Perc_Matrix(i,j)<0.5
            Anestesia.Perc_Matrix_pruned(i,j)=0;
        else
            Anestesia.Perc_Matrix_pruned(i,j)=Anestesia.Perc_Matrix(i,j);
        end
    end
end
%%
d  = graph(Anestesia.Perc_Matrix_pruned,'upper');
figure;
h = plot(d,'LineStyle','--','MarkerSize',5, 'NodeColor', 'g');
highlight(h, 1:93, 'NodeColor', 'r');
%%
figure;
imagesc(Anestesia.Perc_Matrix_pruned)
%%
load('Awakening.mat')
%%
Ch = setdiff(1:192,[66,20,88,174,113,131,133]);
ind = round(size(LFP_Aw,2)/16)-1;
%%
for i=1:16
    ndxT = 1+ind*(i-1):ind+ind*(i-1);
    Matrix(:,:,i) = corr(MUA_Aw(Ch,ndxT)',MUA_Aw(Ch,ndxT)');
    [perc_threshold,Risveglio.Perc_Matrix(:,:,i), Risveglio.Sort_value{i},Risveglio.n_com_size{i},perc_threshold_step] = Percolation(Matrix(:,:,i));
    Risveglio.Perc_Matrix(:,:,i) = tanh(Risveglio.Perc_Matrix(:,:,i));
end
Risveglio.Perc_Matrix = atanh(mean(Risveglio.Perc_Matrix,3));
%%
Risveglio.Perc_Matrix_pruned =[];
for i=1:size(Risveglio.Perc_Matrix,1)
    for j=1:size(Risveglio.Perc_Matrix,2)
        if Risveglio.Perc_Matrix(i,j)<0.5
            Risveglio.Perc_Matrix_pruned(i,j)=0;
        else
            Risveglio.Perc_Matrix_pruned(i,j)=Risveglio.Perc_Matrix(i,j);
        end
    end
end
%%
d  = graph(Risveglio.Perc_Matrix_pruned,'upper');
figure;
h = plot(d,'LineStyle','--','MarkerSize',5, 'NodeColor', 'g');
highlight(h, 1:93, 'NodeColor', 'r');
%%
figure;
imagesc(Risveglio.Perc_Matrix_pruned)

%%
figure
for i=1:16
    plot(Veglia.Sort_value{i},Veglia.n_com_size{i},'r')
    hold on
    plot(Anestesia.Sort_value{i},Anestesia.n_com_size{i},'b')
    plot(Risveglio.Sort_value{i},Risveglio.n_com_size{i},'c')
    % plot(Veglia_finale.Sort_value{i},Veglia_finale.n_com_size{i},'m')
end
%%

cross_veglia =[];
cross_an =[];
cross_ris =[];
for i=1:16
    n=0;
    for j = 80:5:120
        n=n+1;
        nd_veglia = find(Veglia.n_com_size{i}==j,1);
        cross_veglia(i,n) = Veglia.Sort_value{i}(nd_veglia);

        nd_an = find(Anestesia.n_com_size{i}==j,1);
        cross_an(i,n) = Anestesia.Sort_value{i}(nd_an);

        nd_ris = find(Risveglio.n_com_size{i}==j,1);
        cross_ris(i,n) = Risveglio.Sort_value{i}(nd_ris);
    end
end

%%
figure;
errorbar([mean(std(cross_veglia)) mean(std(cross_an)) mean(std(cross_ris))], [std(std(cross_veglia)) std(std(cross_an)) std(std(cross_ris))],'o-')
xlim([0.5 3.5])
