%
%
% Script to perform K-means clustering on Relative Phase Dynamics
%
%
% Code written by Joon-Young Moon
% Final update on 2025-July-28th
%
% If you already have pre-constructed centroids, you don't need to run this
% script.



%% load file and setup
% This is the file produced by "make_movie_rel_phase_v_final.m"
% Change directory names accordingly

file_name = 'michigan_example_data_tm50_tw50_sm20_relp';
load( file_name, 'H', 'time_moving', 'time_window', 'Fs', 'rel_p', 'rel_phase_w_mean',  'chan_coord_xy','time_all', 'time_pt', 'time_moving', 'time_window', 'topo','smooth')


%% Set up vectors
% vectorize data
topo_size = length(topo);

temp1 = topo{1};
topo_idx = isnan(temp1);
[topo_idx_x,topo_idx_y] = find(topo_idx == 0);
topo_vector = zeros(topo_size,length(topo_idx_x));

% topo_vector: time point x frame vector
for i=1:topo_size

    T = topo{i};
    T(isnan(T)) = [] ;
    topo_vector(i,:) = T;

end


%% This step is to test for the optimal K 
% This step takes very long time, skip if you are certain about the number K that you want to try
% tic
% eva = evalclusters(topo_vector,'kmeans','CalinskiHarabasz','KList',[1:10]);
% toc

%% K-means clustering
% Main part for K-means clustering

% set up the parameters

K=4;
test_num = 10;

% k-mean clustering core

    tic
    disp( [ 'number of iterations=' num2str(test_num) ] )
    [IDX_old,C_old,SUMD_old,D_old]=kmeans(topo_vector,K, 'distance', 'sqeuclidean','Replicates',test_num,'Display','final','emptyaction','drop');   


   % [s,h] = silhouette(topo_vector,IDX);
   % S = mean(s);
    toc


figure(3)
hist = histogram(IDX_old,'normalization','probability');
counts = hist.Values;

topo_vector_idx = find(topo_idx'==0);
topo_length = size(topo{1},2);
[IDX, C, SUMD, D, rho] = change_cluster_idx_v_final(K, topo_length, topo_vector_idx, IDX_old, C_old, SUMD_old, D_old);

figure(4)
histogram(IDX,'normalization','probability');

figure(5)
plot(time_all, IDX,'.-')
    

%% get centroids

shiftpreset = 0; % 0 for human , 1 for macaque
f1 = figure(100);
[dataOut, xx, yy,Coord,borderCoords] = topoplot_general_test(rel_p(1,:)', chan_coord_xy(:,1:2),'smooth',smooth,'shiftpreset', shiftpreset, 'scatter', 1);
close(f1)

centroid_K = nan(size(temp1,1),size(temp1,2),K);

for j=1:K
    for i=1:length(topo_idx_x)
        centroid_K(topo_idx_x(i),topo_idx_y(i),j) = C(j,i); 
    end
end

for j=1:K
    figure(100+j)
    topoplot_figure(centroid_K(:,:,j), borderCoords, xx, yy, Coord, 'scatter', 1);
end

centroid_K_vector = zeros(K,length(topo_idx_x)) ;
for i=1:K
    centroid_K_temp = centroid_K(:,:,i);
    centroid_K_temp(isnan(centroid_K_temp)) = [] ;
    centroid_K_vector(i,:) = centroid_K_temp;
end



%% Use TSNE algorithm to see if the clustering is well performed
% skip if unnecessary

tic
 tsne_topo = tsne(topo_vector, 'Algorithm','exact' , 'Distance','seuclidean','NumDimensions',2);
 figure(31)
 gscatter(tsne_topo(:,1),tsne_topo(:,2),IDX);
axis square;
xlim([-125 125]); ylim([-125 125]);
toc

% tic
%     Y = pdist(topo_vector);
%     YS = squareform(Y);
%     Z = linkage(Y);
%     dendrogram(Z);
% toc

%% dPLI and dwPLI computation, and ploting topoplot based on dPLI and dwPLI
% This part computes various other measures
% Utilizes "all_phase_measures_ht","d_PhaseLagIndex2","dw_PhaseLagIndex" made by Heonsoo Lee

[coh, mpc,imc ]=all_phase_measures_ht(H);
coh(1:1+size(coh,1):end) = nan;
fC_coh = double(nanmean(coh));
mpc(1:1+size(mpc,1):end) = nan;
fC_mpc = double(nanmean(mpc));
imc(1:1+size(imc,1):end) = nan;
fC_imc = double(nanmean(imc));


dPLI=d_PhaseLagIndex2(H);
PLI = abs(dPLI);
PLI(1:1+size(PLI,1):end) = nan;
fC_PLI = nanmean(PLI);

dwPLI=dw_PhaseLagIndex(H);
wPLI = abs(dwPLI);
wPLI(1:1+size(wPLI,1):end) = nan;
fC_wPLI = nanmean(wPLI);

figure(41);
topoplot_general_test(fC_coh, chan_coord_xy(:,1:2),'smooth',smooth,'shiftpreset', shiftpreset, 'scatter', 1);
title('Functional Connectivity by coherence')

figure(42);
topoplot_general_test(fC_mpc, chan_coord_xy(:,1:2),'smooth',smooth,'shiftpreset', shiftpreset, 'scatter', 1);
title('Functional Connectivity by phase coherence')

figure(43);
topoplot_general_test(fC_imc, chan_coord_xy(:,1:2),'smooth',smooth,'shiftpreset', shiftpreset, 'scatter', 1);
title('Functional Connectivity by imaginary coherence')



figure(44);
topoplot_general_test(fC_PLI, chan_coord_xy(:,1:2),'smooth',smooth,'shiftpreset', shiftpreset, 'scatter', 1);
title('Functional Connectivity by PLI')

figure(45);
topoplot_general_test(fC_wPLI, chan_coord_xy(:,1:2),'smooth',smooth,'shiftpreset', shiftpreset, 'scatter', 1);
title('Functional Connectivity by wPLI')


%% Looking at the time series of "top-down" vs "bottom-up" fluctuation 

corr_pf = zeros(1,size(rel_p,1));
for i=1:size(rel_p,1)
   corr_pf(i) = corr(rel_p(i,:)',fC_imc','type','Spearman');
end

figure(46)
plot(time_all,corr_pf)
xlabel('time(s)')
ylabel('corr(phase,fC by imc)')
title('Time series for top-down/bottom-up switching, by corr(phase,imc)')


%% Save the data with the centroids produced

save( [ file_name '_centroidK' ] ,'H', 'Fs', 'rel_p', 'rel_phase_w_mean',    'chan_coord_xy','time_all', 'time_pt', 'time_moving', 'time_window', 'topo','smooth','centroid_K','centroid_K_vector', 'K', 'IDX', '-v7.3' );