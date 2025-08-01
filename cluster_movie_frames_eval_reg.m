

%% Script to perform clustering using "universal centroids" on Relative Phase Dynamics
% This script will cluster each frame of Relative Phase Dynamics Movie (topoplots) using linear regression,
% using the universal centroids as the regressor.
% This is actually equivalent to performing PCA with the universal centroids as the PCs (Principal components), or eigenvectors.
%
% Code written by Joon-Young Moon, Jehyup Lee
% Final update on 2025-August-1st
%
% For now, the given universal centroids numbers are 4 (K=4):
% There are 4 canonical modes.
% We will update the code in cases where N>4 centroids are to be used in the near future.  



%% load file and setup
% This is the file produced by "cluster_movie_frames_eval.m"
% Change directory names accordingly

file_name = 'michigan_example_data_tm50_tw50_sm20_relp_centroidK';
load( file_name, 'H', 'time_moving', 'time_window', 'Fs', 'rel_p', 'rel_phase_w_mean',  'chan_coord_xy','time_all', 'time_pt', 'time_moving', 'time_window', 'topo','smooth','centroid_K','centroid_K_vector', 'K' )
% The code assumes K=4 (4 clusters, 4 centroids) for now

% universal centroids
% using eye_close (EO), eye_open data (EC), and combination of these
file_name2='combined_centroids_20240311';
load(file_name2); 


%% Setting up the axis from canonical modes

% "EOEC" data is from EO (eyes-open) and EC (eyes-closed), combined together
% There are 4 canonical modes, and mode 3,4 are direct opposite of mode 1,2
mask=C_comb.EOEC; 

% Taking mean of mode 1 and -1*mode 4 together, and vectorizating the mean:
% Equivalent to constructing the 1st Principal Component (PC1), the first eigenvoctor 
mask([3 4],:,:) = -1.*mask([3 4],:,:);

x1 = squeeze(mean(mask([1 4],:,:),1));
x1 = x1(topo_idx_comb); 
x1 = x1./norm(x1);

% Taking eanof mode 2 and -1*mode 3 together, and vectorizating the mean:
% Equivalent to constructing the 2nd Principal Component (PC2), the second eigenvector
x2 = squeeze(mean(mask([2 3],:,:),1)); 
x2 = x2(topo_idx_comb); 
x2 = x2./norm(x2);

% Making axis from the canonical vectors constructed above
X = cat(2,x1,x2);


%% Setting up the topoplot data (constructed before)
% Set up the topoplot data (each frame of Relative Phase Movie) in order to perform linear regression in the next part

% vectorize the topoplot data
 
U_size=C_comb.EOEC;
U_size1=squeeze(U_size(1,:,:));
f_nan=isnan(U_size1);
mi_size=topo{1};
mi_size(f_nan)=NaN;
temp1 = mi_size;

topo_idx = isnan(temp1);
topo_size = length(topo);
[topo_idx_x,topo_idx_y] = find(topo_idx == 0);
topo_vector = zeros(topo_size,length(topo_idx_x));

% topo_vector: time point x frame vector
for i=1:topo_size

    T = topo{i};
    T(f_nan)=NaN;
    T(isnan(T)) = [] ;
    topo_vector(i,:) = T;
end

%% Linear regression
% Here, linear regression will be performed using the canonical vectors (eigenvectors) produced above as the regressor.
% "beta" is equivalent to the coordinates of each frame of the Relative Phase Movie in the PC space.
% "E" is the error vector: the error of which each frame cannot be explained from given eigenvectors.

y=topo_vector';
tic;
for t = 1:size(y,2)
    [beta(1:2,t),~,E] = mvregress(X,y(:,t));
    beta(3,t) = mean(E).*(10^2);
end
toc;

%% Clustering each frame of the Relative Phase Movie (topoplots) into one of the 4 modes (mode1, mode2, mode3, mode4)
% This is essentially quivalent to K-means clustering.
% Also, this is equivalent to clustering each data point into N groups, by their coordinates in the PC space.
% Again, for now, we assume K=4.

K = 4; 

unit = 0.1; n_nulls = 10000; % unit is the time resolution of the data, in seconds. n_nulls is the # of random rotations for null model of transtions
tic
[IDX, w_prop] = cal_regression_clustering(beta, K, unit, n_nulls); 
toc
prop = cal_transition_prop_v_final(IDX, K, unit);

% "IDX" gives the indices of each frame into one of the K clusters
% "w_prop" gives the transition probability between k clusters, for a null model to compare with "prop" 
% "prop" gives the transition probability between K clusters, amounting to transition matrix.

%% plot
% plotting each frame of the Relative Phase Movie in the PC space,
% constructed by two PCs (eigenvectors, regressiors).
% Each data point represent each topoplot (each frame of the Relative Phase Movie).
% The X-axis corresponds to PC1, Y-axis corresponds to PC2, and z-axis corresponds to |E|.
% Z-axis will present the degree of error for each topoplot: the degree where each frame cannot be explained from the first two PCs.

figure(601);
for i=1:max(IDX)
    temp_idx = find(IDX==i);
    plot3(beta(1,temp_idx), beta(2,temp_idx), beta(3,temp_idx),'o');
    hold on;
    grid on;
    axis on;
    xlabel('PC1');ylabel('PC2');zlabel('|E|X10^2')
    axis equal;
    axis ij;
 %  xlim([-150,150]);ylim([-100,100]);zlim([-100,100]);
 %  xline(0)
 %  yline(0)
    xL = xlim;yL = ylim;zL = zlim;
line([0 0], yL, 'color', 'k' , 'linewidth', 2);  %x-axis
line(xL, [0 0], 'color', 'k' , 'linewidth', 2); 
line([0 0 ], [0 0], zL, 'color', 'k' , 'linewidth', 2 );
end
color_rgb = [ 0 0.4470 0.7410 ; 0.8500 0.3250 0.0980 ; 0.9290 0.6940 0.1250 ; 0.4940 0.1840 0.5560 ; 0.4660 0.6740 0.1880 ; 0.3010 0.7450 0.9330 ; 0.6350 0.0780 0.1840 ] ;

% figure(506);
% for i=1:K
%     hold on;
%     plot3(centroid_K_score(i,1),centroid_K_score(i,2),centroid_K_score(i,3),'pentagram','MarkerFaceColor',color_rgb(i,:), 'MarkerEdgeColor' ,'k' , 'MarkerSize',25);
% end

%% Topoplots of the Universal Centroids
%
% These plots give the shape of universal centroids

f1 = figure(700);
shiftpreset=0;
[dataOut, xx, yy,Coord,borderCoords] = topoplot_general_test( zeros(1,size(chan_coord_xy,1))', chan_coord_xy(:,1:2),'smooth',smooth,'shiftpreset', shiftpreset, 'scatter', 1);
close(f1)
for j=1:K
    figure(700+j)
    topo=squeeze(C_comb.EOEC(j,:,:));
    topoplot_figure(topo, borderCoords, xx, yy, Coord, 'scatter', 1);
    name=['mode',num2str(j)];
    title(name,'FontSize', 14)
end

%% Save the data with the beta, E, IDX and transition probabilities produced

save( [ file_name '_reg' ] ,'beta','prop','H', 'time_moving', 'time_window', 'Fs', 'rel_p', 'rel_phase_w_mean',  'chan_coord_xy','time_all', 'time_pt', 'time_moving', 'time_window', 'topo','smooth','centroid_K','centroid_K_vector', 'K', 'beta', 'E', 'IDX','prop', 'w_prop','-v7.3' );