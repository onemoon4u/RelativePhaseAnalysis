%
%
% Script to perform PCA analysis on Relative Phase Dynamics
%
%
% Code written by Joon-Young Moon
% Final update on 2025-July-30th



%% load file and setup
% This is the file produced by "cluster_movie_frames_eval.m"
% Change directory names accordingly

file_name = 'michigan_example_data_tm50_tw50_sm20_relp_centroidK';
load( file_name, 'H', 'time_moving', 'time_window', 'Fs', 'rel_p', 'rel_phase_w_mean',  'chan_coord_xy','time_all', 'time_pt', 'time_moving', 'time_window', 'topo','smooth','centroid_K','centroid_K_vector', 'K', 'IDX')



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

%% spatial PCA analysis: computation 

[coeff,score,latent,tsquared,explained,mu] = pca(topo_vector);    % Utilizes "pca" from MATLAB function. One can perform PCA with singular value decomposition.


% If you don't need all PCA components, you can save time by reducing number of PCs to be calculated. 
% For example, if you only want to calculate frist 25 PCs:
%
%
% [coeff_r, score_r, latent_r,tsquared_r,explained_r,mu_r] = pca(topo_vector, 'NumComponents', 25);
% similarity = abs(dot(coeff(:,1), coeff_r(:,1)) / (norm(coeff(:,1)) * norm(coeff_r(:,1))));
% disp(['Cosine Similarity: ', num2str(similarity)]);


% Utilizing singular value decomposition function "svd" instead of "pca":
%
% topo_vector_centered = topo_vector - mean(topo_vector);
% 
% [U,S,V] = svd(topo_vector_centered)  or X = USV' where "X" is "topo_vector_centered"
% coeff = V
% score = XV = US


%% spatial PCA analysis: plot each data points (each time points) in PC axises
% each dot in the figure corresponds to each frame of the movie, showing
% how they are located in the PC space, a 3D space constructed using PC1, PC2, PC3 
% The star denotes the position of each centroid in the PC space.
% The colors denote K-means clustered groups.

figure(501);
plot(explained,'o-','LineWidth',2,'MarkerSize',10)
xlabel('principal compnent index')
ylabel('percentage')
title('percentage of the total variance explained')
xlim([1 25])
ylim([0 60])


explained_cdf = zeros(1,length(explained));
explained_cdf(1) = explained(1);
for h=2:length(explained)
    explained_cdf(h)=explained(h)+explained_cdf(h-1);
end

figure(503);
plot(explained_cdf,'o-','LineWidth',2,'MarkerSize',10)
xlabel('principal component index')
ylabel('cumulative percentage')
title('cumulative percentage of the total variance explained')
xlim([1 25])
ylim([0 100])
grid on;


centroid_K_score = centroid_K_vector*coeff;

figure(505);
for i=1:max(IDX)
    temp_idx = find(IDX==i);
    plot3(score(temp_idx,1), score(temp_idx,2), score(temp_idx,3),'o');
    hold on;
    grid on;
    axis on;
    xlabel('PC1');ylabel('PC2');zlabel('PC3')
    axis equal;
    axis ij;
%     xlim([-150,150]);ylim([-100,100]);zlim([-100,100]);
 %   xline(0)
  %  yline(0)
    xL = xlim;yL = ylim;zL = zlim;
line([0 0], yL, 'color', 'k' , 'linewidth', 2);  %x-axis
line(xL, [0 0], 'color', 'k' , 'linewidth', 2); 
line([0 0 ], [0 0], zL, 'color', 'k' , 'linewidth', 2 );
end
color_rgb = [ 0 0.4470 0.7410 ; 0.8500 0.3250 0.0980 ; 0.9290 0.6940 0.1250 ; 0.4940 0.1840 0.5560 ; 0.4660 0.6740 0.1880 ; 0.3010 0.7450 0.9330 ; 0.6350 0.0780 0.1840 ] ;
figure(505);
for i=1:K
    hold on;
    plot3(centroid_K_score(i,1),centroid_K_score(i,2),centroid_K_score(i,3),'pentagram','MarkerFaceColor',color_rgb(i,:), 'MarkerEdgeColor' ,'k' , 'MarkerSize',25);
end





%% plot principal components

L=4; 
% If we set L = K/2, we will get PCs corresponding to the K centroids. 
% Usually, K centroids contains postive and negative versions of PCs, hence the formula L= K/2 results.


shiftpreset = 0;
% 0 if it is human
% 1 if it is monkey

f1 = figure(100);
[dataOut, xx, yy,Coord,borderCoords] = topoplot_general_test( zeros(1,size(chan_coord_xy,1))', chan_coord_xy(:,1:2),'smooth',smooth,'shiftpreset', shiftpreset, 'scatter', 1);
close(f1)

pca_L = nan(size(temp1,1),size(temp1,2),L);

for j=1:L
    for i=1:length(topo_idx_x)
        pca_L(topo_idx_x(i),topo_idx_y(i),j) = coeff(i,j); 
    end
end

for j=1:L
    figure(1000+j)
    topoplot_figure(pca_L(:,:,j), borderCoords, xx, yy, Coord, 'scatter', 1);
end


%% plot of centroids for comparison with principal components
L=4;

f1 = figure(200);
[dataOut, xx, yy,Coord,borderCoords] = topoplot_general_test( zeros(1,size(chan_coord_xy,1))', chan_coord_xy(:,1:2),'smooth',smooth,'shiftpreset', shiftpreset, 'scatter', 1);
close(f1)


cent_L = nan(size(temp1,1),size(temp1,2),L);

for j=1:L
    for i=1:length(topo_idx_x)
        cent_L(topo_idx_x(i),topo_idx_y(i),j) = centroid_K_vector(j,i); 
    end
end

for j=1:L
    figure(2000+j)
    topoplot_figure(cent_L(:,:,j), borderCoords, xx, yy, Coord, 'scatter', 1);
end

%% Save the data with the PCs produced

save( [ file_name '_PCA' ] ,'H', 'Fs', 'rel_p', 'rel_phase_w_mean',    'chan_coord_xy','time_all', 'time_pt', 'time_moving', 'time_window', 'topo','smooth','centroid_K','centroid_K_vector', 'K', 'IDX', 'coeff','score','latent','tsquared','explained','-v7.3' );