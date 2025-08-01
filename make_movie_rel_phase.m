%
%
% Script to run RPA and produce Relative Phase Dynamics Movie
%
%
% Code written by Joon-Young Moon
% Final update on 2025-July-28th


clear all;

%% Define file name and load the file

file_name = 'michigan_example_data';
load(file_name,"data_eye_closed","Coords");
H = data_eye_closed';
Fs = 500;  % sampling frequency


% Make sure that the first column is the x coordinate and the second column is the y coordinate
% The easy way to confirm is:  
% plot( coords(:,1) , coords(:,2), 'o')
% See if the above plot give you left-right symmetry. If not, change the order of the coords(:,1) and coords(:,2)

chan_coord_xy = Coords;


%% Compute relative phase
[rel_phase_w, rel_phase_t, rel_phase, magni, R_theta, T_theta, m_theta] = cal_rel_phase_v_final(H) ;   % "cal_rel_phase_v_final" performs relative phase calculation


%% Set up parameters for the movie
smooth = 20; % spatial smoothness: 20 provides a good balance between local details and global patterns 
tend_temp = 1; % portion of the time which you want to make the movie with
tend = 1./tend_temp;


%% Choosing the temporal resolution by performing time window averaging 
% For the optimal resolution considering computational cost yet capturing all the essential dynamics, resolusion of 20 ms is recommanded 

% Moving time window
time_moving = 50;  %  number of time samples. time_moving / Fs = moving time in seconds. For resolution of 20 ms, and Fs = 500, set the value to 10
time_window = 50;  %  number of time samples. time_window / Fs = window time in seconds. For resolution of 20 ms, and Fs = 500, set the value to 10
dt=1;

[ rel_phase_w_mean ,time_pt] = moving_time_window( rel_phase_w, time_moving, time_window);  % "moving_time_window" performs time window averaging

rel_p = double( rel_phase_w_mean(1:round(size(rel_phase_w_mean,1)./tend),: ) ) ;

%% movie making core code

vidfile = VideoWriter([ file_name '_tm' num2str(time_moving) '_tw' num2str(time_window) '_sm' num2str(smooth) '_relp_test1' ] ,'MPEG-4');
vidfile.FrameRate = round(Fs/time_moving/dt);
vidfile.Quality = 75; % Value of 75 provides a good balance between computational cost and quality. Modify accordingly with respect to the computational power
open(vidfile);    
FrameRate = round(Fs/time_moving/dt);

time_dt = 1/Fs*time_moving;
time_now = 0 + (time_window)/2 /Fs;
time_all = [1:dt:round(size(rel_p,1)./tend )]';

topo = cell(size(time_all,1),1);
formatSpec = '%.2f';

tic
for i=1:dt:(size(rel_p,1))
    
figure(1);
topo{i} = topoplot_general_test(rel_p(i,:)', chan_coord_xy(:,1:2),'smooth',smooth,'scatter',1);   % "topoplot_general_test" is the key code to produce the relative phase topoplot
title(  [ '\rm' num2str(time_now,formatSpec) ' S'; ] ,'fontsize',16 )
drawnow
time_all(i)=time_now;
F(i) = getframe(gcf); 
writeVideo(vidfile,F(i));
time_now = time_now + dt*time_dt;
end
close(vidfile)
toc

% saving necessary data in case one needs to produce the movie again at later time without computing relative phase again  
save( [ file_name '_tm' num2str(time_moving) '_tw' num2str(time_window) '_sm' num2str(smooth)  '_relp' ] ,'H', 'Fs', 'rel_p', 'rel_phase_w_mean',    'chan_coord_xy','time_all', 'time_pt', 'time_moving', 'time_window', 'topo','smooth','-v7.3' );