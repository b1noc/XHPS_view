% init
%run('~/work/tools/XHPS_2021/HPS_simulation/hps_startup.m')
load('states.mat');
state = states;

%% settings

% sim_time %TODO: warning
sim_step = core_params.dt_sim ;
t_start = datenum(HPS_convertMJD2CalendarDate(core_params.start_date(1)));


vidName='time1'

% seconds between body frame plots
step = 20; 

% Vid quality
width = 1280;
height = 720;

% plot length
endframe=length(state);

%% Plot orbit with satellite body frame


% size of space around plot TODO: get from Vector lenght
scale = 1.2;


fig = figure('visible','off');
plot3(state(:,11),state(:,12),state(:,13))


set(fig, 'Position',  [100, 100, width, height])

x = xlim*scale;
y = ylim*scale;
xlim(x)
ylim(y)

%parfor_progress(floor(endframe/step));
m=floor(endframe/step);
fprintf('Progress:\n');
fprintf(['\n' repmat('.',1,m) '\n\n']);


% frame loop (plot state for each frame)
parfor i=1:floor(endframe/step)
	clf
  %D = parallel.pool.DataQueue


	k=i*step - step+1;
	t = t_start + (k-1) * seconds(sim_step);
    
	% Plot orbit line 
	plot3(state(:,11),state(:,12),state(:,13))
	hold on 

	%plot starting (initial) Point
	plot3(state(k,11),state(k,12),state(k,13),'color', 'red', 'Marker', 'o', 'MarkerSize', 7 , 'linewidth',2)

	%plot initial velocity vector
	%quiver3(state(k,11),state(k,12),state(k,13),state(1,8)*500,state(1,9)*500,state(1,10)*500,'linewidth',2)

	%plot inertial kos
	a=6771000;
	axis_length = 1*a/8; %length of coordinate system axis
	x_axis= [axis_length;0;0];
	y_axis= [0;axis_length;0];
	z_axis= [0;0;axis_length];

	quiver3(0,0,0,axis_length*1.5,0,0,'linewidth',3,'color', 'blue')
	quiver3(0,0,0,0,axis_length*1.5,0,'linewidth',3,'color', 'red')
	quiver3(0,0,0,0,0,axis_length*1.5,'linewidth',3,'color', 'green')


	%% Plot body fixed frame 
	q_k = state(k,4:7);

	%transform inertial kos to body fixed by quaternion multiplication
	%(xrvar qnefgryyhat qrf Beovgny Xbs in eci daher vecbyquattransposed)
	x_axis_body = HPS_transformVecByQuatTransposed(x_axis,q_k);
	y_axis_body = HPS_transformVecByQuatTransposed(y_axis,q_k);
	z_axis_body = HPS_transformVecByQuatTransposed(z_axis,q_k);


	% Vektor [1*a/8*1000,0,0] turn with quaternion 

	quiver3(state(k,11),state(k,12),state(k,13),x_axis_body(1),x_axis_body(2),x_axis_body(3),'linewidth',2,'color', 'blue')
	quiver3(state(k,11),state(k,12),state(k,13),y_axis_body(1),y_axis_body(2),y_axis_body(3),'linewidth',2,'color', 'red')
	quiver3(state(k,11),state(k,12),state(k,13),z_axis_body(1),z_axis_body(2),z_axis_body(3),'linewidth',2,'color', 'green')


	%% plot earth 

	grs80 = referenceEllipsoid('grs80','m');

	ax = axesm('globe','Geoid',grs80,'Grid','on', ...
	'GLineWidth',1,'GLineStyle','-',...
	'Gcolor','black','Galtitude',100);
	box off
	view([90 8]);

	%set(gca, 'CameraPosition', [100 5000 2000]);

	xlim(x);
	ylim(y);
	xlabel(datestr(t)) 

	aniVec(i) = getframe(fig, [0 0 width height]);
    
	fprintf('\b#\n');

end

save('aniVec','aniVec', '-v7.3')

wout = VideoWriter(vidName);
wout.FrameRate = 5;
open(wout);
writeVideo(wout,aniVec);
close(wout);



