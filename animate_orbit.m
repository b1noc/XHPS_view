%% settings

debug=0;

% seconds between body frame plots
step = 100; 


%% Plot orbit with satellite body frame


% size of space around plot TODO: get from Vector lenght
scale = 1.2;

state = states;
fig = figure;

% animation frame counter
k=0;

endframe=length(state);
if debug
	endframe= step;	
else
	clear aniVec
end

% frame loop (plot state for each frame)
for i=1:step:endframe

	k=k+1;
	clf

	% Plot orbit line 
	plot3(state(:,11),state(:,12),state(:,13))
	hold on 

	%plot starting (initial) Point
	plot3(state(i,11),state(i,12),state(i,13),'color', 'red', 'Marker', 'o', 'MarkerSize', 7 , 'linewidth',2)

	%plot initial velocity vector
	%quiver3(state(i,11),state(i,12),state(i,13),state(1,8)*500,state(1,9)*500,state(1,10)*500,'linewidth',2)

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
	q_i = state(i,4:7);

	%transform inertial kos to body fixed by quaternion multiplication
	%(xrvar qnefgryyhat qrf Beovgny Xbs in eci daher vecbyquattransposed)
	x_axis_body = HPS_transformVecByQuatTransposed(x_axis,q_i);
	y_axis_body = HPS_transformVecByQuatTransposed(y_axis,q_i);
	z_axis_body = HPS_transformVecByQuatTransposed(z_axis,q_i);


	% Vektor [1*a/8*1000,0,0] turn with quaternion 

	quiver3(state(i,11),state(i,12),state(i,13),x_axis_body(1),x_axis_body(2),x_axis_body(3),'linewidth',2,'color', 'blue')
	quiver3(state(i,11),state(i,12),state(i,13),y_axis_body(1),y_axis_body(2),y_axis_body(3),'linewidth',2,'color', 'red')
	quiver3(state(i,11),state(i,12),state(i,13),z_axis_body(1),z_axis_body(2),z_axis_body(3),'linewidth',2,'color', 'green')


	%% plot earth 

	grs80 = referenceEllipsoid('grs80','m');

	ax = axesm('globe','Geoid',grs80,'Grid','on', ...
	'GLineWidth',1,'GLineStyle','-',...
	'Gcolor','black','Galtitude',100);
	sax.Position = [0 0 1 1];
	box off
	view([90 8])
	%set(gca, 'CameraPosition', [100 5000 2000]);

	if(i==1)
		drawnow;
		pos = get(fig,'Position');
		width = pos(3);
		height = pos(4);
		x = xlim*scale;
	 	y = ylim*scale;
	end

	xlim(x);
	ylim(y);

	aniVec(k) = getframe(fig, [10 10 width-10 height-10]);

end

	wout = VideoWriter('orbit_animation');
	wout.FrameRate = 3;
	open(wout);
	writeVideo(wout,aniVec);
	close(wout);

clear k step x_axis_body y_axis_body z_axis_body q_i i_time i state z_axis y_axis x_axis axis_length
