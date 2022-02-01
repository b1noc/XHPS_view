function frameVec = generateFrames(states, settings, debug)

%% Launching parpool for multithreadding
if isempty(gcp('nocreate')) && debug == 0
   parpool;
end

%% Setting settings
% Defining default values for unspecified user settings
def = struct(	'step', 1, ...
				'resolution', '1080p', ...
				'sim_step', 10, ...
				't_start', 730486, ...
				'borderScale', 1.2, ...
				'vecLength', 500, ... %TODO: factor of zoom
				'camMode', 'fixed', ...
				'zoom', 50, ...
				'earthTransparency', 1, ...
				'satelliteModel', 'none', ...
				'viewAngle', [0 90], ...
				'progress', 1, ...
				'stlfactor', 1e-1, ...
				'stlfile', 'res/jeb.stl', ...
				'stlBaseColor', 'black', ...
				'stlEdgeColor', 'green', ...
				'earthfile', 'res/earth.jpg', ...
				'earthPanels', 180 ...
			);

				%TODO: rotate image file
%TODO: image_file loaded from 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';


% Overwriting default settings with specified user settings
diff = intersect(fieldnames(settings),fieldnames(def));
for n = 1:length(diff)
	def.(diff{n}) = settings.(diff{n});
end

%% Reading in states
% only reading every stepth value of states
endframe=length(states);
state = states(1:def.step:endframe,:);

% optimizing for multithreadded processing
vSatPos = state(:,11:13);
vSatOri = state(:,4:7);
vSatVel = state(:,8:10);

fv = stlread(def.stlfile);

%% Defining vido quality based on user settings
switch def.resolution
	case '720p'
		width = 1280;
		height = 720;
	case '1080p'
		width = 1920;
		height = 1080;
end

%% DEBUG-mode
% Display loaded stl file in figure (DEBUG-mode)
if debug == 1
    points=fv.Points*def.stlfactor;
	ship = trimesh(fv.ConnectivityList, points(:,1),points(:,2),points(:,3));
	set(ship, 'FaceColor', def.stlBaseColor, 'EdgeColor', def.stlEdgeColor);
	view([0 0])
	axis equal
end

% Display video figure in DEBUG-mode
if debug > 0
    fig = figure;
else
    fig = figure('visible','off');
end

%% Prepare looprun
plot3(vSatPos(:,1),vSatPos(:,2),vSatPos(:,3))
set(fig, 'Position',  [100, 100, width, height]);

% scaling axes
xl = xlim*def.borderScale;
yl = ylim*def.borderScale;
xlim(xl)
ylim(yl)

% define last frame 
ie = ceil(endframe/def.step);

% initializing vector to save frames
frameVec = struct('cdata', cell(1, ie), 'colormap', cell(1, ie));


%% printing progress bar (multithread compatible)
if def.progress
	fprintf('Progress:\n');
	fprintf(['\n' repmat('.',1,ie) '\n\n']);
end

%% Main loop
for i=1:ie
    clf
    
%% Calculate time
	% Calculating index k from original states matrix
    k=i*def.step - def.step;
	% Calculating timestamp of current state by adding x*sim_step to the starttime
    tstr = datestr(def.t_start + k * seconds(def.sim_step));
    t = def.t_start + k * def.sim_step;

%% Calculate satellite and camera position vectors
	% velocity vector
    velocity = [vSatVel(i,1)*def.vecLength,vSatVel(i,2)*def.vecLength,vSatVel(i,3)*def.vecLength];
	vn = velocity./norm(velocity);
	vfront = vn*def.zoom;
    % position of satelite
    satPos = [vSatPos(i,1) vSatPos(i,2) vSatPos(i,3)];
    sn = satPos./norm(satPos);
    satView = satPos+sn*def.zoom;

    orbPlaneA = [vSatPos(1,1) vSatPos(1,2) vSatPos(1,3)] ;
    orbPlaneB = [vSatPos(2,1) vSatPos(2,2) vSatPos(2,3)] ;
    orbPlaneC = [vSatPos(3,1) vSatPos(3,2) vSatPos(3,3)] ;
    orbNorm = cross(orbPlaneB-orbPlaneA,orbPlaneC-orbPlaneA);
	clear orbPlaneA orbPlaneB orbPlaneC
    
    % Plot orbit as blue line
    plot3(vSatPos(:,1),vSatPos(:,2),vSatPos(:,3))
    hold on

%% Plot satellite 

%% No sat
	if strcmp(def.satelliteModel, 'none') 
    plot3(satPos(1), satPos(2), satPos(3),'color', 'red', 'Marker', 'o', 'MarkerSize', 7 , 'linewidth',2)

	% 3D sphere when no satellite
    [x, y, z] = ellipsoid(satPos(1),satPos(2),satPos(3), 50, 50, 50, 20);
    ship=surf(x,y,z, 'FaceColor', 'none', 'EdgeColor', 'none');
    
	elseif strcmp(def.satelliteModel, 'stl')
%% STL sat
    points=fv.Points*def.stlfactor+satPos;
	ship = trimesh(fv.ConnectivityList, points(:,1),points(:,2),points(:,3));
	set(ship, 'FaceColor', def.stlBaseColor, 'EdgeColor', def.stlEdgeColor);
	%ship = patch(fv,'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none',        'FaceLighting',    'gouraud', 'AmbientStrength', 0.15);
	end


%% Plot velocity Vec
    %plot velocity vector
    quiver3(satPos(1), satPos(2), satPos(3),velocity(1),velocity(2),velocity(3),'linewidth',2)
     
    %plot inertial kos
    a=6771000;
    axis_length = a; %length of coordinate system axis
    x_axis= [axis_length;0;0];
    y_axis= [0;axis_length;0];
    z_axis= [0;0;axis_length];
    
	%TODO: inertial coordsys
    %quiver3(0,0,0,axis_length*1.5,0,0,'linewidth',3,'color', 'blue')
    %quiver3(0,0,0,0,axis_length*1.5,0,'linewidth',3,'color', 'red')
    %quiver3(0,0,0,0,0,axis_length*1.5,'linewidth',3,'color', 'green')
    
    %% Plot body fixed frame
    q_k = vSatOri(i,:);
    
    %transform inertial kos to body fixed by quaternion multiplication
    %(xrvar qnefgryyhat qrf Beovgny Xbs in eci daher vecbyquattransposed)
    x_axis_body = HPS_transformVecByQuatTransposed(x_axis,q_k);
    y_axis_body = HPS_transformVecByQuatTransposed(y_axis,q_k);
    z_axis_body = HPS_transformVecByQuatTransposed(z_axis,q_k);
    
    
    % Vektor [1*a/8*1000,0,0] turn with quaternion
    %quiver3(state11(i),state12(i),state13(i),x_axis_body(1),x_axis_body(2),x_axis_body(3),'linewidth',2,'color', 'blue')
    %quiver3(state11(i),state12(i),state13(i),y_axis_body(1),y_axis_body(2),y_axis_body(3),'linewidth',2,'color', 'red')
    %quiver3(state11(i),state12(i),state13(i),z_axis_body(1),z_axis_body(2),z_axis_body(3),'linewidth',2,'color', 'green')
    
    
%% Plot earth
    box off
	axis vis3d off
    
    grs80 = referenceEllipsoid('grs80','m');
    
    ax = axesm(	'globe','Geoid',grs80,'Grid','off', ...
				'GLineWidth',1,'GLineStyle','-', ...
        		'Gcolor','black','Galtitude',100);

	% Load Earth image for texture map
    image_file = def.earthfile;
    cdata = imread(image_file);
	% globe transparency level, 1 = opaque, through 0 = invisible

	% rendering earth
    [x, y, z] = ellipsoid(0, 0, 0, a, a, a, def.earthPanels);
    globe = surf(x,y,-z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
    set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', def.earthTransparency, 'EdgeColor', 'none');
    

%% Rotate Earth
    deg_rot = t/3600*15;
    rotate(globe, [0 0 1], deg_rot, [0 0 0]);
   
	%% TODO: enable rotation of ship. 
	%rotate(ship, [1 0 0], 90, satPos);

%% Set Viewpoint 
    campos('manual')
    view(def.viewAngle);
    switch def.camMode 
		case 'fixed'
		case 'follow'
			%TODO: add angle option
            camtarget(satPos)
			camlookat(ship)
			campos(satPos+vfront)
			camup(satPos);
		case 'fov'
			%TODO: add angle option
            camtarget([0 0 0])
            camlookat(ship)
			campos(satView)
			camup(orbNorm)
	end

%% Optimize visuals

	%TODO: Add a camera light, and tone down the specular highlighting
	%camlight('headlight');
	%material('dull');
    %shading faceted; 
    
%% Equalize axes and Framesize
    xlim(xl);
    ylim(yl);
    title(tstr);
    

%% Generate frame from current plot
    if debug == 3
		drawnow;
    else if debug > 0
		drawnow;
		break;
	else
		frameVec(i) = getframe(fig, [0 0 width height]);
	end
    
	if def.progress
		fprintf('\b#\n');
	end
    
end % forloop
end % function
