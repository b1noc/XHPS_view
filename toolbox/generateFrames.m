function frameVec = generateFrames(states, settings, debug)

% getting relative path of toolbox files
[relpath,~,~] = fileparts(mfilename('fullpath'));

%% Setting settings
% Defining default values for unspecified user settings
def = struct(	'step', 1, ...
				'resolution', '1080p', ...
				'sim_step', 10, ...
				't_start', 730486, ...
				'borderScale', 1.2, ...
				'satVecLength', 1, ... 
				'velVecLength', 500, ... 
				'earthVecLength', 1.5, ... 
				'camMode', 'fixed', ...
				'zoom', 1, ...
				'earthTransparency', 1, ...
				'satelliteModel', 'none', ...
				'viewAngle', [0 90], ...
				'progress', 1, ...
				'satFactor', 1, ...
				'stlfile', strcat(relpath,'/res/jeb.stl'), ...
				'encpath', strcat(relpath,'/res/'), ...
				'baseColor', 'yellow', ...
				'edgeColor', 'black', ...
				'showEarth', 1, ...
				'earthfile', strcat(relpath,'/res/earth.jpg'), ...
				'earthPanels', 20, ...
				'showVelocity', 1, ...
				'showSatCoordinates', 1, ...
				'showECEF', 1, ...
				'showECI', 1, ...
				'cRoll', 0, ...
				'cOrb', [0,0] ...
			);

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
vSatPosFull = states(:,11:13);
vSatPos = state(:,11:13);
vSatOri = state(:,4:7);
vSatVel = state(:,8:10);

%% Load 3d model file
if strcmp(def.satelliteModel,'stl')
	fv = stlread(def.stlfile);
end

%% Defining vido quality based on user settings
switch def.resolution
	case '720p'
		width = 1280;
		height = 720;
	case '1080p'
		width = 1920;
		height = 1080;
    case '4k'
        width = 3840;
        height = 2160;
end

%% DEBUG-mode
% Display loaded stl file in figure (DEBUG-mode)
if debug == 1
	figure
	hold on
	if strcmp(def.satelliteModel, 'stl')
		points=fv.Points*def.satFactor;
		ship = trimesh(fv.ConnectivityList, points(:,1),points(:,2),points(:,3));
		set(ship, 'FaceColor', def.baseColor, 'EdgeColor', def.edgeColor);
		view([0 0])
		axis equal
	elseif strcmp(def.satelliteModel, 'enc')
		nodes = plotSat(def.encpath, [0 0 0 1]);
		for j = 1:4:length(nodes)
			nn = (nodes(j:j+3,:))*def.satFactor;
			fill3(nn(:,1),nn(:,2),nn(:,3), def.baseColor, 'EdgeColor', def.edgeColor, 'LineWidth', 2);
		end
		view(3)
		axis equal
	end
end

% Display frame in DEBUG-mode
if debug > 0
    fig = figure;
else
    fig = figure('visible','off');
end

%% Prepare looprun
r_earth=6771000;
fig.Position(3:4)=[width,height];
[x, y, z] = ellipsoid(0, 0, 0, r_earth, r_earth, r_earth, def.earthPanels);
globe = surf(x,y,-z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
hold on
plot3(vSatPos(:,1),vSatPos(:,2),vSatPos(:,3))

% scaling axes
xl = xlim*def.borderScale;
yl = ylim*def.borderScale;
zl = zlim*def.borderScale;
xlim(xl)
ylim(yl)
zlim(zl)

% define last frame 
ie = ceil(endframe/def.step);

% initializing vector to save frames
frameVec = struct('cdata', cell(1, ie), 'colormap', cell(1, ie));


%% printing progress bar 
if def.progress && ~debug > 0
	fprintf('Progress:\n');
	progLength = fprintf(['\n[' repmat('.',1,60) '] 0%%\n']);
end

%% Main loop
for i=1:ie
    clf
    
%% Calculate time and rotation
	% Calculating index k from original states matrix
    k=i*def.step - def.step;
	% Calculating timestamp of current state by adding x*sim_step to the starttime
    tstr = datestr(def.t_start + k * seconds(def.sim_step));
    t = def.t_start + seconds(k * def.sim_step);

	% Calculating earths rotation angle based on time
	dt = t-datetime('2000-01-01 12:00:00') ;
	rot_ang = 360.9856123035484 * days(dt) + 280.46; % [deg]

%% Calculate satellite and camera position vectors
	% velocity vector
    velocity = [vSatVel(i,1)*def.velVecLength,vSatVel(i,2)*def.velVecLength,vSatVel(i,3)*def.velVecLength];
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
    plot3(vSatPosFull(:,1),vSatPosFull(:,2),vSatPosFull(:,3))
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
    points=fv.Points*def.satFactor+satPos;
	ship = trimesh(fv.ConnectivityList, points(:,1),points(:,2),points(:,3));
	set(ship, 'FaceColor', def.baseColor, 'EdgeColor', def.edgeColor);
	%ship = patch(fv,'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none',        'FaceLighting',    'gouraud', 'AmbientStrength', 0.15);
%% ENC sat
	elseif strcmp(def.satelliteModel, 'enc')
		[x, y, z] = ellipsoid(satPos(1),satPos(2),satPos(3), 50, 50, 50, 20);
		ship=surf(x,y,z, 'FaceColor', 'none', 'EdgeColor', 'none');
		nodes = plotSat(def.encpath, vSatOri(i,:));
		for j = 1:4:length(nodes)
			nn = (nodes(j:j+3,:))*def.satFactor+satPos;
			fill3(nn(:,1),nn(:,2),nn(:,3), def.baseColor, 'EdgeColor', def.edgeColor);
		end
	end


%% Plot velocity Vec
    %plot velocity vector
	if def.showVelocity
		quiver3(satPos(1), satPos(2), satPos(3),velocity(1),velocity(2),velocity(3),'linewidth',2, 'color', [0.4902 0.2353 0.5961])
	end
     
%% Plot ECI


	if def.showECI
		axis_length = def.earthVecLength*r_earth; %length of coordinate system axis
		quiver3(0,0,0,axis_length,0,0,'linewidth',3,'color', [0.4980 0.7020 0.8353])
		quiver3(0,0,0,0,axis_length,0,'linewidth',3,'color', [0.9451 0.5804 0.5412])
		quiver3(0,0,0,0,0,axis_length,'linewidth',3,'color', [0.5098 0.8784 0.6667])
	end

%% Plot ECEF
	if def.showECEF
		axis_length = def.earthVecLength*r_earth; %length of coordinate system axis
		ecef = HPS_computeDCMFromRotationAngle( deg2rad(rot_ang), 3 ) * [axis_length, 0, 0; 0, axis_length, 0; 0,0,axis_length];
		quiver3(0,0,0,ecef(1,1),ecef(1,2),ecef(1,3),'linewidth',3,'color', [0 0.4471 0.7412])
		quiver3(0,0,0,ecef(2,1),ecef(2,2),ecef(2,3),'linewidth',3,'color', [0.9059 0.2980 0.2353])
		quiver3(0,0,0,ecef(3,1),ecef(3,2),ecef(3,3),'linewidth',3,'color', [0.4667 0.6745 0.1882])
	end
	
%% Plot satellite coord system

		axis_length = def.satVecLength*r_earth/2; %length of coordinate system axis

		q_k = vSatOri(i,:);
		
		%transform inertial kos to body fixed by quaternion multiplication
		x_axis_bodyUnit = HPS_transformVecByQuatTransposed([1;0;0],q_k);
		y_axis_bodyUnit = HPS_transformVecByQuatTransposed([0;1;0],q_k);
		z_axis_bodyUnit = HPS_transformVecByQuatTransposed([0;0;1],q_k);

		x_axis_body = x_axis_bodyUnit*axis_length;
		y_axis_body = y_axis_bodyUnit*axis_length;
		z_axis_body = z_axis_bodyUnit*axis_length;

		if def.showSatCoordinates
		quiver3(vSatPos(i,1),vSatPos(i,2),vSatPos(i,3),x_axis_body(1),x_axis_body(2),x_axis_body(3),'linewidth',2,'color', [0 0.4471 0.7412])
		quiver3(vSatPos(i,1),vSatPos(i,2),vSatPos(i,3),y_axis_body(1),y_axis_body(2),y_axis_body(3),'linewidth',2,'color', [0.9059 0.2980 0.2353])
		quiver3(vSatPos(i,1),vSatPos(i,2),vSatPos(i,3),z_axis_body(1),z_axis_body(2),z_axis_body(3),'linewidth',2,'color', [0.4667 0.6745 0.1882])
	end
    
%% Plot earth
    
    grs80 = referenceEllipsoid('grs80','m');
    
    ax = axesm(	'globe','Geoid',grs80,'Grid','off', ...
				'GLineWidth',1,'GLineStyle','-', ...
        		'Gcolor','black','Galtitude',100);

	% Load Earth image for texture map
    image_file = def.earthfile;
    cdata = imread(image_file);

	
%% rendering earth
		[x, y, z] = ellipsoid(0, 0, 0, r_earth, r_earth, r_earth, def.earthPanels);
		globe = surf(x,y,-z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);

	if def.showEarth == 1
		set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', def.earthTransparency, 'EdgeColor', 'none');
		

	%% Rotate Earth
		rotate(globe, [0 0 1], rot_ang, [0 0 0]);
	end
   
	%% remove axis
    box off
	axis vis3d off

%% Set Viewpoint 
    campos('manual')
    %view(def.viewAngle);
    switch def.camMode 
		case 'fixed'
		case 'follow'
            camtarget(satPos)
			camlookat(ship)
			campos(satPos+vfront)
			camup(satPos);
		case 'fov'
            camtarget([0 0 0])
            camlookat(ship)
			campos(-satView)
			camzoom(1)
			camup(orbNorm)
		case 'pilot'
			front = satPos+x_axis_body';
			camtarget(front)
			camlookat(ship)
			campos(satPos-x_axis_body'*def.zoom/-5e6)
			camup(z_axis_body');
		case 'earthCentered'
			camtarget([0 0 0])
			%campos([0 -r_earth*def.zoom*10 0])
			campos([0 -1*(1/def.zoom)*1e8 0])
			camorbit(90+def.cOrb(1),def.cOrb(2))
            camroll(def.cRoll)
		case 'satCentered'
			camtarget(satPos)
			campos(satPos+sn*(1/def.zoom)*1e8)
			camup(velocity)
			camroll(def.cRoll)
			camorbit(def.cOrb(1),def.cOrb(2),'camera')
		case 'satFixed'
			camtarget(satPos)
			campos((satPos-z_axis_bodyUnit'*(1/def.zoom)*1e8))
			camup(x_axis_bodyUnit')
			camroll(def.cRoll)
			camorbit(def.cOrb(1),def.cOrb(2),'camera')

			%{ 
			# Other modes 
				- Earthfixed (move with earth rotation)
			## defaults
				- orbit planar
				- 
			%}
    end
       

%% Equalize axes and Framesize
    xlim(xl);
    ylim(yl);
    zlim(zl);

	wh = get(fig, 'Position');
	relWidth = 1/wh(3)*180;
	relHeight = 1/wh(4)*25;
	paddingX = 1/wh(3)*10;
	paddingY = 1/wh(4)*10;
	annotation('textbox', [paddingX, paddingY, relWidth, relHeight], 'EdgeColor', 'black', 'BackgroundColor', 'white', 'string', tstr);

	if def.progress && ~debug > 0
		progress = floor(i*100/ie);
		prog = floor(progress * .6);
		inv = ceil(100*.6 - prog);
		fprintf(repmat('\b',1,progLength));
		progLength = fprintf('\n[%s%s] %d%%', repmat('#',1,prog), repmat('.',1,inv), progress);
	end

%% Generate frame from current plot
    if debug == 2
		drawnow;
    elseif debug > 0
		drawnow;
		break;
	else
		frameVec(i) = getframe(fig, [0 0 width height]);
    end
    
    
end % forloop
fprintf('\n')
end % function
