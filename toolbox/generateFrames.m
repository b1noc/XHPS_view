function frameVec = generateFrames(sat, settings)

runStart = tic;
% getting relative path of toolbox files
[relpath,~,~] = fileparts(mfilename('fullpath'));

%% Setting settings
% Defining default values for unspecified user settings
def = struct(	'debug', 0, ...	
				'nodisplay', 0, ...	
				'step', 1, ...
				'progress', 1, ...
				'resolution', '1080p', ...
				'sim_step', 10, ...
				't_start', 730486, ...
				'camMode', 'earthCentered', ...
				'borderScale', 1, ...
				'viewAngle', [0 0], ...
				'rollAngle', 0, ...
				'zoom', 1, ...
				'showEarth', 1, ...
				'earthTransparency', 1, ...
				'earthPanels', 180, ...
				'showECEF', 1, ...
				'showECI', 1, ...
				'earthVecLength', 1.5, ... 
				'earthfile', strcat(relpath,'/res/earth.jpg'), ...
				'mainSat', 1 ...
			);

defsat = struct('states', 'none', ...
				'satelliteModel', 'none', ...
				'satFactor', 1, ...
				'baseColor', 'yellow', ...
				'edgeColor', 'black', ...
				'showGroundTrack', 1, ...
				'groundTrackLenght', 1, ...
				'showVelocity', 1, ...
				'showSatCoordinates', 1, ...
				'satVecLength', 1, ... 
				'velVecLength', 500, ... 
				'name', '', ...
				'encpath', strcat(relpath,'/res/'), ...
				'stlfile', strcat(relpath,'/res/jeb.stl') ...
			);

% Overwriting default settings with specified user settings
diff = intersect(fieldnames(settings),fieldnames(def));
for n = 1:length(diff)
	def.(diff{n}) = settings.(diff{n});
end

% Overwriting default settings of each satellite with user specified settings
%TODO: simplify
diff = setxor(fieldnames(sat(1)),fieldnames(defsat));
for n = 1:length(diff)
	sat(1).(diff{n}) = defsat.(diff{n});
end
keys = fieldnames(defsat);
for s = 2:length(sat)
	for n = 2:length(keys)
		if isempty(sat(s).(keys{n}))
			sat(s).(keys{n}) = defsat.(keys{n});	
		end
	end
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



%% SATELLITE PREP
% finding the shortest states matrix
lshort = length(sat(1).states);
for s = 2:length(sat)
	lsat = length(sat(s).states);
	if lsat < lshort
		lshort = lsat;
	end
end

%% Reading in states
% define last frame 
endframe=lshort;
ie = ceil(endframe/def.step);

for s = 1:length(sat)
	sat(s).states_full = sat(s).states;
	sat(s).states = sat(s).states(1:def.step:endframe,:);

%% Load 3d model file
	if strcmp(sat(s).satelliteModel,'stl')
		fv = stlread(sat(s).stlfile);
	end

	sat(s).gtlength = ceil(ie * sat(s).groundTrackLenght);
	sat(s).gt=nan(sat(s).gtlength,3);

%% DEBUG-mode
% Display loaded stl file in figure (DEBUG-mode)
	if def.debug == 1
		figure
		hold on

		if isempty(sat(s).name)
			title(['Satellite ' num2str(s,'%1d')])
		else
			title(sat(s).name)
		end
		
		if strcmp(sat(s).satelliteModel, 'stl')
			points=fv.Points*sat(s).satFactor;
			ship = trimesh(fv.ConnectivityList, points(:,1),points(:,2),points(:,3));
			set(ship, 'FaceColor', sat(s).baseColor, 'EdgeColor', sat(s).edgeColor);
			view([0 0])
			axis equal
		elseif strcmp(sat(s).satelliteModel, 'enc')
			nodes = plotSat(sat(s).encpath, [0 0 0 1]);
			for j = 1:4:length(nodes)
				nn = (nodes(j:j+3,:))*sat(s).satFactor;
				fill3(nn(:,1),nn(:,2),nn(:,3), sat(s).baseColor, 'EdgeColor', sat(s).edgeColor, 'LineWidth', 2);
			end
			view(3)
			axis equal
		end
	end

end % satPrep

% Hide figure nodisplay flag
if def.nodisplay == 1
    fig = figure('visible','off');
else
    fig = figure;
end

%% Prepare looprun
r_earth=6771000;
fig.Position(3:4)=[width,height];
[x, y, z] = ellipsoid(0, 0, 0, r_earth, r_earth, r_earth, def.earthPanels);
globe = surf(x,y,-z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
hold on
for i = 1:length(sat)
	plot3(sat(s).states(:,11), sat(s).states(:,12), sat(s).states(:,13));
end

% scaling axes
xl = xlim*def.borderScale;
yl = ylim*def.borderScale;
zl = zlim*def.borderScale;
xlim(xl)
ylim(yl)
zlim(zl)


% initializing vector to save frames
frameVec = struct('cdata', cell(1, ie), 'colormap', cell(1, ie));



%% printing progress bar 
if def.progress && ~def.debug > 0
	fprintf('Progress:\n');
	progLength = fprintf(['\n[' repmat('.',1,60) '] 0%%\n']);
end

%% Main loop
for i=1:ie
	loopStart = tic;
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
	hold on

	if def.showEarth == 1
		set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', def.earthTransparency, 'EdgeColor', 'none');
		
	%% Rotate Earth
		rotate(globe, [0 0 1], rot_ang, [0 0 0]);
	end


%% Calculate satellite and camera position vectors

	for s = 1:length(sat)
		satPos = [sat(s).states(i,11), sat(s).states(i,12), sat(s).states(i,13)];
		velocity = [sat(s).states(i,8), sat(s).states(i,9), sat(s).states(i,10)]*sat(s).velVecLength;
		sn = satPos./norm(satPos);

		if sat(s).showGroundTrack
			rotmat = [cosd(rot_ang) -sind(rot_ang) 0; sind(rot_ang) cosd(rot_ang) 0; 0 0 1];
			
			sat(s).gt = sat(s).gt*inv(rotmat);
			for f = 1:sat(s).gtlength-1
				ff=sat(s).gtlength+1-f;
				sat(s).gt(ff,:) = sat(s).gt(ff-1,:);
			end
			sat(s).gt(1,:) = sn*r_earth+1e5*sn;
			plot3(sat(s).gt(:,1),sat(s).gt(:,2),sat(s).gt(:,3), 'color', 'green', 'linewidth', 2)

			sat(s).gt = sat(s).gt*rotmat;
		end
%% Plot satellite 

		%Plot orbit as blue line
		plot3(sat(s).states_full(:,11),sat(s).states_full(:,12),sat(s).states_full(:,13))
		hold on

		%% No sat
		if strcmp(sat(s).satelliteModel, 'none') 
		plot3(satPos(1), satPos(2), satPos(3),'color', 'red', 'Marker', 'o', 'MarkerSize', 7 , 'linewidth',2)

		% 3D sphere when no satellite
		[x, y, z] = ellipsoid(satPos(1),satPos(2),satPos(3), 50, 50, 50, 20);
		ship=surf(x,y,z, 'FaceColor', 'none', 'EdgeColor', 'none');
		
		%% STL sat
		elseif strcmp(sat(s).satelliteModel, 'stl')
			points = fv.Points*sat(s).satFactor+satPos;
			ship = trimesh(fv.ConnectivityList, points(:,1),points(:,2),points(:,3));
			set(ship, 'FaceColor', sat(s).baseColor, 'EdgeColor', sat(s).edgeColor);

			rotVec = rad2deg(HPS_quat2euler(vSatOri(i,:)));
			rotate(ship, [1 0 0], rotVec(1), satPos);
			rotate(ship, [0 1 0], rotVec(2), satPos);
			rotate(ship, [0 0 1], rotVec(3), satPos);

			%ship = patch(fv,'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none',        'FaceLighting',    'gouraud', 'AmbientStrength', 0.15);

		%% ENC sat
		elseif strcmp(sat(s).satelliteModel, 'enc')
			[x, y, z] = ellipsoid(satPos(1),satPos(2),satPos(3), 50, 50, 50, 20);
			ship=surf(x,y,z, 'FaceColor', 'none', 'EdgeColor', 'none');
			nodes = plotSat(sat(s).encpath, sat(s).states(i,4:7));
			for j = 1:4:length(nodes)
				nn = (nodes(j:j+3,:))*sat(s).satFactor+satPos;
				fill3(nn(:,1),nn(:,2),nn(:,3), sat(s).baseColor, 'EdgeColor', sat(s).edgeColor);
			end
		end


		%% Plot velocity Vec
		if sat(s).showVelocity
			quiver3(satPos(1), satPos(2), satPos(3),velocity(1),velocity(2),velocity(3),'linewidth',2, 'color', [0.4902 0.2353 0.5961])
		end


	
%% Plot satellite coord system

			axis_length = sat(s).satVecLength*r_earth/2; %length of coordinate system axis

			q_k = sat(s).states(i,4:7);
			
			%transform inertial kos to body fixed by quaternion multiplication
			x_axis_bodyUnit = HPS_transformVecByQuatTransposed([1;0;0],q_k);
			y_axis_bodyUnit = HPS_transformVecByQuatTransposed([0;1;0],q_k);
			z_axis_bodyUnit = HPS_transformVecByQuatTransposed([0;0;1],q_k);

			x_axis_body = x_axis_bodyUnit*axis_length;
			y_axis_body = y_axis_bodyUnit*axis_length;
			z_axis_body = z_axis_bodyUnit*axis_length;

			if sat(s).showSatCoordinates
			quiver3(satPos(1), satPos(2), satPos(3), x_axis_body(1), x_axis_body(2),x_axis_body(3), 'linewidth', 2, 'color', [0 0.4471 0.7412])
			quiver3(satPos(1), satPos(2), satPos(3), y_axis_body(1),y_axis_body(2),y_axis_body(3),'linewidth',2,'color', [0.9059 0.2980 0.2353])
			quiver3(satPos(1), satPos(2), satPos(3), z_axis_body(1),z_axis_body(2),z_axis_body(3),'linewidth',2,'color', [0.4667 0.6745 0.1882])
		end
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
	
%% Set Viewpoint 


    satPos = [sat(def.mainSat).states(i,11) sat(def.mainSat).states(i,12) sat(def.mainSat).states(i,13)];
	velocity = [sat(def.mainSat).states(i,8), sat(def.mainSat).states(i,9), sat(def.mainSat).states(i,10)]*sat(def.mainSat).velVecLength;
	sn = satPos./norm(satPos);

	axis_length = sat(def.mainSat).satVecLength*r_earth/2; %length of coordinate system axis

	q_k = sat(def.mainSat).states(i,4:7);
			
	%transform inertial kos to body fixed by quaternion multiplication
	x_axis_bodyUnit = HPS_transformVecByQuatTransposed([1;0;0],q_k);
	y_axis_bodyUnit = HPS_transformVecByQuatTransposed([0;1;0],q_k);
	z_axis_bodyUnit = HPS_transformVecByQuatTransposed([0;0;1],q_k);


    switch def.camMode 
		case 'earthCentered'
			%turning off axis
			set(gca, 'Visible', 'off')
			camtarget([0 0 0])
			%campos([0 -r_earth*def.zoom*10 0])
			campos([0 -1*(1/def.zoom)*1e8 0])
			camorbit(90+def.viewAngle(1),def.viewAngle(2))
            camroll(def.rollAngle)
		case 'satCentered'
			%turning off axis
			axis image vis3d off
			camtarget(satPos)
			campos(satPos+sn*(1/def.zoom)*1e8)
			camup(velocity)
			camroll(def.rollAngle)
			camorbit(def.viewAngle(1),def.viewAngle(2),'camera')
		case 'satFixed'
			%turning off axis
			axis image vis3d off
			camtarget(satPos)
			campos((satPos-z_axis_bodyUnit'*(1/def.zoom)*1e8))
			camup(x_axis_bodyUnit')
			camroll(def.rollAngle)
			camorbit(def.viewAngle(1),def.viewAngle(2),'camera')
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

	if def.progress && ~def.debug > 0
		progress = floor(i*100/ie);
		prog = floor(progress * .6);
		inve= ceil(100*.6 - prog);
		fprintf(repmat('\b',1,progLength));
		progLength = fprintf('\n[%s%s] %d%%', repmat('#',1,prog), repmat('.',1,inve), progress);
	end

%% Generate frame from current plot
    if def.debug == 1
		loopStop = toc(loopStart);
		fprintf('One step ran in %.3f seconds', loopStop)
		drawnow;
		break;
	elseif def.debug == 2
		loopStop(i) = toc(loopStart);
		drawnow;
	else
		loopStop(i) = toc(loopStart);
		frameVec(i) = getframe(fig, [0 0 width height]);
    end
end % forloop
	fprintf('\nIt plotted:\n%03d steps\n%07.3fs total duration\n%07.3fs average step duration\n', length(loopStop), toc(runStart), mean(loopStop))
fprintf('\n')
end % function
