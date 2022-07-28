function frameVec = generateFrames(sat, settings)

% start timer
runStart = tic; 

% getting relative path of toolbox files
[relpath,~,~] = fileparts(mfilename('fullpath'));

% Default values of unspecified user settings
def = struct(	'debug', 0, ...	
				'nodisplay', 0, ...	
				'progress', 1, ...
				'step', 1, ...
				'fps', 1, ...
				'duration', 10, ...
				'sim_step', 10, ...
				't_start', 730486, ...
				'camMode', 'earthCentered', ...
				'mainSat', 1, ...
				'resolution', '1080p', ...
				'viewAngle', [0 0], ...
				'rollAngle', 0, ...
				'zoom', 1, ...
				'earth', 'on', ...
				'earthPanels', 180, ...
				'earthTransparency', 1, ...
				'ecef', 0, ...
				'eci', 0, ...
				'earthVecLength', 1.5, ... 
				'earthfile', strcat(relpath,'/res/earth.jpg') ...
			);

defsat = struct('states', 'none', ...
				'satModel', 'none', ...
				'name', '', ...
				'satScale', 1, ...
				'color', 'yellow', ...
				'edgeColor', 'black', ...
				'plotOrbit', 1, ...
				'orbitColor', '#4DBEEE', ...
				'groundTrack', 0, ...
				'groundTrackLenght', 1, ...
				'satCoordinates', 1, ...
				'velocityVec', 1, ...
				'axisLength', 1, ... 
				'velVecLength', 500, ... 
				'encPath', strcat(relpath,'/res/'), ...
				'stlFile', strcat(relpath,'/res/jeb.stl') ...
			);

% Overwriting default settings with specified user settings
diff = intersect(fieldnames(settings),fieldnames(def));
for n = 1:length(diff)
	def.(diff{n}) = settings.(diff{n});
end

% Creating not specified keys with default keys 
diff = setxor(fieldnames(sat(1)),fieldnames(defsat));
for n = 1:length(diff)
	sat(1).(diff{n}) = defsat.(diff{n});
end

% Overwrite empty sat values with default values for each satellite
keys = fieldnames(defsat);
for s = 1:length(sat)
	for n = 1:length(keys)
		if isempty(sat(s).(keys{n}))
			sat(s).(keys{n}) = defsat.(keys{n});	
		end
	end
end

% Open figure (hidden when nodisplay flag is set)
if def.nodisplay == 1
    fig = figure('visible','off');
else
    fig = figure;
end
hold on


% Defining vido size based on selected quality
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

% setting windowsize
fig.Position(3:4)=[width,height];
axis equal vis3d off

% define last frame by shortest matrix
shortestStates =  sat.states;
endframe = length(shortestStates);

def.step = ceil(endframe/(def.fps * def.duration));
if def.step < 1
	fpsnew = endframe/def.duration;
	fprintf('Warning: not enough datapoints for %d fps\n Reduced framrate to %d fpsnew', def.fps, fpsnew);
	def.step = 1;
end

ie = ceil(endframe/def.step);

% initializing vector to save frames
frameVec = struct('cdata', cell(1, ie), 'colormap', cell(1, ie));

% printing progress bar 
if def.progress && ~def.debug > 0
	fprintf('Progress:\n');
	progLength = fprintf(['\n[' repmat('.',1,60) '] 0%%\n']);
end

% Plot earth
r_earth=6771000;
grs80 = referenceEllipsoid('grs80','m');

if ~strcmp(def.earth, 'off') && ~def.earth == 0
	axesm(	'globe','Geoid',grs80,'Grid','on', ...
			'GLineWidth',1,'GLineStyle','-', ...
			'Gcolor','black','Galtitude',100);
	[x, y, z] = ellipsoid(0, 0, 0, r_earth, r_earth, r_earth, ...
						  def.earthPanels);

	globe = surf(x,y,-z, 'FaceColor', 'none', 'EdgeColor', .5*[1 1 1]);
	globeX = globe.XData;
	globeY = globe.YData;
	globeZ = globe.ZData;

	% Render earth and disable grid
	if ~strcmp(def.earth, 'grid') 
		image_file = def.earthfile;
		cdata = imread(image_file);
		set(globe, 'FaceColor', 'texturemap', 'CData', cdata, ...
			'FaceAlpha', def.earthTransparency, 'EdgeColor', 'none');
	end
end


% Plot ECI
earthAxisLength = def.earthVecLength*r_earth; 
if def.eci
	quiver3(0, 0, 0, earthAxisLength, 0, 0, ...
			'linewidth',3, 'color', [0.4980 0.7020 0.8353])
	quiver3(0, 0, 0, 0, earthAxisLength, 0, ...
			'linewidth', 3, 'color', [0.9451 0.5804 0.5412])
	quiver3(0, 0, 0, 0, 0, earthAxisLength, ...
			'linewidth', 3, 'color', [0.5098 0.8784 0.6667])
end

% Plot ECEF
if def.ecef
	ecefX = quiver3(0, 0, 0, earthAxisLength, 0, 0, ...
			'linewidth', 3, 'color', [0 0.4471 0.7412]);
	ecefY = quiver3(0, 0, 0, 0, earthAxisLength, 0, ...
			'linewidth', 3, 'color', [0.9059 0.2980 0.2353]);
	ecefZ = quiver3(0, 0, 0, 0, 0, earthAxisLength, ...
			'linewidth', 3, 'color', [0.4667 0.6745 0.1882]);
end

% Plot Satellite
for s = 1:length(sat)
	sat(s).states_full = sat(s).states;
	sat(s).states = sat(s).states(1:def.step:endframe,:);

	% Load 3d model file
	if strcmp(sat(s).satModel,'stl')
		fv = stlread(sat(s).stlfile);
	end

	% calculate groundtrack length
	sat(s).gtlength = ceil(ie * sat(s).groundTrackLenght);
	sat(s).gt=nan(sat(s).gtlength,3);

	% Plot orbit
	if sat(s).plotOrbit
		sat(s).plot_orbit = plot3(sat(s).states_full(:,11), ...
								  sat(s).states_full(:,12), ...
								  sat(s).states_full(:,13), ...
								  'color', sat(s).orbitColor, 'LineWidth',.1);
	end
								  %'color', '#D95319', 'LineWidth',.1);

	if strcmp(sat(s).satModel, 'none') 
		sat(s).plot_sat = plot3(0, 0, 0, 'color', 'red', ...
										 'Marker', 'o', ...
										 'MarkerSize', 7 , ...
										 'linewidth',2);

	elseif strcmp(sat(s).satModel, 'enc')
		et = load([sat(s).encPath, 'et.txt.']);
		nt = load([sat(s).encPath, 'nt.txt.']);
		faces_matrix = [et(:,6) et(:,7) et(:,8) et(:,9)];
		sat(s).plot_sat = patch('Vertices', nt, 'Faces', faces_matrix);
		sat(s).plot_sat.FaceColor = sat(s).color;
		sat(s).x = sat(s).plot_sat.XData;
		sat(s).y = sat(s).plot_sat.YData;
		sat(s).z = sat(s).plot_sat.ZData;
	end

	if sat(s).satCoordinates
		sat(s).x_axis = quiver3(0, 0, 0, 1, 0, 0, 'linewidth', 2, ...
								'color', [0 0.4471 0.7412]);
		sat(s).y_axis = quiver3(0, 0 ,0, 0, 1, 0, 'linewidth', 2, ...
								'color', [0.9059 0.2980 0.2353]);
		sat(s).z_axis = quiver3(0, 0, 0, 0, 0, 1, 'linewidth', 2, ...
								'color', [0.4667 0.6745 0.1882]);
	end
	if sat(s).velocityVec
		sat(s).vel_vector = quiver3(0, 0, 0, 1, 0, 0, 'linewidth',2, ...
									'color', [0.4902 0.2353 0.5961]);
	end

	if sat(s).groundTrack
		sat(s).ground_track= [0 0 0];
		sat(s).ground_plot = plot3( 0, 0, 0, 'color', 'green', 'linewidth', 1);
	end
end


%% Main loop
for i=1:ie
	loopStart = tic;
    
	% Calculating timestamp of current state 
    k=i*def.step - def.step;
    tstr = datestr(def.t_start + k * seconds(def.sim_step));
    t = def.t_start + seconds(k * def.sim_step);

	% Calculating earths rotation angle based on time
	dt = t-datetime('2000-01-01 12:00:00') ;
	rot_ang = 360.9856123035484 * days(dt) + 280.46; % [deg]

	%% Rotate Earth
	if ~strcmp(def.earth, 'off') && ~def.earth == 0
		globe.XData = globeX;
		globe.YData = globeY;
		globe.ZData = globeZ;
		rotate(globe, [0 0 1], rot_ang, [0 0 0]);
	end

	% Rotate ECEF coordinate Frame
	if def.ecef
		ecefRot = HPS_computeDCMFromRotationAngle( deg2rad(rot_ang), 3 ) * ...
												   [earthAxisLength, 0, 0; ...
													0, earthAxisLength, 0; ...
													0, 0, earthAxisLength];
		ecefX.UData = ecefRot(1,1);
		ecefX.VData = ecefRot(1,2);
		ecefX.WData = ecefRot(1,3);
		ecefY.UData = ecefRot(2,1);
		ecefY.VData = ecefRot(2,2);
		ecefY.WData = ecefRot(2,3);
	end

	% Position sat
	for s = 1:length(sat)
		satPos = [sat(s).states(i,11), ...
				  sat(s).states(i,12), ...
				  sat(s).states(i,13)];

		if strcmp(sat(s).satModel, 'none') 
			sat(s).plot_sat.XData =satPos(1);
			sat(s).plot_sat.YData =satPos(2);
			sat(s).plot_sat.ZData =satPos(3);

		elseif strcmp(sat(s).satModel, 'enc') 
			r_com = [1.5615 0.9720 -0.3310];

			set(sat(s).plot_sat, ...
				'XData', (sat(s).x-r_com(1))*sat(s).satScale+satPos(1), ...
				'YData', (sat(s).y-r_com(2))*sat(s).satScale+satPos(2), ...
				'ZData', (sat(s).z-r_com(3))*sat(s).satScale+satPos(3)  ...
			);

			rotVec = rad2deg(HPS_quat2euler(sat(s).states(i,4:7)));
			rotate(sat(s).plot_sat, [1 0 0], rotVec(1), satPos);
			rotate(sat(s).plot_sat, [0 1 0], rotVec(2), satPos);
			rotate(sat(s).plot_sat, [0 0 1], rotVec(3), satPos);
		end

		
		q_k = sat(s).states(i,4:7);

		%transform inertial kos to body fixed by quaternion multiplication
		x_axis_bodyUnit = HPS_transformVecByQuatTransposed([1;0;0],q_k);
		y_axis_bodyUnit = HPS_transformVecByQuatTransposed([0;1;0],q_k);
		z_axis_bodyUnit = HPS_transformVecByQuatTransposed([0;0;1],q_k);

		% Move axis vectors
		if sat(s).satCoordinates
			axlen = sat(s).axisLength*sat(s).satScale*10;
			sat(s).x_axis.XData = satPos(1);
			sat(s).x_axis.YData = satPos(2);
			sat(s).x_axis.ZData = satPos(3);
			sat(s).x_axis.UData = x_axis_bodyUnit(1)*axlen;
			sat(s).x_axis.VData = x_axis_bodyUnit(2)*axlen;
			sat(s).x_axis.WData = x_axis_bodyUnit(3)*axlen;
			sat(s).y_axis.XData = satPos(1);
			sat(s).y_axis.YData = satPos(2);
			sat(s).y_axis.ZData = satPos(3);
			sat(s).y_axis.UData = y_axis_bodyUnit(1)*axlen;
			sat(s).y_axis.VData = y_axis_bodyUnit(2)*axlen;
			sat(s).y_axis.WData = y_axis_bodyUnit(3)*axlen;
			sat(s).z_axis.XData = satPos(1);
			sat(s).z_axis.YData = satPos(2);
			sat(s).z_axis.ZData = satPos(3);
			sat(s).z_axis.UData = z_axis_bodyUnit(1)*axlen;
			sat(s).z_axis.VData = z_axis_bodyUnit(2)*axlen;
			sat(s).z_axis.WData = z_axis_bodyUnit(3)*axlen;
		end

		%% Move velocity vector
		if sat(s).velocityVec
			velocity = [sat(s).states(i,8), ...
						sat(s).states(i,9), ...
					   	sat(s).states(i,10)];
			axlen = sat(s).velVecLength*sat(s).satScale/5e2;
			sat(s).vel_vector.XData = satPos(1);
			sat(s).vel_vector.YData = satPos(2);
			sat(s).vel_vector.ZData = satPos(3);
			sat(s).vel_vector.UData = velocity(1)*axlen;
			sat(s).vel_vector.VData = velocity(2)*axlen;
			sat(s).vel_vector.WData = velocity(3)*axlen;
		end

		if sat(s).groundTrack
			sn = satPos./norm(satPos);
			rotmat = [cosd(rot_ang) -sind(rot_ang) 0; ...
					  sind(rot_ang) cosd(rot_ang) 0; ...
					  0 0 1];

			ground_track = sat(s).ground_track*inv(rotmat);
			ground_track = [(r_earth+1e5)*sn; ground_track];

			gt_len = size(ground_track,1)-1;
			sat(s).ground_plot.XData = ground_track(1:gt_len,1);
			sat(s).ground_plot.YData = ground_track(1:gt_len,2);
			sat(s).ground_plot.ZData = ground_track(1:gt_len,3);

			sat(s).ground_track = ground_track*rotmat;
		end
	end


	% Get positional data of mainSat
    satPos = [sat(def.mainSat).states(i,11) ...
			  sat(def.mainSat).states(i,12) ...
			  sat(def.mainSat).states(i,13)];

	velocity = [sat(def.mainSat).states(i,8), ...
				sat(def.mainSat).states(i,9), ...
				sat(def.mainSat).states(i,10)];
	sn = satPos./norm(satPos);
			
	%transform inertial kos to body fixed by quaternion multiplication
	q_k = sat(def.mainSat).states(i,4:7);
	x_axis_bodyUnit = HPS_transformVecByQuatTransposed([1;0;0],q_k);
	z_axis_bodyUnit = HPS_transformVecByQuatTransposed([0;0;1],q_k);

	% Set Viewmode 
	view([0 0])
    switch def.camMode 
		case 'earthCentered'
			camtarget([0 0 0])
			campos([0 -1*(1/def.zoom)*2e8 0])
			camorbit(90+def.viewAngle(1),def.viewAngle(2))
            camroll(def.rollAngle)
		case 'satCentered'
			camtarget(satPos)
			campos(satPos+sn*(1/def.zoom)*2e8)
			camup(velocity)
			camroll(def.rollAngle)
			camorbit(def.viewAngle(1),def.viewAngle(2),'camera')
		case 'satFixed'
			camtarget(satPos)
			campos((satPos-z_axis_bodyUnit'*(1/def.zoom)*2e8))
			camup(x_axis_bodyUnit')
			camroll(def.rollAngle)
			camorbit(def.viewAngle(1),def.viewAngle(2),'camera')
    end
       

	% Calculate position of Textbox
	wh = get(fig, 'Position');
	relWidth = 1/wh(3)*180;
	relHeight = 1/wh(4)*25;
	paddingX = 1/wh(3)*10;
	paddingY = 1/wh(4)*10;
	annotation('textbox', [paddingX, paddingY, relWidth, relHeight], ...
			   'EdgeColor', 'black', 'BackgroundColor', 'white', 'string', tstr);

	% Plot progressbar
	if def.progress && ~def.debug > 0
		progress = floor(i*100/ie);
		prog = floor(progress * .6);
		inve= ceil(100*.6 - prog);
		fprintf(repmat('\b',1,progLength));
		progLength = fprintf('\n[%s%s] %d%%', repmat('#',1,prog), repmat('.',1,inve), progress);
	end

	% Generate frame from current plot
    if def.debug == 1
		drawnow;
		loopStop = toc(loopStart);
		fprintf('One step ran in %.3f seconds', loopStop)
		frameVec = '';
		break;
	elseif def.debug == 2
		drawnow;
		loopStop(i) = toc(loopStart);
	else
		frameVec(i) = getframe(fig, [0 0 width height]);
		loopStop(i) = toc(loopStart);
    end
end % forloop

	% Print statistics
	fprintf('\nPlotted in:\n')
	fprintf('% 7d steps\n', length(loopStop))
	fprintf('% 7.3f total duration [s]\n', toc(runStart))
	fprintf('% 7.3f average step duration [s]\n\n', mean(loopStop))

end % function

