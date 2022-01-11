function aniVec = generateFrames(states, step, sim_step, t_start, angle, res, debug)

if isempty(gcp('nocreate')) && debug ~= 1
   parpool;
end

endframe=length(states);

% size of space around plot TODO: get from Vector lenght
scale = 1.2;

state = states(1:step:endframe,:);
state11= state(:,11);
state12 = state(:,12);
state13 = state(:,13);
state47 = state(:,4:7);
state810 = state(:,8:10);

% Vid quality
switch res
	case '720p'
		width = 1280;
		height = 720;
	case '1080p'
		width = 1920;
		height = 1080;
end

if debug == 1
    fig = figure;
else
    fig = figure('visible','off');
end

plot3(state11,state12,state13);
set(fig, 'Position',  [100, 100, width, height]);

xl = xlim*scale;
yl = ylim*scale;
xlim(xl)
ylim(yl)
ie = floor(endframe/step)+1;
if ie > length(states)
	ie = length(states);
end


m=ie;
fprintf('Progress:\n');
fprintf(['\n' repmat('.',1,m) '\n\n']);

%aniVec(ie)=0;
aniVec = struct('cdata', cell(1, ie), 'colormap', cell(1, ie));

for i=1:ie
    clf
    
    k=i*step - step+1;
    tstr = datestr(t_start + (k-1) * seconds(sim_step));
    t = t_start + (k-1) * sim_step;
       
    % Plot orbit line
    plot3(state11,state12,state13)
    hold on
    
    %plot starting (initial) Point
    %plot3(state11(i),state12(i),state13(i),'color', 'red', 'Marker', 'o', 'MarkerSize', 7 , 'linewidth',2)
    [x, y, z] = ellipsoid(state11(i), state12(i), state13(i), 50, 50, 50, 20);
    shuttle=surf(x,y,-z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
    %set(shuttle, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
		

    
    %plot velocity vector
    quiver3(state11(i),state12(i),state13(i),state810(i,1)*500,state810(i,2)*500,state810(i,3)*500,'linewidth',2)
     
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
    q_k = state47(i,:);
    
    %transform inertial kos to body fixed by quaternion multiplication
    %(xrvar qnefgryyhat qrf Beovgny Xbs in eci daher vecbyquattransposed)
    x_axis_body = HPS_transformVecByQuatTransposed(x_axis,q_k);
    y_axis_body = HPS_transformVecByQuatTransposed(y_axis,q_k);
    z_axis_body = HPS_transformVecByQuatTransposed(z_axis,q_k);
    
    
    % Vektor [1*a/8*1000,0,0] turn with quaternion
    %quiver3(state11(i),state12(i),state13(i),x_axis_body(1),x_axis_body(2),x_axis_body(3),'linewidth',2,'color', 'blue')
    %quiver3(state11(i),state12(i),state13(i),y_axis_body(1),y_axis_body(2),y_axis_body(3),'linewidth',2,'color', 'red')
    %quiver3(state11(i),state12(i),state13(i),z_axis_body(1),z_axis_body(2),z_axis_body(3),'linewidth',2,'color', 'green')
    
    
    %% plot earth
    box off
		axis vis3d off
    
    grs80 = referenceEllipsoid('grs80','m');
    
    ax = axesm('globe','Geoid',grs80,'Grid','off', ...
    	  'GLineWidth',1,'GLineStyle','-',...
        'Gcolor','black','Galtitude',100);

    image_file = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';
% Load Earth image for texture map
    cdata = imread(image_file);
    %cdata = load('Earth_fig_cdata');
    alpha   = .6; % globe transparency level, 1 = opaque, through 0 = invisible
    erad    = a ;%6371008.7714; % equatorial radius (meters)
    prad    = a; %6371008.7714; % polar radius (meters)
    npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels

    [x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
    globe = surf(x,y,-z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
    set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');

    direction = [0 0 1];
    origin = [0 0 0]; 
    deg_rot = t/3600*15;

    rotate(globe, direction, deg_rot, origin);
   
    % velocity vector
    velocity = [state810(i,1)*500,state810(i,2)*500,state810(i,3)*500];
    % position of satelite
    satPos = [state11(i) state12(i) state13(i)];
    satView = satPos./norm(satPos);

		orbPlaneA = [state11(1) state12(1) state13(1)] ;
		orbPlaneB = [state11(2) state12(2) state13(2)] ;
		orbPlaneC = [state11(3) state12(3) state13(3)] ;
		orbNorm = cross(orbPlaneB-orbPlaneA,orbPlaneC-orbPlaneA);


    campos('manual')
    view(angle);
    camlookat(shuttle)
    
    campos(satView*1e8*.9)
    %camtarget(satPos+velocity)
    camtarget([0 0 0])
    camup([orbNorm(1) orbNorm(2) orbNorm(3)])
    
    xlim(xl);
    ylim(yl);
    xlabel(tstr)
    
    if debug == 1
			drawnow;
		else
			aniVec(i) = getframe(fig, [0 0 width height]);
		end
    
    fprintf('\b#\n');
    
end

	if debug ~= 1
		save('aniVec','aniVec', '-v7.3')
	end
end
