function aniVec = generateFrames(states, step, sim_step, t_start, angle, res)

if isempty(gcp('nocreate'))
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

fig = figure('visible','off');
%fig = figure;
plot3(state11,state12,state13)

set(fig, 'Position',  [100, 100, width, height])

x = xlim*scale;
y = ylim*scale;
xlim(x)
ylim(y)
ie = floor(endframe/step)+1;
if ie > length(states)
	ie = length(states)
end


m=ie;
fprintf('Progress:\n');
fprintf(['\n' repmat('.',1,m) '\n\n']);


parfor i=1:ie
    clf
    %D = parallel.pool.DataQueue
    
    
    k=i*step - step+1;
    t = t_start + (k-1) * seconds(sim_step);
    
    % Plot orbit line
    plot3(state11,state12,state13)
    hold on
    
    %plot starting (initial) Point
%    plot3(state11(i),state12(i),state13(i),'color', 'red', 'Marker', 'o', 'MarkerSize', 7 , 'linewidth',2)
    
    %plot initial velocity vector
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
    quiver3(state11(i),state12(i),state13(i),x_axis_body(1),x_axis_body(2),x_axis_body(3),'linewidth',2,'color', 'blue')
    quiver3(state11(i),state12(i),state13(i),y_axis_body(1),y_axis_body(2),y_axis_body(3),'linewidth',2,'color', 'red')
    quiver3(state11(i),state12(i),state13(i),z_axis_body(1),z_axis_body(2),z_axis_body(3),'linewidth',2,'color', 'green')
    
    
    %% plot earth
    
    grs80 = referenceEllipsoid('grs80','m');
    
    ax = axesm('globe','Geoid',grs80,'Grid','on', ...
        'GLineWidth',1,'GLineStyle','-',...
        'Gcolor','black','Galtitude',100);
    box off
    view(angle);
    %view([x_axis_body(2) x_axis_body(3)]);
    
    %set(gca, 'CameraPosition', [100 5000 2000]);
    
    xlim(x);
    ylim(y);
    xlabel(datestr(t))
    
    aniVec(i) = getframe(fig, [0 0 width height]);
    
    fprintf('\b#\n');
    
end

save('aniVec','aniVec', '-v7.3')
end
