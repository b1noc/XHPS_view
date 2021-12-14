%% Plot orbit with satellite body frame

% seconds between body frame plots
step = 200; 
% number of plotted body frames
n_plots = 36;

% number of plotted body frames needs to match with step and simulation time:
if (step*n_plots/core_params.dt_sim) > length(states)
     n_plots = floor(length(states)*core_params.dt_sim/step);
end

state = states;

%%
figure
%plot orbit
plot3(state(:,11),state(:,12),state(:,13))
hold on 
%plot starting (initial) Point
plot3(state(1,11),state(1,12),state(1,13),'color', 'red', 'Marker', 'o', 'MarkerSize', 7 , 'linewidth',2)
axis equal
hold on
%plot initial velocity vector
quiver3(state(1,11),state(1,12),state(1,13),state(1,8)*500,state(1,9)*500,state(1,10)*500,'linewidth',2)
hold on 

%plot inertial kos
a=6771000;

axis_length = 1*a/8; %length of coordinate system axis
x_axis= [axis_length;0;0];
y_axis= [0;axis_length;0];
z_axis= [0;0;axis_length];

quiver3(0,0,0,axis_length*1.5,0,0,'linewidth',3,'color', 'blue')
hold on 
quiver3(0,0,0,0,axis_length*1.5,0,'linewidth',3,'color', 'red')
hold on 
quiver3(0,0,0,0,0,axis_length*1.5,'linewidth',3,'color', 'green')

%plot body kos at time step i_time
i_time = 1;

q_i = state(i_time,4:7);

%transform inertial kos to body fixed by quaternion multiplication
    x_axis_body = HPS_transformVecByQuatTransposed(x_axis,q_i);
    y_axis_body = HPS_transformVecByQuatTransposed(y_axis,q_i);
    z_axis_body = HPS_transformVecByQuatTransposed(z_axis,q_i);

    
% Vektor [1*a/8,0,0] turn with quaternion 

quiver3(state(i_time,11),state(i_time,12),state(i_time,13),x_axis_body(1),x_axis_body(2),x_axis_body(3),'linewidth',2,'color', 'blue')
hold on 
quiver3(state(i_time,11),state(i_time,12),state(i_time,13),y_axis_body(1),y_axis_body(2),y_axis_body(3),'linewidth',2,'color', 'red')
hold on 
quiver3(state(i_time,11),state(i_time,12),state(i_time,13),z_axis_body(1),z_axis_body(2),z_axis_body(3),'linewidth',2,'color', 'green')
hold on


%% Plot body fixed frame for more time steps

for i = 1: n_plots 
    
    i_time = ceil(i * step / core_params.dt_sim);

    q_i = state(i_time,4:7);

%     %transform inertial kos to body fixed by quaternion multiplication
%     %(keine darstellung des Orbital Kos in eci daher vecbyquattransposed)
   x_axis_body = HPS_transformVecByQuatTransposed(x_axis,q_i);
   y_axis_body = HPS_transformVecByQuatTransposed(y_axis,q_i);
   z_axis_body = HPS_transformVecByQuatTransposed(z_axis,q_i);
    

    % Vektor [1*a/8*1000,0,0] turn with quaternion 

    quiver3(state(i_time,11),state(i_time,12),state(i_time,13),x_axis_body(1),x_axis_body(2),x_axis_body(3),'linewidth',2,'color', 'blue')
    hold on 
    quiver3(state(i_time,11),state(i_time,12),state(i_time,13),y_axis_body(1),y_axis_body(2),y_axis_body(3),'linewidth',2,'color', 'red')
    hold on 
    quiver3(state(i_time,11),state(i_time,12),state(i_time,13),z_axis_body(1),z_axis_body(2),z_axis_body(3),'linewidth',2,'color', 'green')
    hold on 
    %plot des normalenvektors
   % quiver3(state(i_time,11),state(i_time,12),state(i_time,13),n_i(i_time,1)*7e5,n_i(i_time,2)*7e5,n_i(i_time,3)*7e5,'linewidth',2,'color', 'black')
   % hold on 


end


%plot earth 
grs80 = referenceEllipsoid('grs80','m');

ax = axesm('globe','Geoid',grs80,'Grid','on', ...
    'GLineWidth',1,'GLineStyle','-',...
    'Gcolor','black','Galtitude',100);
sax.Position = [0 0 1 1];
axis equal
box off
view(3)

clear step n_plots a x_axis_body y_axis_body z_axis_body q_i i_time i state z_axis y_axis x_axis axis_length


    
    
