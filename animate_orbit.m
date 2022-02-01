% init
%run('~/work/tools/XHPS_2021/HPS_simulation/hps_startup.m')
close all;
load('states.mat');

%% local settings
vidName='jeb2';
duration = 10; % seconds
debug = 0;


%% settings
init_parameters
% sim_time %TODO: warning
settings.sim_step = core_params.dt_sim ;
settings.t_start = datenum(HPS_convertMJD2CalendarDate(core_params.start_date(1)));


settings.step = 50; % step includes only every step'th datapoint from states.mat to reduce filesize and rendering time
settings.resolution = '1080p';
settings.viewAngle = [0 20];
settings.camMode = 'fov';
settings.satelliteModel = 'stl';
settings.earthTransparency = 1;
settings.zoom = 10e2;
settings.progress = true;


frameVec = generateFrames(states, settings, debug);

if debug == 0
	% Save frameVec to file system
	save('frameVec','frameVec', '-v7.3')
	saveVid(vidName, duration, frameVec)
end
