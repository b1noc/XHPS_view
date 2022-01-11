% init
%run('~/work/tools/XHPS_2021/HPS_simulation/hps_startup.m')
load('states.mat');

%% settings
vidName='vid2';

% sim_time %TODO: warning
sim_step = core_params.dt_sim ;
t_start = datenum(HPS_convertMJD2CalendarDate(core_params.start_date(1)));

% step includes only every step'th datapoint from states.mat to reduce filesize and rendering time
step = 5;
res = '1080p';
view = [0 20];
duration = 10; % seconds
debug = 1;

frameVec = generateFrames(states, step, sim_step, t_start, view, res, debug);

if debug ~= 1
	saveVid(vidName, duration, frameVec)
end
