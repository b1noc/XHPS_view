% init
%run('~/work/tools/XHPS_2021/HPS_simulation/hps_startup.m')
load('states.mat');

%% settings
vidName='vid1';

% sim_time %TODO: warning
sim_step = core_params.dt_sim ;
t_start = datenum(HPS_convertMJD2CalendarDate(core_params.start_date(1)));

% step includes only every step'th datapoint from states.mat to reduce filesize and rendering time
step = 10;
res = '720p';
view = [0 20];

frameVec = generateFrames(states, step, sim_step, t_start, view, res);

saveVid(vidName, 5, frameVec)



