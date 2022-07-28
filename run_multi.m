close all
clear all
addpath toolbox

run('xhps_results/init_parameters.m')

settings.t_start = datetime(HPS_convertMJD2CalendarDate(core_params.start_date(1))); 
settings.sim_step = core_params.dt_sim; 

settings.debug = 2;
settings.nodisplay = 0;
settings.progress = 1;
settings.duration = 10;
settings.fps = 60;
settings.resolution = '4k'
settings.camMode = 'earthCentered';
settings.viewAngle = [-30 40];
settings.rollAngle = 0;
settings.mainSat = 1;
settings.zoom = 1;
settings.earth = 1;
settings.ecef=1;
settings.eci=1;

load('xhps_results/states1.mat');
s1.states = states;
s1.satModel = 'stl';
s1.satScale = 5e4;
s1.color = 'green';
s1.name = 'Grace';
s1.groundTrack = 1;
s1.axisLength = 1;
s1.velVecLength = 1;

load('xhps_results/states2.mat');
s2.states = states;
s2.satModel = 'enc';
s2.satScale = 1e5;
s2.color = '';
s2.name = '';
s2.groundTrack = 0;
s2.axisLength = 1;
s2.velVecLength = 1;

satellites = [s1 s2];

frameVec = generateFrames(satellites, settings);
if ~isempty(frameVec)
	save('frameVec','frameVec', '-v7.3')
	saveVid('multi_4k_60fps', settings.duration, frameVec)
end
