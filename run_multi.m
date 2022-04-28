close all
clear all
addpath toolbox

init_parameters

settings.t_start = datetime(HPS_convertMJD2CalendarDate(core_params.start_date(1))); 
settings.sim_step = core_params.dt_sim; 

settings.debug = 0;
settings.step = 50; 
settings.camMode = 'earthCentered';
settings.viewAngle = [0 0];
settings.rollAngle = 0;
settings.mainSat = 1;

load('states.mat');
s1.states = states;
s1.satModel = 'enc';
s1.satScale = 1.5e5;
s1.color = 'green';
s1.name = 'Grace';
s1.groundTrack = 1;

load('states-default.mat');
s2.states = states;
s2.satModel = 'enc';
s2.satScale = 1.5e5;
s2.color = '';
s2.name = '';
s2.groundTrack = 0;

satellites = [s1, s2];


frameVec = generateFrames(satellites, settings);
if isempty(frameVec)
	save('frameVec','frameVec', '-v7.3')
	saveVid('orbit-animation', 10, frameVec)
end
