close all
clear all
addpath toolbox

init_parameters

settings.t_start = datetime(HPS_convertMJD2CalendarDate(core_params.start_date(1))); 
settings.sim_step = core_params.dt_sim; 

settings.debug = 2;
settings.step = 250; 
settings.camMode = 'satFixed';
settings.viewAngle = [0 0];
settings.rollAngle = 0;
settings.mainSat = 2;

load('states.mat');
s1.states = states;
s1.satelliteModel = 'enc';
s1.satFactor = 2.5e5;
s1.baseColor = 'green';
s1.name = 'Grace';

load('states-default.mat');
s2.states = states;
s2.satelliteModel = 'enc';
s2.satFactor = 2.5e5;
s2.baseColor = '';
s2.name = '';

satellites = [s1, s2];


frameVec = generateFrames(satellites, settings);
if isempty(frameVec)
	save('frameVec','frameVec', '-v7.3')
	saveVid('orbit-animation', 10, frameVec)
end
