close all
clear s
addpath toolbox

load('states.mat');
init_parameters

s.t_start = datetime(HPS_convertMJD2CalendarDate(core_params.start_date(1))); 
s.sim_step = core_params.dt_sim; 

s.debug = 0;
s.step = 50; 
s.satelliteModel = 'enc';
s.satFactor = 5e5;
s.camMode = 'earthCentered';
s.viewAngle = [0 0];
s.rollAngle = 0;

frameVec = generateFrames(states, s);
if exist('frameVec','var') == 1
	%save('frameVec','frameVec', '-v7.3')
	saveVid('orbit-animation', 10, frameVec)
	%saveGif(videoName, frameVec)
end
