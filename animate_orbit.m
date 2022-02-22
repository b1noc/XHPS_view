%% initialize workspace
close all
addpath toolbox
% Making sure that XHPS toolbox is loaded (path needs to be adjusted accordingly)
% This line can be omitted the hps_startupscript has been run manually or the path was added to Matlab by default
%run('<INSERT_PATH_HERE>/XHPS_2021/HPS_simulation/hps_startup.m')

% Reading in state matrix from Simulation
% This line can be ommitted when the states matrix already exists in the workspace 
% However to be able to create a video different settings at a later time it is advisable to save the states matrix in a file by executing 
% the following command: save('states.mat', 'states')
load('states.mat');

%% local settings
% Setting the output video name. CAREFUL!: This overwrites existing videofiles with that name on the same path
videoName='orbit-animation';

% Setting the video duration (shorter value results in quicker playback)
% This setting should be balancend with the settings.step value below to ensure smooth rendering
% It can be tweaked after the run when nessesary using the saveVid(<name>, <duration>, <frameVector>) function
duration = 10; % [s]

% Turning on the debug mode:
% - 0 = off -> generating full video (might take some time)
% - 1 = oneShot -> generates and displays the first frame (to check settings) and displays the loaded satellite Model if enabled
debug = 0;
gif = 0;

%% settings
% Run init_parameters to read in simulation parameters 
% This line can be omitted when the values already exist in the workspace or are being set manually
init_parameters

% Setting the starttime and the dt of every simulation step to calculate the time (and the earth's rotation if enabled)
settings.t_start = datetime(HPS_convertMJD2CalendarDate(core_params.start_date(1))); % format = [matlab datenum]
settings.sim_step = core_params.dt_sim ; % [s]
% Defining the stepwidth 
% Smaller steps takes longer to compile and cause larger filesizes 
settings.step = 50; % step includes only every step'th datapoint from states.mat to reduce filesize and rendering time

% Setting the video resolution
settings.resolution = '1080p';

% Loading a 3D model of the satellite
% - none -> no satellite model displayed
% - enc	 -> load Element-, Node-, and Reflexion coefficient tables (et.txt, nt.txt, ct.txt) provided by Ansys (path needs to be specified with settings.encFile)
% - stl  -> EXPERIMENTAL: load a local stl file (path needs to be specified with settings.stlfile) e
settings.satelliteModel = 'enc';

% Choosing the view perspective
% - fixed  -> set a fixed view in the ECS coordinate system
% - fov	   -> field of view: fly above the satellite with a view direction to the center of eart
% - follow -> follow the satellite 
% - pilot  -> set view satellite fixed and follow the viewfield of what a pilot would see
settings.camMode = 'fixed';

% Choosing the default view Angle of the scene. It only takes effect if the camMode = 'fixed'
settings.viewAngle = [0 20];

% Choosing the cameras proximity to the satellite
% Negative values result in zooming away (behind the satellite) positive values result in moving in front of the satellite
settings.zoom = -1e8;
settings.satFactor = 5e5;

frameVec = generateFrames(states, settings, debug);
if debug == 0
	% Save frameVec and video to file system
	save('frameVec','frameVec', '-v7.3')
	if ~gif 
		saveVid(videoName, duration, frameVec)
	else
		saveGif(videoName, frameVec)
	end
end

% To regenerate video with a different duration [s], simply run with the frameVec matrix in workspace
% load('frameVec'); % to load the frame Vector into the workspace
% saveVid('videoName', duration, frameVec) % to export the video

