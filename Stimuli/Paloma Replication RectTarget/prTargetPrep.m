function [critInput] = prTargetPrep(Parameters, ScreenData, StimSettings, NumSubframes)
%
% Prepares input parameters for rectTargetDraw

%--------------------------------------------------------------------------
% FlyFly v3.1
%
% Richard Leibbrandt, 2017
%------------------------------------------------------------`--------------
NUM_TRIALS = 3000;
SPEED = 900;
MOVE_DURATION = 0.1;
SUBWIDTH = 200;
SUBHEIGHT = 300;

target_halfwidth = 7;
target_halfheight = 7;

if ScreenData.usePartial
    width   = ScreenData.partial(3);
    height  = ScreenData.partial(4);
else
    screenRes = get(0, 'ScreenSize');   %[1 1 screenWidth screenHeight]
    width   = screenRes(3);
    height  = screenRes(4);  
end
width_adj = floor(width/2 - SUBWIDTH/2);
height_adj = max(height-SUBHEIGHT-120, floor(height/2 - SUBHEIGHT/2));

[Xstart, Ystart, Xend, Yend] = ...
    generateXY(NUM_TRIALS, SPEED, MOVE_DURATION, SUBWIDTH, SUBHEIGHT);
 
fps = 1/ScreenData.ifi;
STATIONARY_FRAMES = round(0.15*fps);
MOVING_FRAMES = round(0.1*fps);
FRAMES_PER_TRIAL = STATIONARY_FRAMES + MOVING_FRAMES;
NUM_FRAMES = FRAMES_PER_TRIAL * NUM_TRIALS;

D = zeros(NUM_FRAMES, 2);

for t = 1:NUM_TRIALS
    disp(t);
    iblock = zeros(FRAMES_PER_TRIAL, 2);
    srange = 1:STATIONARY_FRAMES;
    % need the zero in the following line as we will double-write the first
    % one! hack, hack
    mrange = STATIONARY_FRAMES + (0:MOVING_FRAMES);
    iblock(srange, 1) = Xstart(t);
    iblock(srange, 2) = Ystart(t); 
    iblock(mrange, 1) = round(linspace(Xstart(t), Xend(t), MOVING_FRAMES+1));
    iblock(mrange, 2) = round(linspace(Ystart(t), Yend(t), MOVING_FRAMES+1));
    
    D(((t-1)*FRAMES_PER_TRIAL+1):t*FRAMES_PER_TRIAL, :) = iblock;
	iblock = [iblock(:,1)-target_halfwidth+width_adj,  iblock(:,2)-target_halfheight+height_adj, ...
             iblock(:,1)+target_halfwidth+width_adj,   iblock(:,2)+target_halfheight+height_adj];

    critInput(t).pos = iblock;
    critInput(t).rgb = [0; 0; 0];
end

outfile = 'debug/3000scan.flinders.txt';
D = [zeros(NUM_FRAMES, 4) D zeros(NUM_FRAMES, 7)];
dlmwrite(outfile, D, ' ');

end


