function [critInput] = prTargetPrepOld(Parameters, ScreenData, StimSettings, NumSubframes)
%
% Prepares input parameters for rectTargetDraw

%--------------------------------------------------------------------------
% FlyFly v2
%
% Jonas Henriksson, 2010                                   info@flyfly.se
%------------------------------------------------------------`--------------
NUM_TRIALS = 3000;
SPEED = 900;
MOVE_DURATION = 0.1;
SUBWIDTH = 1280;
SUBHEIGHT = 720;


formatSpec = '%d%d%d%d%f%f%d%d%d%d%d%d%d';
C = textscan(fopen('3000scan.txt'),formatSpec,'HeaderLines', 2, 'Delimiter',' ');
X = C{5};
Y = C{6};


target_halfwidth = 7;
target_halfheight = 7;

fps = 1/ScreenData.ifi;
num_frames = round(0.25*fps);
ls = linspace(1, 90, num_frames);

if ScreenData.usePartial
    width   = ScreenData.partial(3);
    height  = ScreenData.partial(4);
else
    screenRes = get(0, 'ScreenSize');   %[1 1 screenWidth screenHeight]
    width   = screenRes(3);
    height  = screenRes(4);  
end
width_adj = floor(width/2 - SUBWIDTH/2);
height_adj = floor(height/2 - SUBHEIGHT/2);

for k = 1:3000
    v = [X(90*(k-1)+1:90*k) Y(90*(k-1)+1:90*k)];
    iblock = interp1(1:90, v, ls, 'linear');
    iblock = [iblock(:,1)-target_halfwidth+width_adj,  iblock(:,2)-target_halfheight+height_adj, ...
             iblock(:,1)+target_halfwidth+width_adj,   iblock(:,2)+target_halfheight+height_adj];

    critInput(k).pos = iblock;
    critInput(k).rgb = [0; 0; 0];
end

end


