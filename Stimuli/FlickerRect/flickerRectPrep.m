function [critInput] = flickerRectPrep(Parameters, ScreenData, StimSettings, NumSubframes)
%
% Prepares input parameters for flickerRectDraw

%--------------------------------------------------------------------------
% FlyFly v2
%
% Jonas Henriksson, 2010                            info@flyfly.se
%--------------------------------------------------------------------------

if nargin<4
    NumSubframes = 1;
end

[R C] = size(Parameters);

P.height      = Parameters(1,:);
P.width       = Parameters(2,:);
P.xpos        = Parameters(3,:);
P.ypos        = Parameters(4,:);
P.nFrame      = Parameters(5,:)/NumSubframes;
P.bright1     = Parameters(6,:);
P.bright2     = Parameters(7,:);

useX = zeros(1, C);
useY = zeros(1, C);

for k = 1:C
    if StimSettings(k).box1{2} %use x function
        useX(k)    = 1;
        xOffset{k} = StimSettings(k).edit1{2};
    else
        useX(k) = 0;
        xOffset = '1';
    end

    if StimSettings(k).box2{2} %use y function
        useY(k)    = 1;
        yOffset{k} = StimSettings(k).edit2{2};
    else
        useY(k) = 0;
        yOffset = '1';
    end
end

critInput.useX = useX;
critInput.useY = useY;

critInput.xOffset = xOffset;
critInput.yOffset = yOffset;

critInput.RGB1      = [P.bright1; P.bright1; P.bright1];
critInput.RGB2      = [P.bright2; P.bright2; P.bright2];

critInput.pos      = [P.xpos - P.width/2; P.ypos - P.height/2; ...
    P.xpos + P.width/2; P.ypos + P.height/2];

critInput.nFrame = P.nFrame;