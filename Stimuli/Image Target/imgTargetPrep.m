function [critInput] = rectTargetPrep(Parameters, ScreenData, StimSettings, NumSubframes)
%
% Prepares input parameters for rectTargetDraw

%--------------------------------------------------------------------------
% FlyFly v2
%
% Jonas Henriksson, 2010                                   info@flyfly.se
%--------------------------------------------------------------------------

if nargin<4
    NumSubframes = 1;
end


P.height      = Parameters(1,:);
P.width       = Parameters(2,:);
P.xpos        = Parameters(3,:);
P.ypos        = Parameters(4,:);
P.velocity    = Parameters(5,:)/NumSubframes;
P.angle       = Parameters(6,:);

numRuns = size(Parameters,2);

for k = 1:numRuns
    %imagePath = 'louie.jpg';
    imagePath        = StimSettings(k).path1{2};
    [I, ~, alpha] = imread(imagePath);
    if ~isempty(alpha)
        I(:,:,4) = alpha;
    end
    texturePtr(k) = Screen('MakeTexture', ScreenData.wPtr, I);
    aratio = size(I,1)/size(I,2);
    P.srcRect = [0; 0; size(I,1); size(I,2)]
    if P.width(k) == 0
        if P.height(k) == 0
            P.height(k) = size(I, 1);
            P.width(k) = size(I, 2);
        else
            P.width(k) = P.height(k) / aratio;
        end
    else
        if P.height(k) == 0
            P.height(k) = P.width(k) * aratio;
        else
            P.height(k) = min(P.height(k), P.width(k) * aratio);
            P.width(k) = P.height(k) / aratio;
        end
    end;
end

%how many fraction pixels target moves each frame
delta_x = cos(P.angle*pi/180);
delta_y = -sin(P.angle*pi/180);

critInput.RGB      = [0; 0; 0];
critInput.pos      = [P.xpos; P.ypos; ...
    P.xpos + P.width; P.ypos + P.height];

critInput.deltaPos = ScreenData.ifi*...
    [P.velocity.*delta_x; P.velocity.*delta_y; ...
    P.velocity.*delta_x; P.velocity.*delta_y];

critInput.srcRect = P.srcRect;
critInput.texturePtr      = texturePtr;
