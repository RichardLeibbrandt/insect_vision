% Starfield 2 is spherical.
% Fly does not move. Background moves in ONE of the 6 available motions.
% Combined background movements may work but have not been tested yet.

function output = starfieldPrep(Parameters, ScreenData, StimSettings, NumSubframes)
%
% Takes input from the user inteface and makes pre computations.
%

%--------------------------------------------------------------------------
% FlyFly v3.1
%
% Jonas Henriksson, 2010                                   info@flyfly.se
%--------------------------------------------------------------------------

if nargin<4
    NumSubframes = 1;
end

global GL;

DENSITY_CONVERSION = 1e-4;

starfieldUserSettings; %loads user settings from file. loaded files in CAPS

[~, numRuns] = size(Parameters);

P.dotSize     = Parameters(1,1) * ones(1,length(Parameters(1,:)));    % Note we take the first trial's dot size for everything
P.density     = DENSITY_CONVERSION*Parameters(2,1) * ones(1,length(Parameters(2,:)));    % Take the first trial's density for everything

P.rl          = Parameters(3,:)/NumSubframes;   % sideslip
P.ud          = Parameters(4,:)/NumSubframes;   % lift
P.fb          = Parameters(5,:)/NumSubframes;   % thrust
P.pitch       = Parameters(6,:)/NumSubframes;
P.yaw         = Parameters(7,:)/NumSubframes;
P.roll        = Parameters(8,:)/NumSubframes;
P.background_noise = Parameters(9,:)*NumSubframes;
P.retain      = Parameters(10,:);
P.t           = Parameters(11,:)*NumSubframes;

% z-clipping planes
zFar = 200;
%zNear = 3;
zNear = ScreenData.flyDistance;

% Field of view (y)
fovy = 2*atand(0.5*ScreenData.monitorHeight/ScreenData.flyDistance);
% Aspect ratio
ar = ScreenData.rect(3) / ScreenData.rect(4);



P.monHeight    = ScreenData.monitorHeight;

P.monWidth = P.monHeight * ar;
P.viewDistance = zFar;

box_x = P.monWidth*P.viewDistance/ScreenData.flyDistance;
box_y = P.monHeight*P.viewDistance/ScreenData.flyDistance;
box_z = 2*P.viewDistance;

%P.dotSize      = min(max(1, P.dotSize),63);    %min size of dots in drawDots = 1, max = 63

ifi = ScreenData.ifi;

pxPerCm = ScreenData.rect(4) ./ P.monHeight;

%[center(1), center(2)] = RectCenter(ScreenData.rect);
center = ScreenData.flyPos(1:2);

hPx = ScreenData.flyPos(4)-center(2); % px offset from center
h   = hPx/pxPerCm;

cameraLookAt = [(ScreenData.flyPos(3:4)-ScreenData.flyPos(1:2))./pxPerCm ScreenData.flyDistance];
tiltAngleX = atand(((ScreenData.flyPos(3)-ScreenData.flyPos(1))./pxPerCm)/ScreenData.flyDistance);
tiltAngleY = atand(((ScreenData.flyPos(4)-ScreenData.flyPos(2))./pxPerCm)/ScreenData.flyDistance);

center(2) = center(2) + h*pxPerCm;

%gamma1 = (atan(h/ScreenData.flyDistance))/180*pi;
%gamma2 = (pi/2 - gamma1);

box_x = P.monWidth*P.viewDistance/ScreenData.flyDistance;
box_y = P.monHeight*P.viewDistance/ScreenData.flyDistance;
box_z = 2*P.viewDistance;

P.viewDistX = box_x/2;
P.viewDistY = box_y/2;
P.viewDistZ = box_z/2;

seconds = P.t*ifi;  % estimated no of seconds to perform the rotation
full_roll = deg2rad(P.roll.*seconds);
full_yaw = deg2rad(P.yaw.*seconds);
full_pitch = deg2rad(P.pitch.*seconds);
full_rl = -P.rl.*seconds;
full_ud = P.ud.*seconds;
full_fb = -P.fb.*seconds;

% THE FOLLOWING DOESN'T WORK YET, SO FALLING BACK ON THE OLD WAY OF
% AVOIDING MOVING OUT OF VISUAL SPACE - THIS IS STILL TO DO!
%{
% CALCULATE THE SIZE OF THE STARFIELD SPACE
% We need to calculate the size of the required underlying visual space, so that we don't
% "go over the edge" and run out of stars, potentially causing inconsistencies 
% in visual density.

% boxes stores the maximum extent of the visible space in each of the X Y Z dimensions, 
% given how far the origin has moved over the course of several trials.
% (the space is retained from one trial to the next when retain is true)
boxes = cell(1,numRuns);
% we only store something in boxes whenever the star seed is reinitialized (i.e. k==1, or retain is false)
% firstIndex is the index where this occurs each time
firstIndex = 1; 

% we keep track of where the origin moves to, as well as the furthest
% extent that the origin has moved so far in any of the 3 axis directions.
% the value of boxes will be each of these X Y Z displacements, with the
% distance to the far plane (zFar) added to each, multiplied by 2 to
% produce a symmetrical box for use by starSeed2.
ocoords = [0 0 0];  % coordinates of the origin
mins = [0 0 0];     % furthest extent of the origin into the negative x,y,z axes
maxes = [0 0 0];    % furthest extent of the origin into the positive x,y,z axes

% calculate boxes
for k =1:numRuns
    if P.yaw(k) ~= 0
        [ocoords(1), ocoords(3)] = rotateXY(ocoords(1), ocoords(3), full_yaw(k));
    elseif P.pitch(k) ~= 0
        [ocoords(2), ocoords(3)] = rotateXY(ocoords(2), ocoords(3), full_pitch(k));
    elseif P.roll(k) ~= 0
        [ocoords(1), ocoords(2)] = rotateXY(ocoords(1), ocoords(2), full_roll(k));
    else
        ocoords = ocoords + [full_rl(k) full_ud(k) full_fb(k)];
    end
    mins = min([mins; ocoords]);
    maxes = max([maxes; ocoords]);
    if k==numRuns || ~P.retain(k)
        boxes{firstIndex} = 2*(max([abs(mins); abs(maxes)]) + zFar);
        ocoords = [0 0 0];
        mins = [0 0 0];
        maxes = [0 0 0];
        firstIndex = k + 1;
    end
end;
%}

Screen('BeginOpenGL', ScreenData.wPtr);

glViewport(0, 0, ScreenData.rect(3), ScreenData.rect(4));

glDisable(GL.LIGHTING);

glEnable(GL.BLEND);
glBlendFunc(GL.SRC_ALPHA, GL.ONE_MINUS_SRC_ALPHA);

glMatrixMode(GL.PROJECTION);
glLoadIdentity;

gluPerspective(fovy, ar, zNear, zFar);
%glOrtho(-P.monWidth/2, P.monWidth/2, -P.monHeight/2, P.monHeight/2, zNear, zFar);
% glRotatef(-tiltAngleX,1,0,0);
% glRotatef(tiltAngleY,0,1,0);

glMatrixMode(GL.MODELVIEW);
glLoadIdentity;

gluLookAt(0,0,0, cameraLookAt(1),cameraLookAt(2),cameraLookAt(3), 0,1,0);
%gluLookAt(0,0,0,  0,0,1,  0,1,0);

% These are used in gluProject (3d -> 2d coordinates)
% Get the projection matrix
projection = glGetDoublev(GL.PROJECTION_MATRIX);
projectionMatrix = reshape(projection,4,4);
% Get the modelview matrix
modelview = glGetDoublev(GL.MODELVIEW_MATRIX);
modelviewMatrix = reshape(modelview,4,4);
% Get the viewport
viewport = glGetIntegerv(GL.VIEWPORT);
viewport = double(viewport);

MP = projectionMatrix*modelviewMatrix;

output(numRuns) = struct('xymatrix',[],'color',[],'dotsize',[],'center',[]);

origin_pos = [0 0 0];
crash_alert = false;

for k=1:numRuns     % for each trial
    if (k==1) || ~P.retain(k-1)
        % THE FOLLOWING 2 LINES RELATE TO THE APPROACH OF CALCULATING THE
        % STARFIELD VISUAL SPACE SIZE
        %b = boxes{k};
        %[x, y, z, sizes] = starSeed2(P.density(1), P.dotSize(1), b(1), b(2), b(3), zFar);
        [x, y, z, sizes] = starSeed2(P.density(1), P.dotSize(1), box_x, box_y, box_z, zFar);
    end;
    
    N = P.t(k);
    
    output(k).color = cell(1,N);
    output(k).dotsize = cell(1,N);
    output(k).xymatrix = cell(1,N);
    
    cumulative_roll = 0;
    cumulative_yaw = 0;
    cumulative_pitch = 0;
    cumulative_rl = 0;
    cumulative_ud = 0;
    cumulative_fb = 0;
    
    sx = cell(1,N);
    sy = cell(1,N);
    sz = cell(1,N);
    %totaldots = zeros(1,N);
    for n=1:N   % for each frame
        
        noise = 2*P.background_noise(k)*rand() - P.background_noise(k);
        
        %ROTATE
        %yaw - affects x and z
        if P.yaw(k) ~= 0
            inc_yaw = (full_yaw(k) - cumulative_yaw)/(N-n+1);
            inc_yaw = inc_yaw + noise*inc_yaw;
            [x z] = rotateXY(x, z, inc_yaw);
            cumulative_yaw = cumulative_yaw + inc_yaw;
        end
        
        %pitch - affects y and z
        if P.pitch(k) ~= 0
            inc_pitch = (full_pitch(k) - cumulative_pitch)/(N-n+1);
            inc_pitch = inc_pitch + noise*inc_pitch;
            [y z] = rotateXY(y, z, inc_pitch);
            cumulative_pitch = cumulative_pitch + inc_pitch;
        end
        
        %roll - affects x and y
        if P.roll(k) ~= 0
            inc_roll = (full_roll(k) - cumulative_roll)/(N-n+1);
            inc_roll = inc_roll + noise*inc_roll;
            [x y] = rotateXY(x, y, inc_roll);
            cumulative_roll = cumulative_roll + inc_roll;
        end
        
        %TRANSLATE
        if P.rl(k) ~= 0
            inc_rl = (full_rl(k) - cumulative_rl)/(N-n+1);
            inc_rl = inc_rl*(1+noise);
            cumulative_rl = cumulative_rl + inc_rl;
            x = x + inc_rl;
        end;
        
        if P.ud(k) ~= 0
            inc_ud = (full_ud(k) - cumulative_ud)/(N-n+1);
            inc_ud = inc_ud*(1+noise);
            cumulative_ud = cumulative_ud + inc_ud;
            y = y + inc_ud;
        end;
        
        if P.fb(k) ~= 0
            inc_fb = (full_fb(k) - cumulative_fb)/(N-n+1);
            inc_fb = inc_fb*(1+noise);
            cumulative_fb = cumulative_fb + inc_fb;
            z = z + inc_fb;
        end;
        
        % Move stray objects back - AVOIDS MOVING OUT OF
        % THE VISUAL SPACE
        
        x(abs(x)> P.viewDistX) = -(x(abs(x)> P.viewDistX));
        y(abs(y)> P.viewDistY) = -(y(abs(y)> P.viewDistY));
        z(abs(z)> P.viewDistZ) = -(z(abs(z)> P.viewDistZ));
        
              
        % CULLING
        % Frustum culling, [fx; fy; fz] are all the points that are inside the frustum
        % http://www.lighthouse3d.com/tutorials/view-frustum-culling/radar-approach-testing-points-ii/
        h = z*tand(fovy/2);
        w = h*ar;
        
        indices = and(z>zNear,z<zFar);
        indices = and(indices,and(y>-h,y<h));
        indices = and(indices,and(x>-w,x<w));
        
        fx = x(indices);
        fy = y(indices);
        fz = z(indices);
        fs = sizes(indices);
      
        %fx=x;fy=y;fz=z;fs=sizes;
        % CALCULATE ALL DISTANCES FROM THE ORIGIN
        distances = sqrt(fx.^2+fy.^2+fz.^2);                
        
        % CLIPPING
        % Don't draw particles that are too far away
        indices = distances<zFar;
        fx = fx(indices);
        fy = fy(indices);
        fz = fz(indices);
        fs = fs(indices);
        distances = distances(indices);       
        
        % SORT FROM FAR TO NEAR
        % So that the nearer dots are displayed in front of the further
        % ones (drawn last)
        [~, sort_idx] = sort(distances, 'descend');
        fx = fx(sort_idx);
        fy = fy(sort_idx);
        fz = fz(sort_idx);
        fs = fs(sort_idx);
        distances = distances(sort_idx);
        
        % THE FOLLOWING IS FOR DEBUGGING PURPOSES
        if isempty(fx) && ~crash_alert
            fprintf('CRASH COMING in TRIAL %d, FRAME %d', k, n);
            crash_alert = true;
            crash_site = origin_pos + [cumulative_rl cumulative_ud cumulative_fb]
        end;
        
        % PROJECTION
        [sx{n}, sy{n}] = project3d([fx; fy; fz], MP, viewport);
        %[sx{n}, sy{n}] = project3dOriginal([fx; fy; fz],modelviewMatrix,projectionMatrix,viewport);
        %[sx{n}, sy{n}] = project3dBasic([fx; fy; fz],modelviewMatrix,projectionMatrix,viewport);

        output(k).xymatrix{n} = [sx{n}; sy{n}];
        %totaldots(n) = size(sx{n},2);
        % Color and size of dots
        output(k).dotsize{n} = max(min(63, fs./distances),1);
        output(k).color{n} = repmat(255*distances/zFar, 3, 1);
    end;
    origin_pos = origin_pos + [cumulative_rl cumulative_ud cumulative_fb]
    
    output(k).center = [];
end

%fprintf('AVERAGE NUMBER OF DOTS (totaldots) was %f\n', mean(totaldots));

Screen('EndOpenGL', ScreenData.wPtr);



function [x1,y] = rotateXY(x, y, angle)
% rotates coordinate (x, y) angle degrees around origin

x1 = x*cos(angle) - y*sin(angle);
y  = y*cos(angle) + x*sin(angle);



% NOTE I had the following block at the beginning of the target stuff:
%{
        % THIS BLOCK NEEDS TO BE COMMENTED IN/OUT
        % TO TOGGLE WHETHER THE TARGET MOVES WITH THE BACKGROUND OR NOT
        % sort of - except that the meaning of yaw, roll etc has changed -
        % need to modify them by multiplying by ifi as before
        if P.yaw(k) ~= 0
            [target_start(1), target_start(3)] = rotateXY(target_start(1), target_start(3), P.yaw(k));
            [target_end(1), target_end(3)] = rotateXY(target_end(1), target_end(3), P.yaw(k));
            [target_pos(1), target_pos(3)] = rotateXY(target_pos(1), target_pos(3), P.yaw(k));
        end
      
       %pitch - affects y and z
        if P.pitch(k) ~= 0
            [target_start(2), target_start(3)] = rotateXY(target_start(2), target_start(3), P.pitch(k));
            [target_end(2), target_end(3)] = rotateXY(target_end(2), target_end(3), P.pitch(k));
            [target_pos(2), target_pos(3)] = rotateXY(target_pos(2), target_pos(3), P.pitch(k));
        end

         %roll - affects x and y
        if P.roll(k) ~= 0
            [target_start(1), target_start(2)] = rotateXY(target_start(1), target_start(2), P.roll(k));
            [target_end(1), target_end(2)] = rotateXY(target_end(1), target_end(2), P.roll(k));
            [target_pos(1), target_pos(2)] = rotateXY(target_pos(1), target_pos(2), P.roll(k));
        end
        
        %disp('target translate'), tic
        target_start = target_start + [P.rl(k); P.ud(k); P.fb(k)];  % fly translation
        target_end = target_end + [P.rl(k); P.ud(k); P.fb(k)];  % fly translation
        target_pos = target_pos + [P.rl(k); P.ud(k); P.fb(k)];  % fly translation
%}