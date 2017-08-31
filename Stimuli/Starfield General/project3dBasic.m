function [winX, winY] = project3dBasic(objectPos,model,project,viewport)
% MATLAB optimzed version of glhProjectf
% (http://www.opengl.org/wiki/GluProject_and_gluUnProject_code)
% Much faster than calling gluProject in a loop
% (the glhProjectf (works only from perspective projection. With the
% orthogonal projection it gives different results than standard
% gluProject.)
winX = zeros(1,size(objectPos,2));
winY = zeros(1,size(objectPos,2));

for r = 1:size(objectPos,2)
    x = objectPos(1,r);
    y = objectPos(2,r);
    z = objectPos(3,r);
    [ winX(r), winY(r), ~, ~ ] = gluProject( x, y, z, model, project, viewport );
end