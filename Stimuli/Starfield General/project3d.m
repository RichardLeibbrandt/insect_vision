function [winX, winY] = project3d(objectPos,MP,viewport)
% Instead of glh (see "project3dOriginal" below), we are implementing the 
% algorithm from (https://gist.github.com/sixman9/871099).
% This is more standard (see also Mesa: 
% http://collagefactory.blogspot.com.au/2010/03/gluproject-source-code.html)
%
% It divides by the W to take to normalized device coordinates.
% The extra twist with this implementation (as opposed to standard Mesa) is 
% that it also flips Y around (changes sign) in the viewport transform. The default seems to be that Y is pointing 
% downward when we go to the screen.
% As with the Original, it is still "Much faster than calling gluProject in a loop"

%temp = model(:,1:3)*objectPos + model(:,4)*);
temp2 = MP*[objectPos;ones(1,size(objectPos,2))];
temp3 = temp2(4,:);  

r = (temp3~=0);

temp2 = temp2(:,r)./(ones(4,1)*temp3(r));
winX = (temp2(1,:)*.5+.5)*viewport(3)+viewport(1);
winY = (-temp2(2,:)*.5+.5)*viewport(4)+viewport(2);
%winZ = (1+temp2(3,:))*.5;

%{
function [winX winY] = project3dOriginal(objectPos,model,project,viewport)
% MATLAB optimzed version of glhProjectf
% (http://www.opengl.org/wiki/GluProject_and_gluUnProject_code)
% Much faster than calling gluProject in a loop
% (the glhProjectf (works only from perspective projection. With the
% orthogonal projection it gives different results than standard
% gluProject.)

temp = model(:,1:3)*objectPos + model(:,4)*ones(1,size(objectPos,2));

temp2 = project(1:3,:)*temp;
temp3 = -temp(3,:);

r = (temp3~=0);

temp2 = temp2(:,r)./(ones(3,1)*temp3(r));
winX = (temp2(1,:)*.5+.5)*viewport(3)+viewport(1);
winY = (temp2(2,:)*.5+.5)*viewport(4)+viewport(2);
%winZ = (1+temp2(3,:))*.5;
%}
