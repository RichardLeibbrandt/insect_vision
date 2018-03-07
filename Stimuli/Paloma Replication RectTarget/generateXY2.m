function [Xstart, Ystart, Xend, Yend] = generateXY2(num_trials, speed, duration, width, height)
    % TODO: as it is this is slightly biased against near-vertical or
    % near-horizontal movements near the edges of the window.
    % not worth fixing as we just need it for this limited purpose.
    
    % used for random permutations
    hnt = round(num_trials/2);
    p = [zeros(1, hnt) ones(1, num_trials-hnt)];
    
    % distance to move per trial
    dist = floor(speed*duration);
    
    a=linspace(0,2*pi,360);
    Xdiff = round(dist*cos(a));
    Ydiff = round(dist*sin(a));
    diff = [Xdiff' Ydiff'];
    
    % choose random position 1 within maxdist of the window edge (so we don't need
    % to check bounds!)
    X1 = dist + randi(width-2*dist, 1, num_trials);
    Y1 = dist + randi(height-2*dist, 1, num_trials);
    
    % choose position 2's X to be within maxdist from position 1's X
    X2 = X1 - dist + randi(2*dist, 1, num_trials);
    % calculate position 2's Y from that (randomly choose to go above or below position 1)
    flipY = 2*p(randperm(num_trials))-1;
    Y2 = round(Y1 + flipY.*sqrt(dist^2-(X2-X1).^2));
    
    % randomly choose whether to go from pos 1 to pos 2, or vice versa
    % (otherwise we'd be moving outwards more often than inwards)
    flip12 = find(p(randperm(num_trials))==1);
    Xstart = X1; Ystart = Y1; Xend = X2; Yend = Y2;
    Xstart(flip12) = X2(flip12); Ystart(flip12) = Y2(flip12); Xend(flip12) = X1(flip12); Yend(flip12) = Y1(flip12);
end