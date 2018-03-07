function [x, y, z, sizes] = starSeedDebug(xmin, xmax, xspace, ymin, ymax, yspace, zmin, zmax, zspace, size)

xrange = xmin:xspace:xmax;
yrange = ymin:yspace:ymax;
zrange = zmin:zspace:zmax;

x = repmat(xrange, 1, numel(yrange));
y = repmat(yrange, numel(xrange), 1);
y = y(:)';

x = repmat(x, 1, numel(zrange));
y = repmat(y, 1, numel(zrange));
z = repmat(zrange, numel(xrange)*numel(yrange), 1);
z = z(:)';

numDots = numel(x);
sizes = ones(1, numDots) * size;

end



