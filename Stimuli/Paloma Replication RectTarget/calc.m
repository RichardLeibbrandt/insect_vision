formatSpec = '%d%d%d%d%f%f%d%d%d%d%d%d%d';
C = textscan(fopen('3000scan.txt'),formatSpec,'HeaderLines', 2, 'Delimiter',' ');
X = C{5};
Y = C{6};

xstart = X(54:90:end);
xend = X(90:90:end);

ystart = Y(54:90:end);
yend = Y(90:90:end);

dx = (xend-xstart);
dy = (yend-ystart);
a =dx.^2 + dy.^2;
b = sqrt(a);
v = mean(b)*10

