


% output min and max of a vector or matrix;

function OUT = minmax(IN)

minn = nanmin(IN(:));
maxx = nanmax(IN(:));

OUT = [minn maxx];


