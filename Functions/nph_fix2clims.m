
% Since matlab changed the behaviour of contourf, now regions above or
% below the colour scale are not drawn, rather than before where the were
% just saturated. I want them to be saturated so I'm going to limit the
% values to the limits.

function OUT = nph_fix2clims(IN,clims)

OUT = IN;

OUT(OUT < clims(1)) = clims(1);
OUT(OUT > clims(2)) = clims(2);









