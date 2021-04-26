


% outputs plus minus a scalar so you don't have to write it all the time :)

function OUT = pm(IN)

if nargin == 0
    IN = 1;
end

IN = abs(IN);

OUT = [-IN IN];

end

