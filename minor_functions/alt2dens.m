


% CALCULATE APPROXIMATE ATOMSPHERIC DENSITY BY ASSUMING HYDROSTATIC
% EQUILIBRIUM AND HENCE USING THE HYDROSTATIC EQUATION.

% We also use the scale heights used by alt2pres_complex.

% dens(z) = (p0)/(g*H0)*exp(-z/H0);

% OUTPUT IN kg/m3.

function dens = alt2dens(alt)

dens = zeros(size(alt));


% These are the levels from the wiki page for the BAROMETRIC FORMULA

g = 9.8; %ms-2 Could reduce g with height to get a slightly more accurate answer...

% define pressure "base" levels:

p0 = 1013.25; % (hPa) (sea level)
p1 = 226.3210; 
p2 = 54.7489;
p3 = 8.6802;
p4 = 1.1091;
p5 = 0.6694;
p6 = 0.0396;

% and their corresponding altitude levels (km)

z0 = 0;
z1 = 11;
z2 = 20;
z3 = 32;
z4 = 47;
z5 = 51;
z6 = 71;

% now calculate scale heights...

H0 = - (z1-z0) / (ln(p1/p0));
H1 = - (z2-z1) / (ln(p2/p1));
H2 = - (z3-z2) / (ln(p3/p2));
H3 = - (z4-z3) / (ln(p4/p3));
H4 = - (z5-z4) / (ln(p5/p4));
H5 = - (z6-z5) / (ln(p6/p5));

H6 = H5;


% and now calculate the altitudes from pres2alt.m...



for i = 1:length(alt),

    z = alt(i);
    
if z <= z0, % Quick clause to cope with underground things. Just set the density to be surface density and move on.
    
    dens(i) = 1.4; % kg/m3
    
    continue
    
end
    
    
    
if z >= z0 && z < z1, % LAYER 0
    
    ref_pres = p0; ref_alt = z0; 
    
    H = H0;

end

if z >= z1 && z < z2, % LAYER 1
    
    ref_pres = p1; ref_alt = z1; 
    
    H = H1;

end

if z >= z2 && z < z3, % LAYER 2
    
    ref_pres = p2; ref_alt = z2; 
    
    H = H2;

end

if z >= z3 && z < z4, % LAYER 3

    ref_pres = p3; ref_alt = z3; 
    
    H = H3;

end

if z >= z4 && z < z5, % LAYER 4
    
    ref_pres = p4; ref_alt = z4; 
    
    H = H4;

end

if z >= z5 && z < z6, % LAYER 5
    
    ref_pres = p5; ref_alt = z5; 
    
    H = H5;

end

if z >= z6, % LAYER 6 AND ABOVE
    
    ref_pres = p6; ref_alt = z6; 
    
    H = H6;

end

% I think I should be using ref_pres, but when you plot
% >> figure; plot(alt2dens(0:0.1:50),0:0.1:50)
% there's a couple of nasty kinks in the otherwise smooth decrease.
% This is fixed by always using p0, so let's use that.

dens(i) = ((p0*100)/(g*(H*1000))) * exp(-z/H);


end % end cycling through all of alt's element(s)...



