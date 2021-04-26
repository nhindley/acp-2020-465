

% 20181121 - NEW VERSION USING INTERPOLANT

function H = scale_height(alt)


% Pressure and hieght levels from US Standard Atmopshere:
p = [1013.25 226.3210 54.7489 8.6802 1.1091 0.6694 0.0396];
z = [0 11 20 32 47 51 71];

dz = diff(z);
zmid = z(1:6) + dz/2;
zmid(end) = z(end);
zmid = [0 zmid];

% compute scale factor for each height bin:
h = nan(1,6);
for i = 1:6
    h(i) = - dz(i) / ln(p(i+1)/p(i));
end
h = [h(1) h];

% interpolate input values on to this:
F = griddedInterpolant(zmid,h,'pchip','none');
H = F(alt);

% Set upper limit:
H(alt > z(7)) = h(7);



% 
% 
% p0 = 1013.25; % (hPa) (sea level)
% p1 = 226.3210; 
% p2 = 54.7489;
% p3 = 8.6802;
% p4 = 1.1091;
% p5 = 0.6694;
% p6 = 0.0396;
% 
% % and their corresponding altitude levels (km)
% 
% z0 = 0;
% z1 = 11;
% z2 = 20;
% z3 = 32;
% z4 = 47;
% z5 = 51;
% z6 = 71;
% 
% % now calculate scale heights...
% 
% H0 = - (z1-z0) / (ln(p1/p0));
% H1 = - (z2-z1) / (ln(p2/p1));
% H2 = - (z3-z2) / (ln(p3/p2));
% H3 = - (z4-z3) / (ln(p4/p3));
% H4 = - (z5-z4) / (ln(p5/p4));
% H5 = - (z6-z5) / (ln(p6/p5));
% 
% H6 = H5;
% 
% 
% if alt <= z1, H = H0; end
% if alt > z1 && alt <= z2, H = H1; end
% if alt > z2 && alt <= z3, H = H2; end
% if alt > z3 && alt <= z4, H = H3; end
% if alt > z4 && alt <= z5, H = H4; end
% if alt > z5 && alt <= z6, H = H5; end
% if alt > z6, H = H6; end
% 




% % % 
% % % % define pressure "base" levels:
% % % 
% % % p0 = 1013.25; % (hPa) (sea level)
% % % p1 = 226.3210; 
% % % p2 = 54.7489;
% % % p3 = 8.6802;
% % % p4 = 1.1091;
% % % p5 = 0.6694;
% % % p6 = 0.0396;
% % % 
% % % % and their corresponding altitude levels (km)
% % % 
% % % z0 = 0;
% % % z1 = 11;
% % % z2 = 20;
% % % z3 = 32;
% % % z4 = 47;
% % % z5 = 51;
% % % z6 = 71;
% % % 
% % % % now calculate scale heights...
% % % 
% % % H0 = - (z1-z0) / (ln(p1/p0));
% % % H1 = - (z2-z1) / (ln(p2/p1));
% % % H2 = - (z3-z2) / (ln(p3/p2));
% % % H3 = - (z4-z3) / (ln(p4/p3));
% % % H4 = - (z5-z4) / (ln(p5/p4));
% % % H5 = - (z6-z5) / (ln(p6/p5));
% % % 
% % % H6 = H5;
% % % 
% % % 
% % % if alt <= z1, H = H0; end
% % % if alt > z1 && alt <= z2, H = H1; end
% % % if alt > z2 && alt <= z3, H = H2; end
% % % if alt > z3 && alt <= z4, H = H3; end
% % % if alt > z4 && alt <= z5, H = H4; end
% % % if alt > z5 && alt <= z6, H = H5; end
% % % if alt > z6, H = H6; end
% % % 
% % % 
% % % 












