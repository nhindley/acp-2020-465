

% nph_tukeywin

%         ___________
% __---               -------__________
% aaaaaaa-bbbbbbbbbbb-ccccccccccccccccc


% a = length of first section of cosine
% b = length of section equal to unity
% c = length of last section of cosine


function win = nph_tukeywin(a,b,c)

% if isodd(a), alen = 2*a + 1;
% else
% alen = 2*a;
% end
% if iseven(c), clen = 2*c + 1; end


startcosine = tukeywin((2*a)+1,1);
startcosine = startcosine(1:a+1);

endcosine = tukeywin((2*c)+1,1);
endcosine = endcosine(c+1:end);

mid = ones(b,1);

win = [startcosine(1:end-1) ; mid ; endcosine(2:end)];













