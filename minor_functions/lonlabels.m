

%% convert list of longitudes to degrees east/west strings


function labels = lonlabels(lonrange)

degs = repmat(symbol('deg'),length(lonrange),1);

EW = repmat('E',length(lonrange),1);

EW(lonrange < 0) = 'W';

l = [num2str(abs(lonrange)') degs EW];

labels = cell(1,size(l,1));

for i = 1:size(l,1)
    labels{i} = strtrim(l(i,:));
end

if any(lonrange == 0)
labels{lonrange == 0} = ['0' symbol('deg')];
end





