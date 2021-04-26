

%% convert list of latitudes to degrees north/south strings


function labels = latlabels(latrange)

degs = repmat(symbol('deg'),length(latrange),1);

NS = repmat('N',length(latrange),1);

NS(latrange < 0) = 'S';

l = [num2str(abs(latrange)') degs NS];

labels = cell(1,size(l,1));

for i = 1:size(l,1)
    labels{i} = strtrim(l(i,:));
end

if any(latrange == 0)
labels{latrange == 0} = ['0' symbol('deg')];
end









