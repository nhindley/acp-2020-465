

function setfont(value)

axs = findall(gcf,'type','axes');
axs = [axs findall(gcf,'type','polaraxes')]; 

for i = 1:length(axs)
    
    if isa(value,'char'),
        set(axs(i),'fontweight',value);
    end
    
    if isa(value,'double'),
        set(axs(i),'fontsize',value);
    end
    
    
end

end



