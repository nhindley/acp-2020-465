function clim(axes_handle,lims)

if nargin == 1,
    lims = axes_handle;
    set(gca,'clim',lims);
else
    
    set(axes_handle,'clim',lims);
    
end



end



