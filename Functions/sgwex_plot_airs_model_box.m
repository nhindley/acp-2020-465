

% plots the airs/model domain outlines on the 3D plots you've been so
% lovingly making


function sgwex_plot_airs_model_box(airsmodel)

sgbox = sgwexbox;

%%% PLOT BOX OUTLINES =================================================

blon = repmat([sgbox.boxlon sgbox.boxlon(1)],1,2);
blat = repmat([sgbox.boxlat sgbox.boxlat(1)],1,2);

airs_min_height = 18; % km - where do we trust Lars' retrieval to...
airs_max_height = 62; % km

airslinespec = {'linewi',3,'color',[c_color('orange') 0.65],'linest',':'};
airsbalt = [repmat(airs_min_height,1,5) repmat(airs_max_height,1,5)];

modellinespec = {'linewi',3,'color',[c_color('light_blue') 0.65],'linest',':'};
modelbalt = [repmat(1.5,1,5) repmat(75,1,5)];

switch lower(airsmodel)
    
    case 'airs'
        % plot bottom and top:
        hold on; plot3(blon,blat,airsbalt,airslinespec{:});
        % plot sides:
        for i = 1:4
            hold on; plot3([blon(i) blon(i)],[blat(i) blat(i)],[airsbalt(1) airsbalt(end)],airslinespec{:});
        end
        % Box Label:
        text(-45.5,-59,20,'AIRS','fontsize',20,'color',c_color('orange'),'verticalalignment','bottom','horizontalalignment','left');
    case 'model'
        % plot bottom and top:
        hold on; plot3(blon,blat,modelbalt,modellinespec{:});
        % plot sides:
        for i = 1:4
            hold on; plot3([blon(i) blon(i)],[blat(i) blat(i)],[modelbalt(1) modelbalt(end)],modellinespec{:});
        end
        % Box Label
        text(-45.5,-59,3,'MODEL','fontsize',20,'color',c_color('light_blue'),'verticalalignment','bottom','horizontalalignment','left');
end


end












