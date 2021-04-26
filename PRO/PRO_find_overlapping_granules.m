

%% Find 3D AIRS granules (or sequential granule pairs) that completely
% overlap the SGWEX modelling domain box.



% % % % % % % SGWEX BOX:
% % % % % minlon = -42.5; maxlon = -31.71;
% % % % % minlat = -58.55; maxlat = -50.46;
% % % %
% % % % % WARPING-CORRECTED SGWEX BOX:
% % % % % True centre location of box.
% % % % centrelon = -37.1067443847656;
% % % % centrelat = -54.5067496299744;
% % % %
% % % % % Box width ~1200x900, diagonals = 1500km
% % % % % 4 --- 3
% % % % % |     |
% % % % % 1 --- 2
% % % %
% % % % % NOTE: don't use the referenceEllipsoid/referenceSphere functions, I can't
% % % % % seem to get a correct answer with them and i don't know why. Had this
% % % % % trouble at Dstl too.
% % % %
% % % % [boxlat(1),boxlon(1)] = reckon(centrelat,centrelon,km2deg(750),180+45);
% % % % [boxlat(2),boxlon(2)] = reckon(centrelat,centrelon,km2deg(750),180-45);
% % % % [boxlat(3),boxlon(3)] = reckon(centrelat,centrelon,km2deg(750),+45);
% % % % [boxlat(4),boxlon(4)] = reckon(centrelat,centrelon,km2deg(750),-45);
% % % %
% % % %
% % % % % boxlats = [boxlat boxlat(1)];
% % % % % boxlons = [boxlon boxlon(1)];

figure; hold all; grid on;

sgwexbox;

% m_proj('ortho','lat',centrelat,'lon',centrelon,'radius',20)
% m_coast;
% m_grid;



% Allow for granules that just miss out longitudinally:
tolerance = 1;
boxlon([1 4]) = boxlon([1 4]) + tolerance;
boxlon([2 3]) = boxlon([2 3]) - tolerance;



BoundingBox = [boxlon(1) boxlat(1) ; boxlon(2) boxlat(4)];

nph_draw_coastline(BoundingBox,0,'color','b','linewi',2);

hold on; plot(boxlon([1:end 1]),boxlat([1:end 1]),'r');


% On selecting granules: Having been looking around the model box and AIRS
% swaths for some time now, there are generally two times when AIRS does
% an overpass that looks good over South Georgia, which sits on the
% interface between two granules divisions, generally. These occur around
% granules 031-032 and 166-167, or thereabouts due to precession of the
% orbits. It's quite rare that a single granule covers all the box.
%
% Typically, an overpass might have three corners of the modelling box
% covered, at least, with the island itself always being contained within
% the granule. So I think it's safe to take these granules that miss
% sections of the box and 3DST them anyway with zeros for missing values.
%

%% LOAD AIRS 3DST =========================================================

direc = '/Users/neil/data/sgwex/AIRS/all/';
files = dir([direc 'airs*.nc']);

savedirec = '/Users/neil/data/sgwex/AIRS/three_corner_overlapping_granule_pairs/matlab/';

% sort files into ascending order:


for i = 1:length(files)-1
%     for i = 600:663
    
    % Select sequential pairs of granules:
    str1 = strsplit(files(i).name,{'_','.'});
    currentnum = str2double(str1{4});
    
    str2 = strsplit(files(i+1).name,{'_','.'});
    nextnum = str2double(str2{4});
    
    if nextnum ~= currentnum+1
        %         disp('not a pair')
        continue
        %     else
        %         disp({str1{3},str1{4},str2{4}})
    end
    
    % Criteria for granule pair selection:
    % - separately, both granules must each overlap at least one corner of the box.
    % - together, both granules must overlap at least three corners.
    % This places the granule interface inside the sgwex box, and kind of
    % centres the granules over the island.
    
    nc1 = getnet([direc files(i).name]);
    nc2 = getnet([direc files(i+1).name]);
    
    %%
    nc1.ret_temp = permute(nc1.ret_temp,[2 3 1]);
    nc2.ret_temp = permute(nc2.ret_temp,[2 3 1]);
    
    Airs = struct;
    
    Airs.T      = cat(2,nc1.ret_temp,nc2.ret_temp);
    Airs.Lon    = cat(2,  nc1.l1_lon,  nc2.l1_lon);
    Airs.Lat    = cat(2,  nc1.l1_lat,  nc2.l1_lat);
    Airs.Alt = nc2.ret_z;
    
    % Lars' Airs time is seconds since 2000?
    starttimeindayssince2000 = nc1.l1_time(1,1)/60/60/24;
    Airs.StartTime = datenum('01-Jan-2000') + starttimeindayssince2000;
    
    
    % Filenames
    Airs.Files = {files(i).name, files(i+1).name};
    
    % granule 1 coords
    lons1 = [nc1.l1_lon(1:90,1)' nc1.l1_lon(90,1:135) ...
        nc1.l1_lon(90:-1:1,135)' nc1.l1_lon(1,135:-1:1)];
    lats1 = [nc1.l1_lat(1:90,1)' nc1.l1_lat(90,1:135) ...
        nc1.l1_lat(90:-1:1,135)' nc1.l1_lat(1,135:-1:1)];
    
    % granule 2 coords
    lons2 = [nc2.l1_lon(1:90,1)' nc2.l1_lon(90,1:135) ...
        nc2.l1_lon(90:-1:1,135)' nc2.l1_lon(1,135:-1:1)];
    lats2 = [nc2.l1_lat(1:90,1)' nc2.l1_lat(90,1:135) ...
        nc2.l1_lat(90:-1:1,135)' nc2.l1_lat(1,135:-1:1)];
    
    % granule pair coords
    pairlons = [Airs.Lon(1:end,1)' Airs.Lon(end,1:end) ...
        Airs.Lon(end:-1:1,end)' Airs.Lon(1,end:-1:1)];
    pairlats = [Airs.Lat(1:end,1)' Airs.Lat(end,1:end) ...
        Airs.Lat(end:-1:1,end)' Airs.Lat(1,end:-1:1)];
    %%
    if any(inpolygon(boxlon,boxlat,lons1,lats1)) ...
            && any(inpolygon(boxlon,boxlat,lons2,lats2))
        
        if sum(inpolygon(boxlon,boxlat,pairlons,pairlats)) >= 3
            
            % Plot? with interface line shown in dots
            hold on; p1 = plot(lons1,lats1,'linestyle',':');
            hold on; p2 = plot(pairlons,pairlats);
            
            
            % Handle time for the save string:
            % We bin overpasses into the 0300 and 1700 overpasses
            
            hours = round((Airs.StartTime - floor(Airs.StartTime)) * 24);
            
            if hours >= 1 && hours <= 5
                hours = 3;
            end
            if hours >= 15 && hours <= 19
                hours = 17;
            end
            
            % Use datestr to handle zeros and shit:
            roundedtime = floor(Airs.StartTime) + (hours/24);
            
            str = datestr(roundedtime,'yyyymmdd-HHMM');
            disp({Airs.Files,str})
            
            text(min(pairlons),min(pairlats),str,'color',p2.Color)
            drawnow;
            
            % SAVE
            % =========================================================
            % =========================================================
            
            savefilename = [str '_AirsSG.mat'];
            
            save([savedirec savefilename],'Airs');
            
            disp(['Saved "' savefilename '"'])
            
            % Copy original data file pairs to save directory using system command:
            eval(['! cp ' direc files(i).name ' ' savedirec files(i).name])
            eval(['! cp ' direc files(i+1).name ' ' savedirec files(i+1).name])
            
            
        end
        
    end
    
end



return





% Find out if they overlap the box:
%     (at least three corners overlapping)

%     if sum(inpolygon(boxlon,boxlat,airslons,airslats)) >= 3
%
%         hold on; plot(airslons,airslats,'color','k');
%         drawnow;
%
% %     else
%
%         continue
%
%         title(strrep(files(i).name,'_','-'));
%
%         drawnow;
%
%         disp(files(i).name)
%
%         % save([direc '../overlapping_swaths/fullyoverlapping/' files(i).name],'Airs')
%
%         savestr = [direc '../overlapping_swaths/fullyoverlapping/' currentstr(1:end-4) '-' nextstr(end-6:end-4) '.mat'];
%
%         save(savestr,'Airs');
%
%         disp(['Saved granules ' currentstr ' and ' nextstr])
%
%
%
% %     end
%
%
% end

%         % If it's a sequential pair, do stuff:
%         if nextnum == currentnum + 1,
%
%             T1 = Airs.T;
%             Lat1 = Airs.Lat;
%             Lon1 = Airs.Lon;
%             alt = Airs.Alt;
%
%             load([direc nextstr])
%
%             T2 = Airs.T;
%             Lat2 = Airs.Lat;
%             Lon2 = Airs.Lon;
%
%             % Concatenate in along-track direction:
%
%             Airs = struct;
%
%             Airs.T = cat(2,T1,T2);
%             Airs.Lat = cat(2,Lat1,Lat2);
%             Airs.Lon = cat(2,Lon1,Lon2);
%             Airs.Alt = alt;
%
%             %             T1 = Airs.T;
%             %             Lat1 = Airs.Lat;
%             %             Lon1 = Airs.Lon;
%             %
%             %             Airs = struct;
%             %             Airs.T = T1;
%             %             Airs.Lat = Lat1;
%             %             Airs.Lon = Lon1;
%
%
%             sz = size(Airs.T);
%
%             x = [1:sz(1) sz(1)*ones(1,sz(2)) sz(1):-1:1 ones(1,sz(2)) 1];
%             y = [ones(1,sz(1)) 1:sz(2) sz(2)*ones(1,sz(1)) sz(2):-1:1 1];
%
%             airslons = [Airs.Lon(1:sz(1),1)' Airs.Lon(sz(1),1:sz(2)) ...
%                 Airs.Lon(sz(1):-1:1,sz(2))' Airs.Lon(1,sz(2):-1:1)];
%             airslats = [Airs.Lat(1:sz(1),1)' Airs.Lat(sz(1),1:sz(2)) ...
%                 Airs.Lat(sz(1):-1:1,sz(2))' Airs.Lat(1,sz(2):-1:1)];

%     if all(inpolygon([minlon minlon maxlon maxlon],[minlat maxlat maxlat minlat],airslons,airslats)),
%             if all(inpolygon(boxlon,boxlat,airslons,airslats)),
%
%                 hold on; plot(airslons,airslats,'color','k');
%
%                 title(strrep(files(i).name,'_','-'));
%
%                 drawnow;
%
%                 disp(files(i).name)
%
%                 % save([direc '../overlapping_swaths/fullyoverlapping/' files(i).name],'Airs')
%
%                 savestr = [direc '../overlapping_swaths/fullyoverlapping/' currentstr(1:end-4) '-' nextstr(end-6:end-4) '.mat'];
%
%                 save(savestr,'Airs');
%
%                 disp(['Saved granules ' currentstr ' and ' nextstr])
%
%
%
%             end
%
% end
%     end

%     x = [1:sz(1) sz(1)*ones(1,sz(2)) sz(1):-1:1 ones(1,sz(2)) 1];
%     y = [ones(1,sz(1)) 1:sz(2) sz(2)*ones(1,sz(1)) sz(2):-1:1 1];
%
%     lons = [Airs.Lon(1:sz(1),1)' Airs.Lon(sz(1),1:sz(2)) ...
%         Airs.Lon(sz(1):-1:1,sz(2))' Airs.Lon(1,sz(2):-1:1)];
%     lats = [Airs.Lat(1:sz(1),1)' Airs.Lat(sz(1),1:sz(2)) ...
%         Airs.Lat(sz(1):-1:1,sz(2))' Airs.Lat(1,sz(2):-1:1)];
%
%
%     % Check if within SGWEX box:
%
%     if all(inpolygon([minlon minlon maxlon maxlon],[minlat maxlat maxlat minlat],lons,lats)),
%
%         hold on; plot(lons,lats,'color','k');
%
%         title(strrep(files(i).name,'_','-'));
%
%         disp(files(i).name);
%
%         drawnow;
%
%     end
%
%
%
% end
%








% 
% 
% 
% 
% 
% %% Find 3D AIRS granules (or sequential granule pairs) that completely
% % overlap the SGWEX modelling domain box.
% 
% 
% 
% % % % % % % % SGWEX BOX:
% % % % % % minlon = -42.5; maxlon = -31.71;
% % % % % % minlat = -58.55; maxlat = -50.46;
% % % % %
% % % % % % WARPING-CORRECTED SGWEX BOX:
% % % % % % True centre location of box.
% % % % % centrelon = -37.1067443847656;
% % % % % centrelat = -54.5067496299744;
% % % % %
% % % % % % Box width ~1200x900, diagonals = 1500km
% % % % % % 4 --- 3
% % % % % % |     |
% % % % % % 1 --- 2
% % % % %
% % % % % % NOTE: don't use the referenceEllipsoid/referenceSphere functions, I can't
% % % % % % seem to get a correct answer with them and i don't know why. Had this
% % % % % % trouble at Dstl too.
% % % % %
% % % % % [boxlat(1),boxlon(1)] = reckon(centrelat,centrelon,km2deg(750),180+45);
% % % % % [boxlat(2),boxlon(2)] = reckon(centrelat,centrelon,km2deg(750),180-45);
% % % % % [boxlat(3),boxlon(3)] = reckon(centrelat,centrelon,km2deg(750),+45);
% % % % % [boxlat(4),boxlon(4)] = reckon(centrelat,centrelon,km2deg(750),-45);
% % % % %
% % % % %
% % % % % % boxlats = [boxlat boxlat(1)];
% % % % % % boxlons = [boxlon boxlon(1)];
% 
% figure; hold all; grid on;
% 
% sgwexbox;
% 
% % m_proj('ortho','lat',centrelat,'lon',centrelon,'radius',20)
% % m_coast;
% % m_grid;
% 
% 
% 
% % Allow for granules that just miss out longitudinally:
% tolerance = 1;
% boxlon([1 4]) = boxlon([1 4]) + tolerance;
% boxlon([2 3]) = boxlon([2 3]) - tolerance;
% 
% 
% 
% BoundingBox = [boxlon(1) boxlat(1) ; boxlon(2) boxlat(4)];
% 
% nph_draw_coastline(BoundingBox,0,'color','b','linewi',2);
% 
% hold on; plot(boxlon([1:end 1]),boxlat([1:end 1]),'r');
% 
% 
% % On selecting granules: Having been looking around the model box and AIRS
% % swaths for some time now, there are generally two times when AIRS does
% % an overpass that looks good over South Georgia, which sits on the
% % interface between two granules divisions, generally. These occur around
% % granules 031-032 and 166-167, or thereabouts due to precession of the
% % orbits. It's quite rare that a single granule covers all the box.
% %
% % Typically, an overpass might have three corners of the modelling box
% % covered, at least, with the island itself always being contained within
% % the granule. So I think it's safe to take these granules that miss
% % sections of the box and 3DST them anyway with zeros for missing values.
% %
% 
% %% LOAD AIRS 3DST =========================================================
% 
% direc = '/Users/neil/data/AIRS/SGWEX-SGeorgia/';
% files = dir([direc 'airs*.nc']);
% 
% savedirec = '/Users/neil/data/AIRS/SGWEX-SGeorgia/three_corner_overlapping_granule_pairs/';
% 
% % sort files into ascending order:
% 
% 
% for i = 1:length(files)-1
%     % for i = 16:17
%     
%     % Select sequential pairs of granules:
%     str1 = strsplit(files(i).name,{'_','.'});
%     currentnum = str2double(str1{4});
%     
%     str2 = strsplit(files(i+1).name,{'_','.'});
%     nextnum = str2double(str2{4});
%     
%     if nextnum ~= currentnum+1
%         %         disp('not a pair')
%         continue
%         %     else
%         %         disp({str1{3},str1{4},str2{4}})
%     end
%     
%     % Criteria for granule pair selection:
%     % - separately, both granules must each overlap at least one corner of the box.
%     % - together, both granules must overlap at least three corners.
%     % This places the granule interface inside the sgwex box, and kind of
%     % centres the granules over the island.
%     
%     nc1 = getnet([direc files(i).name]);
%     nc2 = getnet([direc files(i+1).name]);
%     
%     %%
%     nc1.ret_temp = permute(nc1.ret_temp,[2 3 1]);
%     nc2.ret_temp = permute(nc2.ret_temp,[2 3 1]);
%     
%     Airs = struct;
%     
%     Airs.T      = cat(2,nc1.ret_temp,nc2.ret_temp);
%     Airs.Lon    = cat(2,  nc1.l1_lon,  nc2.l1_lon);
%     Airs.Lat    = cat(2,  nc1.l1_lat,  nc2.l1_lat);
%     Airs.Alt = nc2.ret_z;
%     
%     % Lars' Airs time is seconds since 2000?
%     starttimeindayssince2000 = nc1.l1_time(1,1)/60/60/24;
%     Airs.StartTime = datenum('01-Jan-2000') + starttimeindayssince2000;
%     
%     
%     % Filenames
%     Airs.Files = {files(i).name, files(i+1).name};
%     
%     % granule 1 coords
%     lons1 = [nc1.l1_lon(1:90,1)' nc1.l1_lon(90,1:135) ...
%         nc1.l1_lon(90:-1:1,135)' nc1.l1_lon(1,135:-1:1)];
%     lats1 = [nc1.l1_lat(1:90,1)' nc1.l1_lat(90,1:135) ...
%         nc1.l1_lat(90:-1:1,135)' nc1.l1_lat(1,135:-1:1)];
%     
%     % granule 2 coords
%     lons2 = [nc2.l1_lon(1:90,1)' nc2.l1_lon(90,1:135) ...
%         nc2.l1_lon(90:-1:1,135)' nc2.l1_lon(1,135:-1:1)];
%     lats2 = [nc2.l1_lat(1:90,1)' nc2.l1_lat(90,1:135) ...
%         nc2.l1_lat(90:-1:1,135)' nc2.l1_lat(1,135:-1:1)];
%     
%     % granule pair coords
%     pairlons = [Airs.Lon(1:end,1)' Airs.Lon(end,1:end) ...
%         Airs.Lon(end:-1:1,end)' Airs.Lon(1,end:-1:1)];
%     pairlats = [Airs.Lat(1:end,1)' Airs.Lat(end,1:end) ...
%         Airs.Lat(end:-1:1,end)' Airs.Lat(1,end:-1:1)];
%     %%
%     if any(inpolygon(boxlon,boxlat,lons1,lats1)) ...
%             && any(inpolygon(boxlon,boxlat,lons2,lats2))
%         
%         if sum(inpolygon(boxlon,boxlat,pairlons,pairlats)) >= 3
%             
%             % Plot? with interface line shown in dots
%             hold on; plot(lons1,lats1,'color','k','linestyle',':');
%             hold on; plot(pairlons,pairlats,'color','k');
%             
%             drawnow;
%             
%             % Handle time for the save string:
%             % We bin overpasses into the 0300 and 1700 overpasses
%             
%             hours = round((Airs.StartTime - floor(Airs.StartTime)) * 24);
%             
%             if hours >= 1 && hours <= 5
%                 hours = 3;
%             end
%             if hours >= 15 && hours <= 19
%                 hours = 17;
%             end
%             
%             % Use datestr to handle zeros and shit:
%             roundedtime = floor(Airs.StartTime) + (hours/24);
%             
%             str = datestr(roundedtime,'yyyymmdd-HHMM');
%             disp({Airs.Files,str})
%             
%             savefilename = [str '_AirsSG.mat'];
%             
%             save([savedirec savefilename],'Airs');
%             
%             disp(['Saved "' savefilename '"'])
%             
%             % Copy original data file pairs to save directory using system command:
%             eval(['! cp ' direc files(i).name ' ' savedirec files(i).name])
%             eval(['! cp ' direc files(i+1).name ' ' savedirec files(i+1).name])
%             
%             
%         end
%         
%     end
%     
% end
% 
% 
% 
% 
% 
% 
% % Find out if they overlap the box:
% %     (at least three corners overlapping)
% 
% %     if sum(inpolygon(boxlon,boxlat,airslons,airslats)) >= 3
% %
% %         hold on; plot(airslons,airslats,'color','k');
% %         drawnow;
% %
% % %     else
% %
% %         continue
% %
% %         title(strrep(files(i).name,'_','-'));
% %
% %         drawnow;
% %
% %         disp(files(i).name)
% %
% %         % save([direc '../overlapping_swaths/fullyoverlapping/' files(i).name],'Airs')
% %
% %         savestr = [direc '../overlapping_swaths/fullyoverlapping/' currentstr(1:end-4) '-' nextstr(end-6:end-4) '.mat'];
% %
% %         save(savestr,'Airs');
% %
% %         disp(['Saved granules ' currentstr ' and ' nextstr])
% %
% %
% %
% % %     end
% %
% %
% % end
% 
% %         % If it's a sequential pair, do stuff:
% %         if nextnum == currentnum + 1,
% %
% %             T1 = Airs.T;
% %             Lat1 = Airs.Lat;
% %             Lon1 = Airs.Lon;
% %             alt = Airs.Alt;
% %
% %             load([direc nextstr])
% %
% %             T2 = Airs.T;
% %             Lat2 = Airs.Lat;
% %             Lon2 = Airs.Lon;
% %
% %             % Concatenate in along-track direction:
% %
% %             Airs = struct;
% %
% %             Airs.T = cat(2,T1,T2);
% %             Airs.Lat = cat(2,Lat1,Lat2);
% %             Airs.Lon = cat(2,Lon1,Lon2);
% %             Airs.Alt = alt;
% %
% %             %             T1 = Airs.T;
% %             %             Lat1 = Airs.Lat;
% %             %             Lon1 = Airs.Lon;
% %             %
% %             %             Airs = struct;
% %             %             Airs.T = T1;
% %             %             Airs.Lat = Lat1;
% %             %             Airs.Lon = Lon1;
% %
% %
% %             sz = size(Airs.T);
% %
% %             x = [1:sz(1) sz(1)*ones(1,sz(2)) sz(1):-1:1 ones(1,sz(2)) 1];
% %             y = [ones(1,sz(1)) 1:sz(2) sz(2)*ones(1,sz(1)) sz(2):-1:1 1];
% %
% %             airslons = [Airs.Lon(1:sz(1),1)' Airs.Lon(sz(1),1:sz(2)) ...
% %                 Airs.Lon(sz(1):-1:1,sz(2))' Airs.Lon(1,sz(2):-1:1)];
% %             airslats = [Airs.Lat(1:sz(1),1)' Airs.Lat(sz(1),1:sz(2)) ...
% %                 Airs.Lat(sz(1):-1:1,sz(2))' Airs.Lat(1,sz(2):-1:1)];
% 
% %     if all(inpolygon([minlon minlon maxlon maxlon],[minlat maxlat maxlat minlat],airslons,airslats)),
% %             if all(inpolygon(boxlon,boxlat,airslons,airslats)),
% %
% %                 hold on; plot(airslons,airslats,'color','k');
% %
% %                 title(strrep(files(i).name,'_','-'));
% %
% %                 drawnow;
% %
% %                 disp(files(i).name)
% %
% %                 % save([direc '../overlapping_swaths/fullyoverlapping/' files(i).name],'Airs')
% %
% %                 savestr = [direc '../overlapping_swaths/fullyoverlapping/' currentstr(1:end-4) '-' nextstr(end-6:end-4) '.mat'];
% %
% %                 save(savestr,'Airs');
% %
% %                 disp(['Saved granules ' currentstr ' and ' nextstr])
% %
% %
% %
% %             end
% %
% % end
% %     end
% 
% %     x = [1:sz(1) sz(1)*ones(1,sz(2)) sz(1):-1:1 ones(1,sz(2)) 1];
% %     y = [ones(1,sz(1)) 1:sz(2) sz(2)*ones(1,sz(1)) sz(2):-1:1 1];
% %
% %     lons = [Airs.Lon(1:sz(1),1)' Airs.Lon(sz(1),1:sz(2)) ...
% %         Airs.Lon(sz(1):-1:1,sz(2))' Airs.Lon(1,sz(2):-1:1)];
% %     lats = [Airs.Lat(1:sz(1),1)' Airs.Lat(sz(1),1:sz(2)) ...
% %         Airs.Lat(sz(1):-1:1,sz(2))' Airs.Lat(1,sz(2):-1:1)];
% %
% %
% %     % Check if within SGWEX box:
% %
% %     if all(inpolygon([minlon minlon maxlon maxlon],[minlat maxlat maxlat minlat],lons,lats)),
% %
% %         hold on; plot(lons,lats,'color','k');
% %
% %         title(strrep(files(i).name,'_','-'));
% %
% %         disp(files(i).name);
% %
% %         drawnow;
% %
% %     end
% %
% %
% %
% % end
% %
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
