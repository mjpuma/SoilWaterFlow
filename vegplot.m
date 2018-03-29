clear all;

load veg_matrix.txt
veg = reshape(veg_matrix,360,180);
veg = veg.*(veg>0);

contourf
% % NOAA/NASA Pathfinder AVHRR SST product
% % http://podaac.jpl.nasa.gov/sst/
% 
% [P,map]=imread('../m_mapWK/199911h54ma-gdm.hdf');
% 
% % Documentation for the 54km dataset gives
% % this formula for temperature
% P=0.15*double(P)-3; % deg C
% 
% %...and defines this Lat/Long grid for the data
% Plat=90-.25-[0:359]*.5;Plon=-180+.25+[0:719]*.5;
% 
% % Since the grid is rectangluar in lat/long (i.e. not
% % really a projection at all, althouhg it is included in
% % m_map under the name 'equidistant cyldindrical'), we
% % don't want to use the 'image' technique. Instead...
% % Create a grid, offsetting by half a grid point to account
% % for the flat pcolor
% [Plg,Plt]=meshgrid(Plon-0.25,Plat+0.25);
% 
% m_proj('hammer-aitoff','clongitude',-150);
% 
% % Rather than rearranging the data so its limits match the
% % plot I just draw it twice (you can see the join at 180W
% % because of the quirks of flat pcolor) (Note that
% % all the global projections have 360 deg ambiguities)
% m_pcolor(Plg,Plt,P);shading flat;colormap(map);
% hold on;
% m_pcolor(Plg-360,Plt,P);shading flat;colormap(map);
% 
% m_coast('patch',[.6 1 .6]);
% m_grid('xaxis','middle');
% 
% % add a standard colorbar.
% h=colorbar('h');
% set(get(h,'title'),'string','AVHRR SST Nov 1999');