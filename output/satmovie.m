clear; 
%load plant1d.txt
%plant1d_ba = plant1d; 
%load plant1droot.pl2
load plant1d.pl2
%load retention.pl2

%plant1d = plant1droot;%retention;%
maxdepth = 300
maxtime = 399%999;
water(maxdepth,3) = zeros;
start = 1;

%axis tight
set(gca,'nextplot','replacechildren');

for j = 1:maxtime
   depth(start:maxdepth-1+start) = -plant1d(start:maxdepth-1+start,2);
   sat15(start:maxdepth-1+start) = plant1d(start:maxdepth-1+start,4);
   %satm(start:maxdepth-1+start) = plant1d_ba(start:maxdepth-1+start,4);

    
   water(1:maxdepth,3)=(1-sat15(start:maxdepth-1+start))';
   water(1:maxdepth,2)= 0.1;
   colormap(water)
   
   satur(1:maxdepth) = sat15(start:maxdepth-1+start);
%   saturm(1:maxdepth) = satm(start:maxdepth-1+start);
   
   elevation(1:maxdepth) = depth(start:maxdepth-1+start);
   plot(satur, elevation, 'g')%,saturm, elevation)   
   
   x1(1:maxdepth) = 1; 
   y1(1:maxdepth) = -1:-1:-maxdepth;
   w(1:maxdepth) = 1;
   h(1:maxdepth) = 1;

   for i = 1:maxdepth
      rectangle('Position',[x1(i),y1(i),w(i),h(i)],'FaceColor',...
         [0 water(i,2) water(i,3)],'Edgecolor', [0 water(i,2) water(i,3)])
   end 
   axis([0.8 1 -maxdepth 0]);
   xr = xlabel('Relative Soil Moisture');
   yr = ylabel('Elevation [cm]');
   tr = title('Soil Moisture Profile for alpha = 1.5 cm and lamba = 1/6 days');
   
   %xmax = max(max(rt),max(sat));
	
   satmov(j) = getframe;
   start = start+maxdepth;
end

movie(satmov,1,10)
        

