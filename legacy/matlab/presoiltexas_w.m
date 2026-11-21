% PREPROCESSOR FOR THE BUCKET MODEL AND THE RICHARDS MODEL
% A.Guswa (2001); M.J.Puma (Jan. - April 2003, Feb 2004)

clear;
%----------------------------------------------------------------------
% BUCKET MODEL PREPROCESSOR 

% Input parameters
%
etmax(1) = 0.442;%0.476;	% maximum daily ET over entire root zone [cm/day]
etwilt(1) = 0.02; %0.013;%   % max. evaporation rate (starts at wilting)
                            % [cm/day]
shygr(1) = 0.14;            % hygroscopic saturation 
swilt(1) = 0.29;% 0.28;%  % wilting point saturation 
sstar(1) = 0.44;% 0.46;%    % saturation at stomata closure 
sfc(1) = 0.56;              % saturation at field capacity 

expn(1) = 4.9;%1/0.378;  % constant b in Campbell model (data from Clapp & 
                    % Hornberger (1978)) or 1/lambda in Brooks and 
                    % Corey model (1964) data from Handbook of Hydrology
rksat(1) = 82.2;    % saturated hydraulic conductivity [cm/day]
poros(1) = 0.43;	% porosity or saturated moisture content

init2 = 0.35;	        % initial saturation or head if iccode==1
zroot = 100.0;%40.0;%      % depth of root zone (need not be an integer 
                        % # of blocks)

tday = 1;           % length of a day in units consistent w/ 
                    % temporal discretization
tmax = 660;		% maximum simulation time
ntmax = 20000000*4;	% maximum number of time steps allowed
dtmax = 2e-03;      % maximum time step size

% Top Boundary Conditions
%
lambda = 0.28;	    oldlambda=lambda; % number of storms per day
meandepth = 1.742;	% mean storm depth (exponential) [L]
lengthcm = 1.0;     % length of a centimeter in units consistent with L
Delta = 0.2;%0.1;  % % interception [L] 0.1 for grass and 0.2 for tree
lambda = lambda*exp(-Delta/meandepth);  % modify lambda to account 
                                        % for interception

% GENERATE OUTPUT FILES (Point Model Input Files)
% Generate storm sequence
%
lambdamod = lambda/tday;  % number of storms per unit time
ntstart = 0;	    % vector of storm start times 
                    % (in terms of time step #)

while max(ntstart)<ntmax & (max(ntstart)*dtmax)<tmax
   delta = round(randexp(1./lambdamod,1,1)/dtmax);
   newtime = delta + max(ntstart);
   ntstart = [ntstart; newtime];
end

% Number of storms, which produce throughfall, that occur
nstorm = size(ntstart,1)-1;

% Create vector of storm depths and subtract vegetation interception
olddepth = randexp(meandepth,2*nstorm,1);
j=1;
for i = 1:nstorm
    test = olddepth(j,1)-Delta;
    while test <= 0
        j=j+1;
        test = olddepth(j,1)-Delta;
    end
    depth(i,1)=olddepth(j,1)-Delta;
    i=i+1; j= j+1;
end

actlambda = nstorm/tmax;
actmeandepth = mean(depth);

% Storm Duration
%
% Case 1: From a Uniform Distribution
%mind = 0.05;     % minimum storm duration [days]
%maxd = 0.15;	 % maximum storm duration [days]
%dur = round((rand(nstorm,1)*(maxd-mind)+mind)/dtmax);
%dur(dur==0)=1;
%
% Case 2: From a Beta Distribution
mind = 0.3/24;       % minimum storm duration [days]
maxd = 6/24;    % maximum storm duration [days]
Abd = 2;          % shape of pdf
Bbd = 4.67;       % based on mean storm duration of 1.5 hours
dur = round((betarnd(Abd,Bbd,[nstorm 1])*(maxd-mind)+mind)/dtmax);
dur(dur==0)=1;

% Prepare SoilBox input
%
boxbc = [0,0; ntstart(2:nstorm+1)*dtmax,depth];
ndays = ceil(boxbc(nstorm+1,1)/tday);
boxbcnew = zeros(ndays,2);
for j=1:ndays
   boxbc(:,1) = boxbc(:,1)-tday;
   boxbcnew(j,:) = [j,sum(boxbc(boxbc(:,1)<0,2))];
   boxbc = boxbc(boxbc(:,1)>=0,:);
end
totrain = sum(boxbcnew(:,2));
Tseas = size(boxbcnew,1);


% INPUT FILES FOR BUCKET MODEL
%  Output for rain.in
%
save rain.in boxbcnew -ascii

%  Output for soil.in
%
SBpar = [poros init2 shygr swilt sstar sfc zroot etmax ...
         etwilt rksat*tday expn];
save soil.in SBpar -ascii

%  Output for rainpar.in
%
rainpar = [meandepth lambda lengthcm tday];
save rainpar.in rainpar -ascii
%
%thsat = poros;
%thzero = shygr*poros;
%thwilt = swilt*poros;
%thstar = sstar*poros;
%int = depth./(dur*dtmax);

%
%----------------------------------------------------------------------
% RICHARDS MODEL PREPROCESSOR

% Input Filenames
%
infile = 'plant1d.in';
%
titlename = 'Prosopis glandulosa - 660 days 05/04';%'Paspalum setaceum';
geofile = 'geometry.in';
controlfile = 'control.in';
matfile = 'material.in';
initfile = 'initial.in';
bcfile = 'bc.in';

% Geometry parameters
%
nblocks = 200;	% number of blocks in soil column
ztop = 0.0;		% depth of top boundary
zevap = 20.0;   % depth of zone of evaporation
%zroot          % See POINT MODEL INPUT 
ngrid = 1;      % if ngrid==1 then dz is quasi-constant
dzconst = 1.0;  % block size
nnotdz = 0;     % number of blocks not dz in size

% The distribution of root weights is assigned according  
% to a beta distribution with parameters Ar and Br
nrweight = 0;  % nrweight=1 implies uniform weighting
Ar = 1.1;    
Br = 2.7;    
[rweight,rcdfcheck] = betaweight(zroot,dzconst,Ar,Br);
rweight=rweight/rcdfcheck;
% The distribution of evap weights is assigned according 
% to a beta distribution with parameters Ae and Be
neweight = 0;   %neweight==1 implies uniform weighting
Ae = 0.9;
Be = 5.0;
[eweight,ecdfcheck] = betaweight(zevap,dzconst,Ae,Be);
eweight=eweight/ecdfcheck;
% Control parameters
%
tstart = 0;		% starting time for the simulation
dt = 0.001;		% initial time step size
%tmax           See POINT MODEL INPUT	
%ntmax          See POINT MODEL INPUT	
%tday           See POINT MODEL INPUT 

etstart = .29;  % start time for evapotranspiration (within a day)
etfin = .79;    % end time for evapotranspiration (within a day)

ipond = 1;      % indicator of what to do under ponded conditions
                % ipond == 0 ; force the water in
                % ipond == 1 ; saturate the surface, excess runs off
                % ipond == 2 ; saturate the surface, excess ponds
%phi = 0;       % angle of column away from vertical

dtwrite = 0.1;  % print interval in terms of time for 
                % integrated saturation
mprint = 100;   % print saturation profile every mprint times 
iscrn = 1;		% flag to write information to screen


imax = 11;		% maxmimum number of nonlinear iterations allowed
itmin = 5;		% in iterations are < itmin, time step is modified
itmax = 11;		% if iterations are > itmax, time step is modified
dtredu = 0.5;	% reduction factor for time step
dtincr = 2.8;	% magnification factor for time step
dtmin = 1e-07;  %minimum time step size
%dtmax          See POINT MODEL INPUT

errpr = 1e-05;	% error tolerance on pressure head
errres = 1e-06; % error tolerance on the residual


% Material parameters
%
nmatl = 1;      % number of different materials (only 1 allowed)
jj(1) = 1;

%poros(1)           See POINT MODEL INPUT 
%expn(1)            See POINT MODEL INPUT  
%rksat(1)           See POINT MODEL INPUT  
alpha(1) = 7.18;%21.8; % entry SUCTION HEAD for Pc-S curve [cm]
                       % value from Clapp and hornberger data within the 
                       % range given in Handbook of Hydrology (p. 5.14)
stor(1) = 1e-10;    % specific storativity
sres(1) = shygr(1);      % value on soil moisture at which all ET = 0
Se(1) = 0.98;       % value of saturation at ~ entry pressure head
%
rtavgdens(1) = 0.02    % NOT RELEVANT IN CURRENT FORMULATION
                       % root density [cm of roots/cm^3 of soil]
frac(1) = 0.5;      % fraction of roots at full saturation 
                    % needed for Tmax
%etmax(1)           See POINT MODEL INPUT 
%etwilt(1)          See POINT MODEL INPUT 
hwilt(1) = -32700;%-45900;%  % plant wilting point [cm]
hstar(1) = -1220;% -920;%  % moisture potential at which 
                    % stomata begin to close [cm]
%swilt(1)           See POINT MODEL INPUT  
%
nparam = 1;		% indicator of how to assign materials
matconst = 1;   % material type to assign to the entire domain
nnotmat = 0;    %
imatl(1) = 0;   % not used unless material is not homogeneous

% Initial conditions
% (if iccode==3, then init1 is the z reference
% and init2 is the reference head)

iccode = 1;	% iccode==1 -> constant initial value
            % iccode==3 -> hydrostatic distribution
isat = 1;	% initial value is in terms of saturation 
            % if isat==1
init1 = 0;	% number of initial values different from 
            % value if iccode==1
%init2      See POINT MODEL INPUT
 
% Top Boundary condition  (ALSO SEE POINT MODEL INPUT)
%
%lday = 1;	    % length of a day in units consistent with dt
fracramp = 0.1; % fraction of the storm duration 
                % for each ramp [-]
nbctop = 2;		% indicates type of b.c. (2=flux)
isattop = 1;	% irrelevant for flux b.c.


% Bottom bounday condition
%
nbotbc = 1;	    % number of different bottom bcs
nbcbot = 1;		% fixed condition
bcbot = swilt(1);    % for comparison with Burkea africana simulations 
isatbot = 1;	% 0=pressure head condition, 
                % 1=saturation condition
timeb = tmax;	% end time for bottom condition


% Prepare Richards Model Input
%
int = depth./(dur*dtmax);
bcchange = zeros(2*nstorm,2);

actmeandur = mean(dur*dtmax);
actmeanint = actmeandur/actmeandepth;

for j=1:nstorm
   ramp = ceil(fracramp*dur(j));

   bcchange(4*j-3,1) = ntstart(j+1);
   bcchange(4*j-2,1) = ntstart(j+1)+ramp;
   bcchange(4*j-1,1) = ntstart(j+1)+dur(j);
   bcchange(4*j,1) = ntstart(j+1)+dur(j)+ramp;
   
   bcchange(4*j-3,2) = int(j)/2;
   bcchange(4*j-2,2) = int(j)/2;
   bcchange(4*j-1,2) = -int(j)/2;
   bcchange(4*j,2) = -int(j)/2;
end

clear bc
bc(1,:) = [0 0.0];
j = 1;

while bc(j,1)<max(bcchange(:,1)) 
   nextntime = min(bcchange(bcchange(:,1)>bc(j,1),1));
   newbc = bc(j,2) + sum(bcchange(bcchange(:,1)==nextntime,2));
   newbc(newbc<0)=0;
   bc = [bc ; nextntime newbc];
   j = j+1;
end  

ntopbc = size(bc,1)-1;
bcout = [ones(ntopbc,1)*nbctop bc(1:ntopbc,2) ...
      ones(ntopbc,1)*isattop bc(2:(ntopbc+1),1)*dtmax]';
%
bcp = plot_bc(bcout,totrain,lambda,meandepth,Tseas,tday);


%----------------------------------------------------------------------
%  INPUT FILES FOR THE RICHARDS MODEL

%  Output for plant1d.in
%
fid = fopen(infile,'w');
fprintf(fid,titlename);
fprintf(fid,'\n');
fprintf(fid,geofile);
fprintf(fid,'\n');
fprintf(fid,controlfile);
fprintf(fid,'\n');
fprintf(fid,matfile);
fprintf(fid,'\n');
fprintf(fid,initfile);
fprintf(fid,'\n');
fprintf(fid,bcfile);
fprintf(fid,'\n');
status = fclose(fid);

%  Output for geometry.in
%
fidg = fopen(geofile,'w');
fprintf(fidg,'%4d %6.3e %6.3e %6.3e\n',nblocks,ztop,zevap,zroot);
fprintf(fidg,'%4d %4d %4d\n',ngrid,nrweight,neweight);
fprintf(fidg,'%6.3e %4d\n',dzconst,nnotdz);
if nrweight~=1
  nroot = size(rweight,1);
  for jjj = 1:nroot
     fprintf(fidg,'%6.3e\n',rweight(jjj));
 end
end 
if neweight~=1
  nevap = size(eweight,1);
  for jjj = 1:nevap
     fprintf(fidg,'%6.3e\n',eweight(jjj));
 end
end
status = fclose(fidg);


%  Output for control.in
%
fidc = fopen(controlfile,'w');
fprintf(fidc,'%6.3e\n',tstart);
fprintf(fidc,'%6.3e\n',dt);
fprintf(fidc,'%6.3e\n',tmax);
fprintf(fidc,'%8d\n',ntmax);
fprintf(fidc,'%6.3e\n',tday);
fprintf(fidc,'%6.3e\n',etstart);
fprintf(fidc,'%6.3e\n',etfin);
fprintf(fidc,'%1d\n',ipond);
fprintf(fidc,'%6.3e\n',dtwrite);
fprintf(fidc,'%4d\n',mprint);
fprintf(fidc,'%1d\n',iscrn);
fprintf(fidc,'%4d\n',imax);
fprintf(fidc,'%4d\n',itmin);
fprintf(fidc,'%4d\n',itmax);
fprintf(fidc,'%6.3e\n',dtredu);
fprintf(fidc,'%6.3e\n',dtincr);
fprintf(fidc,'%6.3e\n',dtmin);
fprintf(fidc,'%6.3e\n',dtmax);
fprintf(fidc,'%6.3e\n',errpr);
fprintf(fidc,'%6.3e\n',errres);
status = fclose(fidc);


%  Output for material.in
%
fidm = fopen(matfile,'w');
fprintf(fidm,'%4d\n',nmatl);
for i=1:nmatl   
   fprintf(fidm,'%4d %6.3e %6.3e %6.3e %6.3e %6.3e %6.3e %6.3e\n',...
      jj(i),poros(i),alpha(i),expn(i),rksat(i),stor(i),sres(i),Se(i));
   fprintf(fidm,'%6.3e %6.3e %6.3e %6.3e %6.3e %6.3e %6.3e %6.3e\n',...
      frac(i),rtavgdens(i),etmax(i),etwilt(i),hwilt(i),hstar(i),swilt(i),...
      sstar(i));
end
fprintf(fidm,'%1d\n',nparam);
if nparam==1
   fprintf(fidm,'%4d %4d\n',matconst,nnotmat);
else
   for i=1:nblocks
      fprintf(fidm,'%4d\n',imatl(i));
  end
end
status = fclose(fidm);


%  Output for initial.in
%
fidi = fopen(initfile,'w');
fprintf(fidi,'%1d %1d\n',iccode,isat);
fprintf(fidi,'%1d %1d\n',init1,init2);
status = fclose(fidi);


%  Output for bc.in
%
%fidb = fopen(bcfile,'w');
%fprintf(fidb,'%4d %4d\n',ntopbc,nbotbc);
%fprintf(fidb,'%1d %9.5e %1d %9.5e\n',bcout);
%fprintf(fidb,'%1d %9.5e %1d %9.5e\n',nbcbot,bcbot,isatbot,timeb);
%status = fclose(fidb);

% Putput for climate parameters
%climate = [oldlambda lambda actlambda meandepth actmeandepth actmeandur actmeanint]
%save climate climate -ascii;
%----------------------------------------------------------------------
% PLOTS OF EVAPORATION AND ROOT WEIGHTS 
%
figure(9)
rd1 = plot(rweight(1:(zroot/dzconst)),-(dzconst/2:dzconst:zroot)/zroot,'g-');
hold on
ed1 = plot(eweight(1:(zevap/dzconst)),-(dzconst/2:dzconst:zevap)/zevap,'r--');
xmax = max(max(rweight),max(eweight));
axis([0 xmax -1 0]);
xr = xlabel('Weight [-]');
yr = ylabel('Depth/(Zroot or Zevap) [-]');
tr = title('ET weights as a function of depth');
tx_leg = text(0.1,-0.1,...
    'Root weights - solid, Evap weights - dashed');
tx_rcdf = text(0.1,-0.2,['Integral of Root weights = ',...
        num2str(rcdfcheck)]);
tx_ecdf = text(0.1,-0.3,['Intgeral of Evap weights = ',...
        num2str(ecdfcheck)]);
hold off


%----------------------------------------------------------------------
% Run Point Model
%SoilBox_puma;