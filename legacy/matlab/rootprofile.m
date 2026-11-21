zroot = 100;
zevap = 100;
dzconst = 1.0;  % block size

% The distribution of root weights is assigned according  
% to a beta distribution with parameters Ar and Br
nrweight = 0;  % nrweight=1 implies uniform weighting
Ar = 2.7; %1.3; %
Br = 1.1;  %4.0; %
[rweight,rcdfcheck] = betaweight(zroot,dzconst,Ar,Br);
rweight=rweight/rcdfcheck;
% The distribution of evap weights is assigned according 
% to a beta distribution with parameters Ae and Be
neweight = 0;   %neweight==1 implies uniform weighting
Ae = 2;
Be = 2;
[eweight,ecdfcheck] = betaweight(zevap,dzconst,Ae,Be);
eweight=eweight/ecdfcheck;

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
yr = ylabel('Depth/Zroot [-]');
tr = title('ET weights as a function of depth');
tx_leg = text(0.1,-0.1,...
    'Burkea africana - solid, Alternate distribution - dashed');

hold off