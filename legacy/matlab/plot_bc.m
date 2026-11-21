function bcp = plot_bc(bcout,totrain,lambda,alpha,Tseas,tday)
% function to plot top boundary condition
%
lamat = lambda*alpha*Tseas;
figure(2)
hold off
bc2 = bcout([2 4],:)';
sz = size(bc2,1);
t1 = 0;
for i=1:sz
   val = bc2(i,1);
   t2 = bc2(i,2);
   bcp = plot([t1/tday t2/tday],[val val]);
   hold on
   t1 = t2;
end
x1 = xlabel('Day');
y1 = ylabel('rainfall rate [cm/day]');
title(['Total rain = ',num2str(totrain),...
      '; \lambda*\alpha*T = ',num2str(lamat)]);
temp = axis;
text(temp(1,2)*0.6,temp(1,4)*0.8,['alpha = ',num2str(alpha,2),' [depth/storm]']);
text(temp(1,2)*0.6,temp(1,4)*0.72,['lambda = ',num2str(lambda,2),' [storms/day]']);



   
