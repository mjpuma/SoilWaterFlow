clear all
load sparsat.out;
sat=sparsat;
shortS = sat(1:50000,:); [m,n] = size(shortS);
midtime = floor(m/2);
timesat = zeros(m/2,n);% varsat = zeros(m/2,n);

timesat(1,:)=sat(midtime,:);
%varsat(1,:)=0;
for i=1:(m/2-1)    
   timesat(i+1,:)=mean(sat(midtime-i:midtime+i,:));
  % varsat(i+1,:)=var(sat(midtime-i:midtime+i,:));
end

