function [c,f]=elastic_plastic(m,s,a,sgmc,sgm0,sgmi)
%   Elastic-Plastic Case
pr0=sgm0+1/8*m*sgmc-1/8*sqrt(sgmc*(m^2*sgmc+16*sgm0.*m+16*s*sgmc));
C1=m*pr0/sgmc+s;
pe=pr0+(2*(sgm0-pr0)-sgmc*C1^a)/(2+a*m*C1^(a-1));
%sgmi=0.0*pe;
k=1+(sgmc*((m*sgmi/sgmc+s)^a-(m*pe/sgmc+s)^a))/(sgmi-pe);
sgma=(((k-1)/a/m)^(1/(a-1))-s)*sgmc/m;
C=(1-k)*(sgma+sgmi)/2+sgmc*((m*sgmi/sgmc+s)^a+(m*sgma/sgmc+s)^a)/2;
fei=asin((k-1)/(k+1));
c=C*(1-sin(fei))/cos(fei)/2;
f=rad2deg(fei);

end

