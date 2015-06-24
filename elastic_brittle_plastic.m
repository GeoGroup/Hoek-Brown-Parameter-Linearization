function [cp,fp,cr,fr]=elastic_brittle_plastic(m,s,a,sgmc,mr,sr,ar,sgmcr,sgm0,sgmi)
% Elastic-Brittle-Plastic Case
pr0=sgm0+1/8*m*sgmc-1/8*sqrt(sgmc*(m^2*sgmc+16*sgm0.*m+16*s*sgmc));
C1=m*pr0/sgmc+s;
pe=pr0+(2*(sgm0-pr0)-sgmc*C1^a)/(2+a*m*C1^(a-1));
%Peak parameter
kp=1+a*m*(m*pe/sgmc+s)^(a-1);
Cp=((m*pe/sgmc+s)^a)*(sgmc-a*m*pe/(m*pe/sgmc+s));
fei=asin((kp-1)/(kp+1));
cp=Cp*(1-sin(fei))/cos(fei)/2;
fp=rad2deg(fei);
%Residual parameter
kr=1+(sgmcr*((mr*sgmi/sgmcr+sr)^ar-(mr*pe/sgmcr+sr)^ar))/(sgmi-pe);
sgma=(((kr-1)/ar/mr)^(1/(ar-1))-sr)*sgmcr/mr;
C=(1-kr)*(sgma+sgmi)/2+sgmcr*((mr*sgmi/sgmcr+sr)^ar+(mr*sgma/sgmcr+sr)^ar)/2;
fei=asin((kr-1)/(kr+1));
cr=C*(1-sin(fei))/cos(fei)/2;
fr=rad2deg(fei);
end

