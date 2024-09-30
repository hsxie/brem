% Hua-sheng Xie, huashengxie@gmail.com, 2024-09-28 12:35
% Accurate fitting of Bremsstrahlung Radiation Gaunt factor (error <1%)
% t=kB*Te/(me*c^2), P=P0*g,  P0 given in Eq.(38) of [Xie24]
% g=Z*gei+gee=sum(Zi^2*ni/ne*gei(t,Zi))+gee(t)
% Ref: 
% [Xie24] H.S.Xie, Bremsstrahlung Radiation Power in Fusion Plasmas 
% Revisited: Towards Accurate Analytical Fitting, PPCF, 2024
% CB=32*pi/3^1.5*sqrt(2*pi)/(me^1.5*c^3*h)*(qe^2/(4*pi*epsilon0))^3*sqrt(1e3*qe)=4.86e-37
% Pbrem=CB*(ne)^2.*sqrt(Te).*g [W/m^3], Te in keV, ne in m^-3
function [g,gei,gee]=fgfit(t,Z)

iFee=1;
gei=fgeifit(t,Z);
gee=fgeefit(t,iFee);
g=Z*gei+gee;

end

function gei=fgeifit(t,Z)

ceff=2*sqrt(3)/pi;
x = [0.4365    2.3857    0.7952    0.5305    0.3257];
fnr= 1+x(1)*(1-exp(-(x(2)*1e-4./(t/Z^2)).^x(4)))- ...
    (x(1)+(1-1/ceff))*exp(-((t/Z^2)./(x(3)*1e-5)).^x(5));

xx=sqrt(t./(t+1));
cc=9/8*sqrt(6/pi);

cr=[1.4502   -2.6772    2.9998   -0.9198];
fr=cr(1)*xx+cr(2)*xx.^2+...
    cc.*(1+cr(3)*(xx.^1-1)+cr(4)*(xx.^3-1)).*sqrt(t).*(log(2*t+1)+3/2-0.5772);

x0= [57601.4561174080	3.44046808898792	16.8063152455324	0.133254253876019];
fz=(Z/10)*x0(1)*(t/sqrt(Z/10)/1e-2).^x0(2)./(exp(x0(3)*(t/sqrt(Z/10)/1e-2).^x0(4))-1);

gei=ceff*(fnr-1*fz)+fr;

end

function gee=fgeefit(t,iFee)

gnr0=2*sqrt(3)/pi*(3/sqrt(2)*t);
cc=3/4*sqrt(pi);

if(iFee==1)
    Fee=0.5*(tanh(0.602*(log10(t)+5.06))+1);
else
    Fee=1;
end

xx=sqrt(t./(t+0.7));
x=[-0.106    3.347   -2.642]; % to 0.3%
gee= 1./sqrt(t+1).*(1.0+x(1)*xx.^2+...
    cc.*(0.295+x(2)*xx.^2+x(3)*xx.^3).*(log(2*t+1))).*Fee.*gnr0;

end