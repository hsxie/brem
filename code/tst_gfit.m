% 24-09-30 08:07 compare fitting of gei and gee
close all; clear; clc;

figure('unit','normalized','DefaultAxesFontSize',16,...
    'DefaultAxesFontWeight','bold','DefaultAxesLineWidth',2,...
    'position',[0.01,0.05,0.85,0.5]);

ceff=2*sqrt(3)/pi;

t=10.^(-3:0.01:1);
%%
Z=3;
% gei=fgtz_hoof(t,Z);

g_rider_ei=ceff.*(1+0.7936*t+1.874*t.^2);
g_svensson_ei=ceff.*(1+1.781*t.^1.34);
g_khvesyuk_ei=(1.10+0.59*t+3.06*t.^2-2.56*t.^3+0.85*t.^4);
g_chirkov_ei=ceff.*(0.68+0.32*exp(-4.4*t)+2.07*t+Z*(0.39*(1-exp(-0.008./sqrt(t))-0.49*exp(-505*sqrt(t)))));

%%
% gee = ceff*fgt(t).*t;
% gee = fgt(t);

g_rider_ee=2*sqrt(3)/pi*(3/sqrt(2)*t);
g_svensson_ee=2*sqrt(3)/pi*(5/24/sqrt(2)*(44-3*pi^2)).*t.*(1+1.1*t+t.^2-1.25*t.^2.5);
g_khvesyuk_ee=(1.78*t-0.15*t.^2+0.58*t.^3);
g_chirkov_ee=2*sqrt(3)/pi*(3/sqrt(2)*t).*(1+0.64*t+6.6*t.^2-22.6*t.^3+33.8*t.^4-24.7*t.^5+7.1*t.^6);

gnr=ceff*3/sqrt(2)*t; gnr(t>0.1)=NaN;
% ger=9/8*sqrt(6/pi*t).*(log(2*t)-0.5772+5/4); ger(t<1.0)=NaN;
ger=ceff*9/pi*sqrt(t).*(log(2*t)-0.5772156649+5/4); ger(t<1.0)=NaN;

%%
% g=Z*gei+gee;
[g,gei,gee]=fgfit(t,Z);

g_rider=g_rider_ee+Z*g_rider_ei;
g_svensson=g_svensson_ee+g_svensson_ei*Z;
g_khvesyuk=g_khvesyuk_ee+g_khvesyuk_ei*Z;
g_chirkov=g_chirkov_ee+g_chirkov_ei*Z;
%%

ax1 = axes('Position',[0.07 0.18 0.25 0.75]);
semilogx(t,g_svensson_ei./gei,'-',t,g_rider_ei./gei,'--',...
    t,g_khvesyuk_ei./gei,'--',t,g_chirkov_ei./gei,':','LineWidth',2); hold on;
ylabel('g_{ei}/g_{ei,hoof}'); xlabel('t=k_BT_e/m_ec^2'); xlim([1e-3,1e0]);
legend('Svensson','Rider,McNally','Khvesyuk',...
    'Chirkov','Location','best','Fontsize',12);
legend('boxoff');
text(1e-1,0.8,['Z=',num2str(Z)]);
text(0.6,0.75,'(a)');

ax2 = axes('Position',[0.39 0.18 0.25 0.75]);
semilogx(t,g_svensson_ee./gee,'-',t,g_rider_ee./gee,'--',...
    t,g_khvesyuk_ee./gee,'--',t,g_chirkov_ee./gee,':','LineWidth',2); hold on;
ylabel('g_{ee}/g_{ee,nozawa}'); xlabel('t=k_BT_e/m_ec^2'); xlim([1e-3,1e0]);
legend('Svensson','Rider,Johner,Gould,McNally','Khvesyuk',...
    'Chirkov','Location','best','Fontsize',12);
legend('boxoff');
text(0.6,0.55,'(b)');

ax3 = axes('Position',[0.71 0.18 0.25 0.75]);
semilogx(t,g_svensson./g,'-',t,g_rider./g,'--',...
    t,g_khvesyuk./g,'--',t,g_chirkov./g,':','LineWidth',2); hold on;
ylabel('g/g_{accurate}'); xlabel('t=k_BT_e/m_ec^2'); xlim([1e-3,1e0]);
legend('Svensson','Rider,McNally','Khvesyuk',...
    'Chirkov','Location','best','Fontsize',12);
legend('boxoff');
text(0.6,0.8,'(c)');


set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits',...
    'Inches','PaperSize',[screenposition(3:4)]);

% print(gcf,'-dpng',['cmpgeimodel.png']);
print(gcf,'-dpng',['cmpfitting_Z=',num2str(Z),'.png']);
% print(gcf,'-dpdf',['cmpfitting_Z=',num2str(Z),'.pdf']);
