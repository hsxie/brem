% Hua-sheng XIE, huashengxie@gmail.com, 2024-03-31 19:17
% Merge e-i Brem radiation factors similar to Hoof15
% Ref: Greene59 (non-rel), Chluba20 (reletivistic)
% BH Z=1-36, and non-rel Sommerfeld

close all; clear; clc;

load('WEei_NR_-8.3_3_0.05.mat');
load('WEei_BH_-6_3_0.02.mat');

Zmax=36;

ekmin=-5; ekmax=3; dek=0.01;
ekm=10.^(ekmin:dek:ekmax); % mc^2
pstr0=[num2str(ekmin),'_',num2str(ekmax),'_',num2str(dek)];

weimZ=zeros(length(ekm),Zmax);
for Z=1:Zmax
    wbh=wbhZ(:,Z);
    wsf=interp1(ekk*Z^2,ws,ekm); wsf(isnan(wsf))=1;
    wbhf=interp1(erkk,wbh,ekm);
    dwf=wbhf-wsf;
    ind=find(dwf==min(dwf));
    % wmf0=[wsf(1:(ind-1)),sqrt(wbhf(ind)*wsf(ind)),wbhf((ind+1):end)];

    fw=0.5*(1+tanh(-((log10(ekm)-log10(ekm(ind)))/0.4))); % width of weight
    % wmf=fw.*wsf+(1-fw).*wbhf;
    wmf=(wsf.^fw).*(wbhf.^(1-fw));
    wmf(isnan(wmf))=1;
    weimZ(:,Z)=wmf;
end

pstr1='mergedWEei_';
save([pstr1,pstr0,'.mat'],'ekm','weimZ');

figure('unit','normalized','DefaultAxesFontSize',16,...
    'DefaultAxesFontWeight','bold','DefaultAxesLineWidth',2,...
    'position',[0.01,0.05,0.55,0.5]);

dist=3;
cmap = jet(floor(Zmax/dist)+1);
j=0;
for Z=1:dist:Zmax
    j=j+1;
    semilogx(ekm,weimZ(:,Z),'-','Color', cmap(j, :),'LineWidth',2); hold on;
    lgdstr{j}=['Z=',num2str(Z),', BH'];
end
legend(lgdstr,'Location','best','FontSize',8);
legend('boxoff');

xlabel('E_{k0}/m_ec^2'); ylabel('W(E_{k0})');
ylim([0,5]);
% 
% 
% semilogx(ekm,wmf,ekm,wmf0,'--',ekm,fw,ekm,dwf,':','LineWidth',2);
% legend('Merged','Sommerfeld+BH','Weight','BH-Sommerfeld dist','Location','best');
% legend('boxoff');
% xlabel('E_{k0}/m_ec^2'); ylabel('W(E_{k0})');
% ylim([0,3]);

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits',...
    'Inches','PaperSize',[screenposition(3:4)]);

print(gcf,'-dpng',[pstr1,pstr0,'.png']);

