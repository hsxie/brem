% Hua-sheng XIE, huashengxie@gmail.com, 2024-02-22 10:35
% Calculate e-i Brem radiation factors
% Ref: Greene59 (non-rel), Chluba20 (reletivistic)
% 24-02-27 12:15 update
% 24-03-31 17:53 output of BH Z=1-36, and non-rel

close all; clear; clc;

iout=1;

ekmin=1; ekmax=3; dek=0.05;
ekk=10.^(ekmin:dek:ekmax); % mc^2

pstr0=[num2str(ekmin),'_',num2str(ekmax),'_',num2str(dek)];
if(iout==1) % non-relativistic

    Z=1; % only for Z=1 is sufficient, NR can be normalized by x=E/Z^2
    a=1/137; % fine structure constant
    b=(Z*a)^2/2; % 24-02-25 11:44

    wk=0.*ekk; wb=0.*ekk; we=0.*ekk; ws=0.*ekk;
    for j=1:length(ekk)

        ek=ekk(j) % mc^2

        % ek for kinetic energy E0
        fk=@(x) 1/ek+0.*x; % Kramers
        fb=@(x) 1/ek*sqrt(3)/pi*log((sqrt(ek)+sqrt(ek-x))./(sqrt(ek)-sqrt(ek-x))); % Born
        fe=@(x) sqrt(3)/pi./sqrt(ek.*(ek-x)).*(1-exp(-2*pi*sqrt(b/ek)))./(1-exp(-2*pi*sqrt(b./(ek-x)))...
            ).*log((sqrt(ek)+sqrt(ek-x))./(sqrt(ek)-sqrt(ek-x))); % Elwert
        fs=@(x) funfs(x,ek,b); % Sommerfeld

        wk(j)=integral(fk,0,ek);
        wb(j)=integral(fb,0,ek);
        we(j)=integral(fe,0,ek);
        ws(j)=integral(fs,0,ek);
        % ws(j)=integral(fs,0,0.999*a,'RelTol',1e-6,'AbsTol',1e-6);
    end
    
    eZ2kk=ekk*Z^2; % eZ2kk=ekk*Z^2
    pstr1='WEei_NR_';
    save([pstr1,pstr0,'.mat'],'eZ2kk','ekk','wk','wb','we','ws');

elseif(iout==2) % non-relativistic

    Zmax=36;
    wbhZ=zeros(length(ekk),Zmax);
    for j=1:length(ekk)
        ek=ekk(j) % mc^2
        for Z=1:Zmax
            % ek for kinetic energy E0
            fbh=@(x) funfbh(x,ek,Z); % Bethe–Heitler
            wbhZ(j,Z)=integral(fbh,0,ek);
        end
    end

    erkk=ekk; % relativistic Ek
    pstr1='WEei_BH_';
    save([pstr1,pstr0,'.mat'],'erkk','wbhZ');
end

%%
close all;
figure('unit','normalized','DefaultAxesFontSize',16,...
    'DefaultAxesFontWeight','bold','DefaultAxesLineWidth',2,...
    'position',[0.01,0.05,0.55,0.5]);
if(iout==1)

    semilogx(ekk,wk,ekk,wb,ekk,we,ekk,ws,':','LineWidth',2);
    legend('Kramers','Born','Elwert','Sommerfeld','Location','best');
    legend('boxoff');
    xlabel('E_{k0}/m_ec^2'); ylabel('W(E_{k0})');
    ylim([0,2]);
    xlim([min(ekk),max(ekk)]);
    % grid on;

elseif(iout==2)

    dist=2;
    cmap = jet(floor(Zmax/dist)+1);
    j=0;
    for Z=1:dist:Zmax
        j=j+1;
        semilogx(ekk,wbhZ(:,Z),'-','Color', cmap(j, :),'LineWidth',2); hold on;
        lgdstr{j}=['Z=',num2str(Z),', BH'];
    end
    legend(lgdstr,'Location','best','FontSize',8);
    legend('boxoff');

    xlabel('E_{k0}/m_ec^2'); ylabel('W(E_{k0})');
    ylim([0,10]);
    % grid on;
end

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits',...
    'Inches','PaperSize',[screenposition(3:4)]);

print(gcf,'-dpng',[pstr1,pstr0,'.png']);

%%
function fs=funfbh(x,ek,Z) % Bethe–Heitler cross-section
% ek=gamma0-1, is kinetic energy, normalized by mc^2
% x=h\nu, photon energy, normalized by mc^2
a=1/137; % fine structure constant

fs=0.*x;
for j=1:length(x)
    xj=x(j);
    gam0=ek+1;
    gamf=gam0-xj;
    p0=sqrt(gam0^2-1);
    pf=sqrt(gamf^2-1);
    eta0=a*Z*gam0/p0;
    etaf=a*Z*gamf/pf;
    lmd0=log(gam0+p0);
    lmdf=log(gamf+pf);
    L=log((gam0*gamf+p0*pf-1)/xj);
    gbh=p0*pf/4-3/8*gam0*gamf*(p0/pf+pf/p0)+gam0*gamf*L+...
        3/8*xj*L*((1+gam0*gamf/p0^2)*lmd0/p0-(1+gam0*gamf/pf^2)*lmdf/pf+...
        xj*(1+gam0*gamf*(gam0*gamf+1)/(p0^2*pf^2)))+...
        3/8*(gamf*pf/p0^2*lmd0+gam0*p0/pf^2*lmdf-2*lmd0*lmdf);

    fs(j)=sqrt(3)/pi/p0^2*etaf/eta0*(1- ...
        exp(-2*pi*eta0))/(1-exp(-2*pi*etaf))*gbh*2; % 24-02-25 15:17
end
fs(isnan(fs))=0;
% fs=real(fs);
fs(imag(fs)~=0)=0;
end

function fs=funfs(x,ek,b) % Sommerfeld
% x=h*nu, ek=E0, b=Z^2*R
% eta0=sqrt(Z^2*R/E0)=sqrt(b/ek);
% etaf=sqrt(Z^2*R/(E0-hv))=sqrt(b/(ek-x));
% y=4*eta0*etaf/(eta0-etaf)^2
% y*d|F|^2/dy
fs=0.*x;
for j=1:length(x)
    eta0=sqrt(b/ek);
    etaf=sqrt(b./(ek-x(j)));
    y=4*eta0.*etaf./(eta0-etaf).^2;
    ydFy2absdy=eta0.*etaf.*abs((hypergeom([-1i*eta0+1,-1i*etaf],1,-y)).^2-...
        (hypergeom([-1i*etaf+1,-1i*eta0],1,-y)).^2)./(etaf-eta0);
    fs(j)=sqrt(3)*pi/ek./(exp(2*pi*eta0)-1)./(1-exp(-2*pi*etaf)).*ydFy2absdy;
end
fs(isnan(fs))=0;
fs(imag(fs)~=0)=0;
end