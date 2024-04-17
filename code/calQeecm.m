% 2024-03-29 14:39 compare Qcm of Haug89 (11) fiting & integral (A1)
% 20:27 Add Elwert factor and integral low k with (A2) to reduce numerical
% error
% 24-03-31 23:01 add low E with Itoh01, and high E with ER
% 24-04-03 11:50 Update for more accurate integral

close all; clear; clc;

a=1/137;%fine structure constant
r0=2.8179e-15;
mec2=511;

Ekk=10.^(-8:0.05:3);
ee=Ekk+1;
% Ek=(e-1)*mec2;

Q8911=fQcm8911(ee);
Q89A1=fQcm89A1(ee,0); % without Elwert
Q89A1Fee=fQcm89A1(ee,1); % with Elwert

QNR=fQcmNR(ee);
QER=fQcmER(ee);
Qitoh=fQcmitoh(ee,0); % without Elwert
QitohFee=fQcmitoh(ee,1); % with Elwert
Q75b=fQcm75b(ee);

%%
figure('unit','normalized','DefaultAxesFontSize',16,...
    'DefaultAxesFontWeight','bold',...
    'DefaultAxesLineWidth',2,...
    'position',[0.01,0.05,0.85,0.45]);

Qmerge=Q89A1; Qmerge(Ekk<=5e-3)=Qitoh(Ekk<=5e-3);
QmergeFee=Q89A1Fee; QmergeFee(Ekk<=5e-3)=QitohFee(Ekk<=5e-3);
subplot(121);
loglog(Ekk,Q89A1,Ekk,Q89A1Fee,'-.',Ekk,Q8911,'--',Ekk,Q75b,'--',Ekk,Qitoh,':',...
    Ekk,QitohFee,':',Ekk,QNR,':',Ekk,QER,':','LineWidth',2);
% semilogx(ee,(a*r0^2)*Q,ee,QA1,'--','LineWidth',2);
xlabel('E_k/mc^2'); ylabel('Qcm');
legend('Haug89(A1)','Haug89(A1) Fee','Haug89(11)','Haug75b',...
    'Itoh','Itoh Fee','NR','ER','Location','best');
legend('boxoff');
xlim([min(Ekk),max(Ekk)]);

subplot(122);
loglog(Ekk,Q8911./QmergeFee,'-',Ekk,Q75b./QmergeFee,'-','LineWidth',2);
xlabel('E_k/mc^2'); ylabel('Qcm(11)/Qcm(A1)');
legend('Haug89(11)','Haug75b','Location','best');
legend('boxoff');
ylim([0.9,1.1]);

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits',...
    'Inches','PaperSize',[screenposition(3:4)]);

print(gcf,'-dpng',['calQcm_accurate.png']);
save('Haug89Qcmdata_accurate','ee','Ekk','Q89A1','Q89A1Fee','QNR','QER','Qitoh',...
    'QitohFee','Q8911','Q75b','Qmerge','QmergeFee');

%%
function Q=fQcm8911(ee)
mec2=511;

Q=0.*ee;
for j=1:length(ee)
    e=ee(j);
    Ek=(e-1)*mec2;
    % e1=Ek1/mec2+1; % gamma
    p=sqrt(e.^2-1); % normalized by mec2
    p2=p.^2;

    if(Ek<100)
        if(Ek<=10)
            a0=7.31; a1=3.74; a2=1.78; a3=-23.8; a4=-106;
        else
            a0=7.83; a1=-0.27; a2=-0.08; a3=0.24; a4=0.68;
        end
        Q(j)=p2*(a0+a1*p+a2*p^2+a3*p^3+a4*p^4);
    elseif(Ek>=100 && Ek<=10e3)
        b0=8.59; b1=-12.08; b2=-0.738; b3=-0.254; b4=5.77;
        Q(j)=p2/e*(b0+b1*(p/e)+b2*(p/e)^2+b3*(p/e)^3+...
            (32/3+b4*(p/e)^4)*log(e+p));
    elseif(Ek>10e3)
        Q(j)=8*p^2/e*(1-4/3*(p/e)+2/3*(2+(p/e)^2)*log(e+p));
    end
end
Q=real(Q);
end

function Q=fQcmER(ee)
Q=16*ee.*(log(2*ee)-1/6);
end

function Q=fQcmNR(ee)
Q=8*(ee.^2-1);
end

function Q=fQcm75b(ee)
p=sqrt(ee.^2-1);
Q=8*p.^2./ee.*(1-4/3*(p./ee)+2/3*(2+(p./ee).^2).*log(ee+p));
end

function Q=fQcm89A1(ee,iFee)
mec2=511;

Q=0.*ee;
for j=1:length(ee)
    e=ee(j);

    Ek1=(e-1)*mec2;
    % e1=Ek1/mec2+1; % gamma
    p=sqrt(e.^2-1); % normalized by mec2

    % 24-04-03 11:55
    kmax=p^2/e;
    funx=@(k)k.*frdsgmdk89(k,Ek1,iFee);
    funx2=@(k)k.*frdsgmdk89A2(k,Ek1,iFee);
    if((e-1)<1e-3) % 24-04-03 13:20
        Q(j)=NaN;
    else
        if(e-1<1e-2)
            kmid=0.3*kmax; %
            Q(j)=integral(funx2,0,kmid)+integral(funx,kmid,kmax);
        else
            Q(j)=integral(funx,0,kmax);
        end
    end

end
Q=real(Q);
end


function Q=fQcmitoh(ee,iFee)
mec2=511;

Q=0.*ee;
for j=1:length(ee)
    e=ee(j);

    Ek1=(e-1)*mec2;
    % e1=Ek1/mec2+1; % gamma
    p=sqrt(e.^2-1); % normalized by mec2

    kmax=p^2/e;
    funx=@(k)k.*frdsgmdkitoh(k,Ek1,iFee);
    Q(j)=integral(funx,0,kmax);

end
Q=real(Q);
end

function dsdk=frdsgmdkitoh(kk,Ek1,iFee)

a=1/137;%fine structure constant
mec2=511;

% k=hnu/mec2;
dsdk=0.*kk;
for jk=1:length(kk)
    % k=kk(jk)/mec2;
    k=kk(jk);

    e=Ek1/mec2+1; % gamma
    p=sqrt(e^2-1); % normalized by mec2

    if(iFee==1)
        a1=a*(e^2+p^2)/(2*e*p); % 24-03-29 20:50 Elwert factor
        a2=a*(e^2+p^2-2*e*k)/(2*sqrt((e^2-e*k)*(p^2-e*k)));
        Fee=a2/a1*(exp(2*pi*a1)-1)/(exp(2*pi*a2)-1);
    else
        Fee=1;
    end
    Fy=@(y)((17-3*y.^2./(2-y).^2).*sqrt(1-y)+...
        ((12*(2-y).^4-7*(2-y).^2.*y.^2-3*y.^4)./(2-y).^3).*log(1./sqrt(y)+ ...
        sqrt(1./y-1)));

    dsdk(jk)=4/15*Fee/(k)*Fy(k/p^2);
end
end

function dsdk=frdsgmdk89(kk,Ek1,iFee)
% 24-03-29 14:00 should be ok, Haug89(A1)

a=1/137;%fine structure constant
% r0=2.8179e-15;
mec2=511;

% k=hnu/mec2;
dsdk=0.*kk;
for jk=1:length(kk)
    % k=kk(jk)/mec2;
    k=kk(jk);

    e=Ek1/mec2+1; % gamma
    p=sqrt(e^2-1); % normalized by mec2

    if(iFee==1)
        a1=a*(e^2+p^2)/(2*e*p); % 24-03-29 20:50 Elwert factor
        a2=a*(e^2+p^2-2*e*k)/(2*sqrt((e^2-e*k)*(p^2-e*k)));
        Fee=a2/a1*(exp(2*pi*a1)-1)/(exp(2*pi*a2)-1);
    else
        Fee=1;
    end

    L1=log(sqrt(e^2-e*k)+sqrt(p^2-e*k));
    L3=log((e*sqrt(p^2-e*k)+p*sqrt(e^2-e*k))^2/(e*k));

    AA=sqrt((p^2-e*k)/(e^2-e*k))*...
        (1/p*log(e+p)*(2*e+1/e+(e-k)*(1+(4-2*e*k)/p^2+1/p^4+k/(e*p^2+k))...
        +(p^2/(p^2-e*k))*(4*k-e+2*k/(p^2-e*k)))-(e-k)*(16/3*e+2*(2*e-k)/p^2+...
        1/(e*p^4))-4*k^2-k/e-2*p^2/(p^2-e*k))...
        +k/(2*p)*sqrt(p^2-e*k)*log((2*p*sqrt(e^2-e*k)+k)/(2*p*sqrt(e^2-e*k)-k))*...
        ((e/(e*p^2+k))-(e^2+p^2)/(p^2-e*k)-2*e*k/((p^2-e*k)^2))-...
        (e^2*k^2)/(2*p*sqrt((p^2-e*k)*(e^2-e*k)))*log((2*e*p+(e-p)*k)/(2*e*p-(e+p)*k))...
        *(2*(e-k)/(p^2-e*k)+(e^2+1)/(e*p^2+k))+(18-2*k/e+1/e^2)*L1+L3/(e*p)*...
        (32/3*e^3*(e-k)+8*e^2*k^2+k^2-14*e^2-14/3*p^2+19/3*e*k...
        -k/e-k/p^2*(e-k)-k/(2*e*p^4)-(e^2+p^2)/(p^2-e*k)+...
        1/p*log(e+p)*(16*e*p^2-10*p^2*k-4*e-k+11/(2*e)+(4*k-e)/p^2+k/p^4+...
        +(3*(e^2+p^2)-4*e^3*k)/(2*(p^2-e*k))*k+(p^2/2)*((e+k)/((p^2-e*k)^2))-...
        2*((e^2+p^2)^3)/(e^2*k))+...
        1/p*log((2*p^2-(e-p)*k)/(2*p^2-(e+p)*k))*(2*e^3-4*e+k+3/(4*e)...
        +1/(p^2-e*k)*(e^3*k^2-3/4*k*(e^2+p^2)-(p^2/4)*(e+k)/(p^2-e*k))))+...
        2*L1*L3/(p*sqrt((e^2-e*k)*(p^2-e*k)))*(6*p^2*k-16*e*p^2-e*k^2-2*e-3/e+...
        (e^2+p^2)^3/(e^2*k));

    % 24-04-03 12:10
    fkap1=@(t) k*(e-p*t);
    fkap2=@(t) k*(e+p*t);

    fR1=@(t) 4*(p^2-e*k)+4*fkap1(t)+fkap1(t).^2/(e^2-e*k);
    fW2=@(t) sqrt(e*fkap2(t)/(e-k).*(2*k^2+(p^2+k^2-2*e*k)*fkap2(t)));
    fW4=@(t) 2*p*sqrt(4*p^2*(p^2-e*k)+fkap1(t).*fkap2(t)/(e^2-e*k));
    fL2=@(t) log(1+2*(e-k)./(k*fkap1(t)).*((p^2-e*k).*fkap2(t)+sqrt(p^2-e*k).*fW2(t)));
    fL4=@(t) log(1+2*(e^2-e*k)./(fkap1(t).*fkap2(t)).*(4*p^2*(p^2-e*k)+sqrt(p^2-e*k).*fW4(t)));
    fL=@(t) log(sqrt(e^2-e*k)./fkap1(t).*(2*p^2-fkap2(t)+sqrt((p^2-e*k)*fR1(t))));

    fBB=@(t) (2*fL(t)./sqrt(fR1(t))).*...
        ((2-e^2)./fkap2(t)+2./fkap1(t).^2+(e*p^2)./(fkap1(t).*fkap2(t))*(4*e-k)+...
        (4*p^2-3*e*k+1./fkap1(t))./(fkap1(t).*(2*p^2-fkap2(t)))+e./(fkap2(t).*fR1(t)).*...
        (e*fkap1(t)-(e-k)*fkap2(t)).*((2*p^2)./fkap1(t)+(e^2+1)/(e^2-e*k)))...
        +(fL2(t)./fW2(t)).*(2*(e^2+p^2)./fkap2(t)-3/2+ ...
        (e^2+p^2-2*e*k)/(e*k).*(e^2-e*k+fkap1(t)/4)+...
        p^2./(2*p^2-fkap2(t))+1./(fkap1(t).*(2*p^2-fkap2(t))).*(...
        (e^2+p^2-2*e*k)*((e^2+p^2-e*k).*fkap2(t)-(e^2+p^2))...
        -fkap2(t)-fkap2(t).^2/2)+1./(4*e*k*(2*p^2-fkap2(t))).*(...
        2*(e^2+p^2-2*e*k)*(e*k+1-4*p^2*(p^2-e*k))+...
        fkap1(t).^2-3*fkap1(t)+(e^2+p^2).*fkap2(t)))...
        -(fL4(t)./fW4(t)).*(1-e*k/(2*p^2)-1/(4*p^4)+ ...
        1./(fkap1(t).*fkap2(t)).*(1-fkap1(t)-fkap1(t).^2)+... %
        (e^2+p^2-2*e*k)./(fkap1(t).*fkap2(t))*(4*p^2*(p^2-e*k)+e^2*(2*k^2-2+e*k/p^2)));

    fdsgm=@(t) fBB(t)*k^2/sqrt(e^2-e*k);
    dsdk(jk)=Fee/(e*p*k)*(AA+integral(fdsgm,-1,1));

end
end

function dsdk=frdsgmdk89A2(kk,Ek1,iFee)
% 24-03-29 14:00 should be ok, Haug89(A2)

a=1/137;%fine structure constant
% r0=2.8179e-15;
mec2=511;

% k=hnu/mec2;
dsdk=0.*kk;
for jk=1:length(kk)
    % k=kk(jk)/mec2;
    k=kk(jk);

    e=Ek1/mec2+1; % gamma
    p=sqrt(e^2-1); % normalized by mec2

    fx=@(x)log(e-p*x)./(e+p*x);
    Jep=integral(fx,-1,1);

    if(iFee==1)
        a1=a*(e^2+p^2)/(2*e*p);
        a2=a*(e^2+p^2-2*e*k)/(2*sqrt((e^2-e*k)*(p^2-e*k)));
        Fee=a2/a1*(exp(2*pi*a1)-1)/(exp(2*pi*a2)-1);
    else
        Fee=1;
    end

    dsdk(jk)=Fee/k*(log(4*e*p^2/k)*...
        (32/3+8/(3*e^2*p^2)-(e^2+p^2)^3/(e^4*p^4)+(e^2+p^2)^2/(e^5*p^5)*log(e+p))-...
        16/3+4/e^2-2*e^2/p^4+...
        2/(e*p)*(1+3/e^2+e^4/p^4)*log(e+p)+Jep/(e*p^2)*(1/e^2-4-1/(2*e^2*p^2))+...
        k*(log(4*e*p^2/k)*...
        (1/e*(-32/3+1/(6*p^2)+1/(e^2*p^4)-3*e^4/(4*p^6)-3*p^2/(4*e^4))...
        +1/(4*p)*log(e+p)*(1/p^4-1/e^4+3/p^6+3/e^6))-16*p^2/(3*e^3)-...
        2/(e*p^2)+2/(3*e*p^4)-1/(2*e^5)-1/(e*p^6)...
        +1/(e^2*p)*log(e+p)*(-2+3/e^2-1/(e^2*p^2)-1/(2*e^4)+...
        1/(2*e^2*p^4)+1/p^6)-Jep/(2*e^2*p^4)*(1+3/(4*e^2)+3/(4*p^2)))...
        +k^2*(log(4*e*p^2/k)*(49/(4*p^2)+1/(6*p^4)-29/(24*p^6)+...
        7/(8*p^8)-17/(4*e^2)-37/(12*e^4)+5/(8*e^6)+1/(e*p)*log(e+p)*...
        (7/(4*e^2)-1/(2*e^4)+5/(8*e^6)-39/(4*p^2)+1/(4*p^4)+5/(8*p^6)-...
        7/(8*p^8)))+4/(3*e^2)+7/(2*e^4)-9/(16*e^6)-29/(48*e^2*p^4)-...
        5/(4*e^4*p^4)-6/p^6-31/(16*e^2*p^8)...
        +1/(e*p)*log(e+p)*(2*e^4/p^8-e^2/(16*p^8)+19/(8*p^6)+4*e^2/p^4+...
        9/(4*e^4*p^4)-2/e^2-9/(16*e^6))+Jep/(e*p^2)*(2-33/(32*p^2)-...
        9/(4*p^4)-23/(32*p^6)-15/(32*e^2)+9/(32*e^4))));

end
end