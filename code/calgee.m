% 2024-03-28 12:10 Calculate Haug89 e-e Brem radiation power
% 24-04-03 14:56 cal Maxwellian gee from Qcm

close all; clear; clc;

a=1/137;%fine structure constant
r0=2.8179e-15;
mec2=511;
if(1==1)
    figure('unit','normalized','DefaultAxesFontSize',16,...
        'DefaultAxesFontWeight','bold',...
        'DefaultAxesLineWidth',2,...
        'position',[0.01,0.05,0.75,0.45]);

    load('Haug89Qcmdata_accurate.mat');

    tt=10.^(-4:0.05:2);
    gee=0.*tt;
    gee8911=0.*tt;
    geewo=0.*tt;
    gee75b=0.*tt;
    geeitoh=0.*tt;
    geeitohwo=0.*tt;
    geeer=0.*tt;

    for j=1:length(tt)
        t=tt(j);
        ct=3/4*sqrt(3/(2*pi))./sqrt(t);
        fg8911=@(e)e.^2.*(e.^2-1).*besselk(2,2*e/t,1).*exp(-2*(e-1)/t).*fQcm8911(e);
        gee8911(j)=ct*integral(fg8911,1,Inf)/(t*(besselk(2,1/t,1))^2);

        fg=@(e)e.^2.*(e.^2-1).*besselk(2,2*e/t,1).*exp(-2*(e-1)/t).*fQcm(e,Ekk,QmergeFee);
        gee(j)=ct*integral(fg,1,Inf)/(t*(besselk(2,1/t,1))^2); % with Elwert

        fgwo=@(e)e.^2.*(e.^2-1).*besselk(2,2*e/t,1).*exp(-2*(e-1)/t).*fQcm(e,Ekk,Qmerge);
        geewo(j)=ct*integral(fgwo,1,Inf)/(t*(besselk(2,1/t,1))^2); % with out Elwert

        fg75b=@(e)e.^2.*(e.^2-1).*besselk(2,2*e/t,1).*exp(-2*(e-1)/t).*fQcm75b(e);
        gee75b(j)=ct*integral(fg75b,1,Inf)/(t*(besselk(2,1/t,1))^2);

        fgitoh=@(e)e.^2.*(e.^2-1).*besselk(2,2*e/t,1).*exp(-2*(e-1)/t).*fQcm(e,Ekk,QitohFee);
        geeitoh(j)=ct*integral(fgitoh,1,Inf)/(t*(besselk(2,1/t,1))^2); % Itoh01, with Elwert

        fgitohwo=@(e)e.^2.*(e.^2-1).*besselk(2,2*e/t,1).*exp(-2*(e-1)/t).*fQcm(e,Ekk,Qitoh);
        geeitohwo(j)=ct*integral(fgitohwo,1,Inf)/(t*(besselk(2,1/t,1))^2); % Itoh01, without Elwert

        fger=@(e)e.^2.*(e.^2-1).*besselk(2,2*e/t,1).*exp(-2*(e-1)/t).*fQcm(e,Ekk,QER);
        geeer(j)=ct*integral(fger,1,Inf)/(t*(besselk(2,1/t,1))^2); % ER
    end
    % gt=fgt_xie(tt).*tt;
    geenozawa = fgtnozawa(tt);
    geeitohnr=fgitoh01(tt);

    subplot(121);
    loglog(tt,gee,tt,gee8911,':',tt,geenozawa,'--',...
        tt,geewo,'-.',tt,gee75b,'-.',tt,geeitohnr,':','LineWidth',2);
    xlabel('t'); ylabel('gee');
    legend('Haug89(A1)','Haug89(11)','Nozawa09',...
        'Haug89(A1) w/o Elwert','Haug75b','Itoh01','Location','best');
    legend('boxoff');
    xlim([min(tt),max(tt)]);

    subplot(122);
    loglog(tt,geenozawa./gee,'-',tt,gee8911./gee,'-',...
        tt,geewo./gee,'-.',tt,gee75b./gee,'-.','LineWidth',2);
    xlabel('t'); ylabel('gee/gfit');
    legend('Nozawa09/Haug89(A1)','Haug89(11)/Haug89(A1)',...
        'Haug89(A1)wo/Haug89(A1)','Haug75b/Haug89(A1)','Location','best');
    legend('boxoff');
    xlim([min(tt),max(tt)]);

    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits',...
        'Inches','PaperSize',[screenposition(3:4)]);

    print(gcf,'-dpng',['cmpgee.png']);
    % print(gcf,'-dpdf',['cmpgee.pdf']);

    save('gee_accurate.mat','tt','gee','geewo','gee8911','gee75b',...
        'geenozawa','geeitohwo','geeitoh','geeer','geeitohnr');
end
%%
% g=fgitoh01(1e-4)

%%
% ee=10.^(-3:0.01:2);
% Q=fQcm(ee);
% loglog(ee,Q);

function g=fgitoh01(tt)

g=0.*tt;
for j=1:length(tt)
    t=tt(j);
    a=1/137.036;
    fI=@(y,u) u.^2./sqrt(1-y).*(exp(pi*a*sqrt(y./(t.*u)))-1)./(exp(pi*a*sqrt(y./(t.*u.*(1-y))))-1 ...
        ).*(exp(-u./y)./y.^3).*((17-3*y.^2./(2-y).^2).*sqrt(1-y)+...
        ((12*(2-y).^4-7*(2-y).^2.*y.^2-3*y.^4)./(2-y).^3).*log(1./sqrt(y)+sqrt(1./y-1)));
    if(t<1e-3)
        It=integral2(fI,0,1,0,1e2,'AbsTol',1e-10,'RelTol',1e-10);
    else
        It=integral2(fI,0,1,0,Inf);
    end

    g(j)=sqrt(3/2)/(10*pi).*It.*t;
end

end


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

function Q=fQcm75b(ee)
p=sqrt(ee.^2-1);
Q=8*p.^2./ee.*(1-4/3*(p./ee)+2/3*(2+(p./ee).^2).*log(ee+p));
end

function Q=fQcm(ee,Ekk,Qin) % accurate Qcm
% 24-04-03 14:41
% mec2=511;
% Ek=(ee-1)*mec2;
Ek=(ee-1);
Q=interp1(Ekk,Qin,Ek,'pchip');
Q(Ek>max(Ekk))=0; % to ER
Q(Ek<min(Ekk))=0;

end

%%
function gtt = fgtnozawa(tt)

emass = 510.998902;

gtt=0.*tt;
for j=1:length(tt)
    t=tt(j);

    if ((t>=0.05/emass)&& (t<1.0/emass)) % GI
        gt = GI(t);
    elseif((t>=1/emass) && (t<300/emass)) % GII
        gt = GII(t);
    elseif((t>=300/emass) && (t<=7000/emass)) % GIII
        gt = GIII(t);
    elseif (t>7000/emass) % GIV
        gt = GIV(t);
    else
        gt=0.*t;
    end
    gtt(j)=sqrt(3*pi/8)*gt.*t; % 24-04-03 14:35

end

end

%% Function GI
function GI = GI(tau)

pi = 3.141592654;
emass = 510.998902;
b = [2.21564, 1.83879e-1, -1.33575e-1, 5.89871e-2, -1.45904e-2, ...
    -7.10244e-4, 2.80940e-3, -1.70485e-3, 5.26075e-4, 9.94159e-5, -1.06851e-4];
logth = log10(tau);
theta = (logth + 2.65) / 1.35;
jfit = b(1) + b(2) * theta + b(3) * theta^2 + b(4) * theta^3 + ...
    b(5) * theta^4 + b(6) * theta^5 + b(7) * theta^6 + b(8) * theta^7 + ...
    b(9) * theta^8 + b(10) * theta^9 + b(11) * theta^10;
GI = sqrt(8 / (3 * pi)) * jfit;
end

%% Function GII
function GII = GII(tau)

A2 = [0.9217, -13.4988, 76.4539, -217.8301, 320.9753, -188.0667, -82.4161, 163.7191, -60.0248];
A1 = [-9.3647, 95.9186, -397.0172, 842.9376, -907.3076, 306.8802, 291.2983, -299.0253, 76.3461];
A0 = [-37.3698, 380.3659, -1489.8014, 2861.4150, -2326.3704, -691.6118, 2853.7893, -2040.7952, 492.5981];
B1 = [-8.6991, 63.3830, -128.8939, -135.0312, 977.5838, -1649.9529, 1258.6812, -404.7461, 27.3354];
B0 = [-11.6281, 125.6066, -532.7489, 1142.3873, -1156.8545, 75.0102, 996.8114, -888.1895, 250.1386];
C2 = [-5.7752, 46.2097, -160.7280, 305.0070, -329.5420, 191.0770, -46.2718, 0, 0];
C3 = [30.5586, -248.2177, 874.1964, -1676.9028, 1828.8677, -1068.9366, 260.5656, 0, 0];
C4 = [-54.3272, 450.9676, -1616.5987, 3148.1061, -3478.3930, 2055.6693, -505.6789, 0, 0];
C5 = [36.2625, -310.0972, 1138.0531, -2260.8347, 2541.9361, -1525.2058, 380.0852, 0, 0];
C6 = [-8.4082, 74.7925, -282.9540, 576.3930, -661.9390, 404.2930, -102.2330, 0, 0];

AA = zeros(1, 3);
% BB = zeros(1, 2);
% CC = zeros(1, 5);
BB = zeros(1, 3);
CC = zeros(1, 6);

for i = 0:8
    AA(3) = AA(3) + A2(i+1) * tau^(i/8);
    AA(2) = AA(2) + A1(i+1) * tau^(i/8);
    AA(1) = AA(1) + A0(i+1) * tau^(i/8);
    BB(2) = BB(2) + B1(i+1) * tau^(i/8);
    BB(1) = BB(1) + B0(i+1) * tau^(i/8);
    CC(6) = CC(6) + C6(i+1) * tau^(i/6);
    CC(5) = CC(5) + C5(i+1) * tau^(i/6);
    CC(4) = CC(4) + C4(i+1) * tau^(i/6);
    CC(3) = CC(3) + C3(i+1) * tau^(i/6);
    CC(2) = CC(2) + C2(i+1) * tau^(i/6);
end

G0 = 0;
G1 = 0;
G2 = 0;

for Xi = 0:2
    i = floor(Xi);
    for Xj = 2:6
        j = floor(Xj);
        Z = Xi + Xj/8 + 1;
        GAM = gamma(Z);
        G1 = G1 + GAM * AA(i+1) * CC(j);
        G2 = G2 + GAM / Z * BB(i+1) * CC(j);
    end
end

G0 = 2 * AA(3) + AA(2) + AA(1) + 0.5 * BB(2) + BB(1);

GII = G0 + G1 + G2;
end

%% Function GIII
function GIII = GIII(tau)

A2 = [64.7512, -213.8956, 174.1432, 136.5088, -271.4899, 89.3210, 58.2584, -46.0807, 8.7301];
A1 = [49.7139, -189.7746, 271.0298, -269.7807, 420.4812, -576.6247, 432.7790, -160.5365, 23.3925];
A0 = [52.1633, -257.0313, 446.8161, -293.0585, 0, 77.0474, -23.8718, 0, 0.1997];
B1 = [376.4322, -1223.3635, 628.6787, 2237.3946, -3828.8387, 2121.7933, -55.1667, -349.4321, 92.2059];
B0 = [-8.5862, 34.1348, -116.3287, 296.5451, -393.4207, 237.5497, -30.6000, -27.6170, 8.8453];

AA = zeros(1, 3);
BB = zeros(1, 2);

for i = 0:8
    AA(3) = AA(3) + A2(i+1) * tau^(i/8);
    AA(2) = AA(2) + A1(i+1) * tau^(i/8);
    AA(1) = AA(1) + A0(i+1) * tau^(i/8);
    BB(2) = BB(2) + B1(i+1) * tau^(i/8);
    BB(1) = BB(1) + B0(i+1) * tau^(i/8);
end

G0 = 2 * AA(3) + AA(2) + AA(1) + 0.5 * BB(2) + BB(1);

GIII = G0;
end
%% Function GIV
function GIV = GIV(tau)

GUM = exp(0.5772156649);
CON1 = 9 / pi / sqrt(tau);
Gt = log(2 * tau / GUM) + 5 / 4;

GIV = CON1 * Gt;
end
