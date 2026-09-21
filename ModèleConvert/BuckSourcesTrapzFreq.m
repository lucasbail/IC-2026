C1=120e-12;
C2=80e-12;
C3=100e-12;
Cf=1.5e-3;
Lf=10e-6;
Rp=5e-3;
Lp=20e-9;
Rmc=10;

F0=200e3;
Fend=40e6;

%% Tension et ocurant trapezoidales dans le domaine de la fréquence
% tension MC
T0=1/F0;
tm=0.06e-6;
td=0.08e-6;
a=0.52;
V=327;

[~,Vk] = FourierTrapeze(V,T0,a,tm,td,Fend);

% Courant MD
tm=0.06e-6;
td=0.12e-6;
a=0.51;
I=10;
[f,Id] = FourierTrapeze(I,T0,a,tm,td,Fend);

%% Solution du circuit
w=2*pi*f;

x2=zeros(2,length(f));
Zl=1i.*w*Lf;
Zc=Rp*ones(size(w))+1i.*w*Lp+1./(1i.*w*Cf);
for k=1:length(f)
    x2(:,k)=[(-1/Zl(k)-1/Zc(k)-1i*w(k)*C3) (1/Zl(k)+1/Zc(k)) ;...
        (1/Zc(k)+1/Zl(k)) (-1/Zc(k)-1/Zl(k)-1i*w(k)*(C2+C1)) ]...
        \[(Id(k)); (1i*w(k)*C1*Vk(k)-Id(k))];
end
V1=x2(1,:).';
V2=x2(2,:).';
iC1=(Vk+V2)*1i.*w*C1;
iC2=(V2)*1i.*w*C2;
iC3=(V1)*1i.*w*C3;
% V1(k)-s(1)
% V2(k)-s(2)
% iC1(k)-Ic1
% iC2(k)-Ic2
% iC3(k)-Ic3

zero=iC2+iC3+iC1;
disp(max(abs(zero)))

%% Données de simulation sur simplorer
addpath('Données Simplorer')

[f1_MHz,magC1I_A] = importSim('Ic1.csv');
[f2_MHz,magC2I_A] = importSim('Ic2.csv');
[f3_MHz,magC3I_A] = importSim('Ic3.csv');


figure(4)
subplot 311
loglog(f,abs(iC1),'*',f1_MHz*1e6,magC1I_A,'o')
ax=gca;ax.XLim=[f(1)/2 5e7];ax.YLim=[1e-4 1e-1];
title('mag(I_{C1})');ylabel('(A)');xlabel('(Hz)');
legend('Laplace solution','Simulation''s FFT','Location','bestoutside')
ax.YMinorGrid='on';ax.XMinorGrid='on';
subplot 312
loglog(f,abs(iC2),'*',f2_MHz*1e6,magC2I_A,'o')
ax=gca;ax.XLim=[f(1)/2 5e7];
title('mag(I_{C2})');ylabel('(A)');xlabel('(Hz)')
legend('Laplace solution','Simulation''s FFT','Location','bestoutside')
ax.YMinorGrid='on';ax.XMinorGrid='on';
subplot 313
loglog(f,abs(iC3),'*',f3_MHz*1e6,magC3I_A,'o')
ax=gca;ax.XLim=[f(1)/2 5e7];ax.YLim=[1e-4 2e-2];
title('mag(I_{C3})');ylabel('(A)');xlabel('(Hz)')
legend('Laplace solution','Simulation''s FFT','Location','bestoutside')
ax.YMinorGrid='on';ax.XMinorGrid='on';