addpath('D:\ownclound\Thèse\Rassemblement des modèles\Frequence2time')
%%
A=100;
T0=5e-5;
a=0.5;
tm=10e-6;
td=10e-6;
DCb=-0.5;
Fend=500/tm;
f=linspace(0,1/tm,1e5);

[ft,Trap,Cn,Phase] = f_FourierTrapeze(A,T0,a,tm,td,Fend,DCb);
[fa,TrapA,CnA,PhaseA] = f_FourierTrapezeArrondi(A,T0,a,tm,td,Fend,DCb);

t=linspace(0,T0,1e4);
[Trapt2] = f_ifseries(ft,Trap,t);
[TraptA] = f_ifseries(fa,TrapA,t);

Trapt=zeros(size(t));
ind1=t<tm;
Trapt(ind1)=A/tm*t(ind1)+DCb;
ind2=logical((t>=tm).*(t<a*T0+tm/2-td/2));
Trapt(ind2)=A+DCb;
ind3=logical((t>=a*T0+tm/2-td/2).*(t<a*T0+tm/2+td/2));
Trapt(ind3)=A/td*(a*T0+tm/2+td/2-t(ind3))+DCb;
ind4=(t>=a*T0+tm/2+td/2);
Trapt(ind4)=+DCb;

figure(11)
plot(t,Trapt,t,Trapt2,t,TraptA,'--')
xlabel('Time (s)');ylabel('Signal');
legend('Input','Reconstructed','Arrondi');

figure(12)
subplot 211
semilogy(ft,Cn,'.',fa,CnA,'.')
xlabel('Frequency (Hz');ylabel('Magnitude')
legend('Trapz','Trapz arrondi','Location','best')
subplot 212
plot(ft,Phase*180/pi,'.',fa,PhaseA*180/pi)
xlabel('Frequency (Hz');ylabel('Angle (deg)')


figure(13)
p=plot(t,Trapt,t,TraptA);grid on
xlabel('Time (s)');ylabel('Signal');
legend('Input','Arrondi');
Marker=['o';'+'];
for k=1:length(p)
    pts=10;
	lx=length(p(k).XData);
    p(k).MarkerIndices =1+k*floor(lx/pts/length(p)):floor(lx/pts):lx; % pts : nombre de markers 
    p(k).LineStyle='-';			
    p(k).Marker=Marker(k,1);				
    p(k).MarkerSize=5;
end