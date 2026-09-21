addpath('D:\ownclound\Thèse\Rassemblement des modèles\Frequence2time')
%%
A=100;
T0=5e-5;
a=0.5;
tm=1e-7;
td=1e-7;
DCb=-0.5;
Fend=1/tm;
f=linspace(0,1/tm,1e5);

[ft,Trap,Cn,Phase] = f_FourierTrapeze(A,T0,a,tm,td,Fend,DCb);
[TrapEnv,CnEnv,PhaseEnv] = f_FourierTrapeze_Env(A,T0,a,tm,td,f,DCb);

t=linspace(0,T0,1e4);
[Trapt2] = f_ifseries(ft,Trap,t);

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
plot(t,Trapt,t,Trapt2)
xlabel('Time (s)');ylabel('Signal');
legend('Input','Reconstructed');

figure(12)
subplot 211
semilogy(ft,Cn,'.',f,CnEnv)
xlabel('Frequency (Hz');ylabel('Magnitude')
legend('Harmonics','Envelope','Location','best')
subplot 212
plot(ft,Phase*180/pi,'.',f,PhaseEnv*180/pi)
xlabel('Frequency (Hz');ylabel('Angle (deg)')
legend('Harmonics','Envelope','Location','best')