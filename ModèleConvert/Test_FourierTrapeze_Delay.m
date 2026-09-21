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
t=linspace(0,T0,1000);

[ft,Trap] = f_FourierTrapeze(A,T0,a,tm,td,Fend,DCb);
[Trapt] = f_ifseries(ft,Trap,t);
delay=T0/5;
[ftd,Trapd] = f_FourierTrapeze(A,T0,a,tm,td,Fend,DCb,delay);
[Trapdt] = f_ifseries(ftd,Trapd,t);

figure(1)
plot(t,Trapt,t,Trapdt)