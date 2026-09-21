addpath('ModèleConvert','ModèleLT','ModèleLT\IdentificationPont2','PSpice sim');
%% Source
% Changer les paramètres pour correspondre à Vd
F0=20e3;
T0=1/F0;
tm=60e-9;
td=45e-9;
a=0.4025;
V=120;
Fend=10/td;
[f,Vd] = f_FourierTrapeze(V,T0,a,tm,td,Fend);

figure(10)
h=gca;
loglog(f,abs(Vd),'.')
h.XLim=[f(1) f(end)];
xlabel('Frequency (Hz)');ylabel('Magnitude (V)');
title('Trapezoidal source Fourier decomposition')
%% Solution fréquentielle
ligne='separatexLargeBande';
Lc=12;
Zs=0;typeS='R';
Zl=[10 500e-6 2.2e-9];typeL='RLsCp';
[zmes] =[0,Lc];

[Vz,Iz]= fZ_VzIzCompletSansOscillo(Vd,Zs,typeS,zmes,f,Lc,Zl,typeL,ligne);

%% Passage au domaine du tem
t=linspace(0,1/F0,10000);

[Izt] = f_ifseries(f,Iz,t);
[Vzt] = f_ifseries(f,Vz,t);
Vdt=f_ifseries(f,Vd,t);

fr=1/min([td,tm]);
ind=logical(f<=fr);
[Iztfr] = f_ifseries(f(ind),Iz(ind,1),t);
[Vztfr] = f_ifseries(f(ind),Vz(ind,2),t);

figure(1);h=gca;
plot(t,Vzt(:,2),t,Vztfr,'--')
h.XLim=[t(1) t(end)]; 
title('Voltage at z=l');xlabel('Time (s)');ylabel('Voltage (V)')
legend('Up to 10/t_r','Up to 1/t_r')

figure(2);h=gca;
plot(t,Izt(:,1),t,Iztfr,'--')
h.XLim=[t(1) t(end)]; 
title('Current at z=0');xlabel('Time (s)');ylabel('Current (A)')
legend('Up to 10/t_r','Up to 1/t_r')

fr2=1/max([td,tm]);
ind2=logical(f<=fr2);
[Iztfr2] = f_ifseries(f(ind2),Iz(ind2,1),t);
[Vztfr2] = f_ifseries(f(ind2),Vz(ind2,2),t);

figure(3);h=gca;
plot(t,Vzt(:,2),t,Vztfr2,'--')
h.XLim=[t(1) t(end)]; 
title('Voltage at z=l');xlabel('Time (s)');ylabel('Voltage (V)')
legend('Up to 10/t_m','Up to 1/t_m')

figure(4);h=gca;
plot(t,Izt(:,1),t,Iztfr2,'--')
h.XLim=[t(1) t(end)]; 
title('Current at z=0');xlabel('Time (s)');ylabel('Current (A)')
legend('Up to 10/t_m','Up to 1/t_m')


fr3=0.35/max([td,tm]);
ind3=logical(f<=fr3);
[Iztfr3] = f_ifseries(f(ind3),Iz(ind3,1),t);
[Vztfr3] = f_ifseries(f(ind3),Vz(ind3,2),t);

figure(5);h=gca;
plot(t,Vztfr2,t,Vztfr3,'--')
h.XLim=[t(1) t(end)]; 
title('Voltage at z=l');xlabel('Time (s)');ylabel('Voltage (V)')
legend('Up to 1/t_r','Up to 0.35/t_r')

figure(6);h=gca;
plot(t,Iztfr2,t,Iztfr3,'--')
h.XLim=[t(1) t(end)]; 
title('Current at z=0');xlabel('Time (s)');ylabel('Current (A)')
legend('Up to 1/t_r','Up to 0.35/t_r')
