function [f,Phasor,Cn,theta] = f_FourierTrapeze(A,T0,a,tm,td,Fend,DCb,varagin)
%UNTITLED Summary of this function goes here
%  ATTENTION: ne marche pas bien pour les signaux carrés. 
if nargin==8
    delay=varagin;
else
    delay=0;
end

Bm=tm/T0;
Bd=td/T0;
F0=1/T0;
n=1:Fend/F0;
n=n(:);
f=[0; F0*n];
% sinc(t)=sin(pi*t)/(pi*t) pour matlab
Cn=zeros(size(f));
theta=zeros(size(f));
Cn(1)=A*a+DCb;
theta(1)=0;
Cn(2:end)=(A./(pi*n)).*sqrt((sinc(n*Bm)).^2+(sinc(n*Bd)).^2-2*sinc(n*Bm).*sinc(n*Bd).*cos(pi*n*(2*a)));
phDelay=-2*pi*f(2:end)*delay;
if (Bm~=0 && Bd~=0)
    theta(2:end)= -atan2((Bd*sin(2*n*pi*Bm)+Bm*(sin(n*pi*(2*a-Bd+Bm))-sin((n*pi*(2*a+Bd+Bm))))),(-Bd+Bd*cos(2*n*pi*Bm)+Bm*cos(n*pi*(2*a-Bd+Bm))-Bm*cos(n*pi*(2*a+Bd+Bm))))+phDelay;
else
    if (Bm==0 && Bm==0)
        theta=[0 -a*n*pi]+phDelay;
    else
        warning('La phase n''a pas pu être calculée')
        theta=zeros(size(n));
    end
end
Cn=Cn(:);
theta=theta(:);
Phasor=Cn.*exp(1i*theta);
% verif
end

