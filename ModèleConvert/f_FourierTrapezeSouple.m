function [f,Phasor,Cn,theta] = f_FourierTrapezeSouple(A,T0,a,tm,td,Fend)
%UNTITLED Summary of this function goes here
%  ATTENTION: ne marche pas bien pour les signaux carrés. 
Bm=tm/T0;
Bd=td/T0;
F0=1/T0;
n=1:Fend/F0;
f=[0 F0*n];
f=f(:);
tau=min([td tm])/5;

% sinc(t)=sin(pi*t)/(pi*t) pour matlab
Cn=zeros(size(f));
Cn(1)=A*(a-Bm/2+Bd/2);
Cn(2:end)=(A./(pi*n)).*sqrt((sinc(n*Bm)).^2+(sinc(n*Bd)).^2-2*sinc(n*Bm).*sinc(n*Bd).*cos(pi*n*(2*a+Bd-Bm))).*sinc(tau*n);
if (Bm~=0 && Bd~=0)
    theta=[0 -atan2((Bd*sin(2*pi*n*Bm)+Bm*(sin(n*pi*a*2)-sin(2*n*pi*(a+Bd)))),(-2*Bd*sin(pi*n*Bm).^2+Bm*(cos(n*pi*a*2)-cos(2*n*pi*(a+Bd)))))];
else
    if (Bm==0 && Bm==0)
        theta=[0 -a*n*pi];
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

