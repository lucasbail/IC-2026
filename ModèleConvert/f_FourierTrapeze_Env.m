function [Phasor,Cn,theta] = f_FourierTrapeze_Env(A,T0,a,tm,td,f,DCb)
%UNTITLED Summary of this function goes here
%  ATTENTION: ne marche pas bien pour les signaux carrés. 
% sinc(t)=sin(pi*t)/(pi*t) pour matlab
Cn=(A./(pi*f*T0)).*sqrt((sinc(f.*tm)).^2+(sinc(f*td)).^2-2*sinc(f*tm).*sinc(f*td).*cos(pi*f*T0*(2*a)));
if (tm~=0 && td~=0)
    theta= -atan2((td*sin(2*pi*f*tm)+tm*(sin(pi*f*(T0*2*a-td+tm))-sin(pi*f*(T0*2*a+td+tm)))),(-td+td*cos(2*f*pi*tm)+tm*(cos(pi*f*(T0*2*a-td+tm))-cos(pi*f*(T0*2*a+td+tm)))));
else
    if (tm==0 && tm==0)
        theta=[0 -a*f*T0*pi];
    else
        warning('La phase n''a pas pu être calculée')
        theta=zeros(size(f));
    end
end
Cn(1)=A*a+DCb;
theta(1)=0;
Cn=Cn(:);
theta=theta(:);
Phasor=Cn.*exp(1i*theta);
% verif
end

