function [tout,sigout] = f_PasConstant(t,sig)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
dt=min(diff(t));

tout=t(1):dt:t(end);
tout=tout(:);
sigout=zeros(size(tout));
sigout(1)=sig(1);
for j=1:length(t)-1
   ind=logical(logical(tout<=t(j+1)).*logical(tout>t(j)));
   sigout(ind)=sig(j)+diff(sig(j:j+1))/diff(t(j:j+1))*(tout(ind)-t(j));
end

end