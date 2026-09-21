function [s] = f_ifseries(f,X,t)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
r=size(X);
if r(1)<r(2)
    X=transpose(X);
end
t=t(:);
s=zeros(length(t),size(X,2));
for i=1:length(f)
    s=s+(ones(size(t))*abs(X(i,:))).*cos(2*pi*f(i)*t+ones(size(t))*angle(X(i,:)));
end
end

