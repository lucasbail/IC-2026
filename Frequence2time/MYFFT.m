function [f,fft_x,ang_x,ts,xs]=MYFFT(t,x,Te,N,TT,interp,detrend)
%MYFFT     Version de FFT avec reechantillonnage et selection de fenetre de temps
%
% SYNTAXE:
%   [f,fft_x,ang_x]=myfft(t,x,Te,N,TiDur,Interp,detrend)
% DESCRIPTION:
%   Calcule la FFT avec conservation des amplitudes.
%   Il affiche aussi la courbe temporelle de la partie traitee.
% ENTREES:
%   t,x    : Donnees a traiter.
%   Te     : (optional) Periode d'echantillonnage.
%   N      : (optional) Nombre de points avec lesquels on realise la FFT.
%            Par defaut - N=length(x).
%   TiDur  : (optional) Vecteur donnant le temps initial et la duree 
%            a partir desquels on realise la FFT.
%            TT(1): Temps initial.
%            TT(2): Duree. Si l'on donne cette valeur, N est ignore.
%            Par defaut - TT=[min(t)]
%   Interp : (optional) Metode d'interpolation.
%            0 : Lineaire (valeur par defaut)
%            1 : Bloqueur d'Ordre Zero (pour les signaux discretes)
%   Detrend: (optional) Enlever la partie moyenne ou une derive lineaire.
%            0 : Utiliser le signal tel quel (valeur par defaut)
%            1 : Enlever la valeur moyenne
%            2 : Enlever une derive lineaire
% SORTIES:
%   f    : Frequence 
%   fft_x: amplitude de la FFT de la partie selectionnee de x
%   ang_x: angle de la FFT de la partie selectionnee de x
%
% NOTES:
%   Le temps doit etre monotone croissant.
%   Cette fonction essaie de detecter la periodicité du signal
%   d'entree et prends les valeurs a une periode constante
%   pour que ca marche. Cependant s'il existe beacoup de points
%   entre deux periodes d'echantillonnage (MLI par exemple)
%   Il ne peut plus detecter la periode. Dans ce cas il
%   realise un surechantillonnage avec une interpolation.
%
%   Si le nombre de points n'est pas specifie, il enleve toujours 
%   le dernier point de la sequence donnée. Si l'on specifie la durée
%   du signal egalement il enleve un point.
%
% AUTEUR: Eñaut Muxika Olasagasti
% DATE  : Le 29 Janvier 1999.

% Traitement des valeurs par defaut
lent=length(t);
if nargin<2,
   error('Il faut donner au moins le signal d''entree');
end;
if min(size(t))~=1 | min(size(x))~=1,
   error('Il faut fournir deux vecteurs [t,x]');
end;
x=x(:);t=t(:);gette=0;
if nargin<3,
   gette=1;
elseif isempty(Te),
   gette=1;
end;
if nargin<4,
   N=[];
end;
if nargin<5,
   TT=min(t);
end;
zoh=0;
if nargin>=6 & ~isempty(interp),
   zoh=all(all(interp == 1));
end;
if isempty(TT),
   ti=min(t);
else
   ti=TT(1);
end;
if length(TT)>1,
   flag=1;
   tf=ti+TT(2);
else
   flag=0;
   tf=max(t);
end;
if nargin<7,detrend=0;end;
quiet=0;
global MYFFT_SILENT
if ~isempty(MYFFT_SILENT) & MYFFT_SILENT == 1,quiet=1;end;

% Verification des erreurs
if(ti>tf),disp('Duree negative (ti>tf)!!!');end;
if(ti<min(t) & ~quiet),disp('Temps initiale du fft < Temps initiale du signal');end;

% Temps final specifiee
if flag==1,
   % Verifier s'ils sont correctes
   if(tf>max(t) & ~quiet),disp('Temps finale du fft < Temps finale du signal');end;
   % Il ne prend pas le point t==tf
   idx=find(t>=ti & t<tf);
else
   idx=find(t>=ti);
end;
% Prendre les valeurs qui seront traites.
t=t(idx);x=x(idx);

% obtenir la periode d'echantillonnage
dt=diff(t); 	% ne pas utiliser les points au debut et a la fin
if any(dt==0),warning('Le temps n''est pas monotone croissant!!');end;
Npt=[];
if gette,
   Tex=max(dt);Tem=min(dt);
   DT=max(t)-min(t);
   Te=Tex;
   % select the nearest bigger power of two number of points
   % when calculation Te for interpolation
   if abs(Tem-Tex)>Tex*1e-7,
      Npt=4*2^ceil(log(DT/Tex)/log(2));
      if Npt<length(t),
         Npt=2^ceil(log(length(t))/log(2));
      end;
      Te=DT/Npt;
   end;
end;
dt=[dt;inf];  	% ne pas utiliser les points au debut et a la fin
if ~quiet,disp(['Periode d''echantillonnage :',num2str(Te)]);end;

% resampling & interpolation to a constant data rate
ts=[ceil(min(t)/Te-eps):floor(max(t)/Te+eps)]'*Te;
xs=zeros(size(ts));

% Recherche des indexes utilises pour l'interpolation
idx=zeros(size(ts));
in=1;
for ins=1:length(ts),
   while ts(ins) >t(in),in=in+1;end;
   if    ts(ins)==t(in),idx(ins)=in;
   else                 idx(ins)=in-1;end;
end;
if ~quiet,
	if min(diff(idx))==0,
	   disp('Interpolation requise');
	else
	   disp('Interpolation non requise');
   end;
end;
if zoh>0,
   % Interpolation discrete (Bloqueur d'Ordre Zero)
   xs=x(idx);
else
   % Interpolation lineaire (Bloqueur d'Ordre Un)
   dx=[diff(x);0];dxs=dx(idx);
   dts=dt(idx);res=(ts-t(idx))./dts;
   ix=find(~isfinite(res));
   if ~isempty(ix) & ~quiet,
      disp('MYFFT: residus == NaN or Inf !!');res(ix)=1;
   end;
   % Realiser l'interpolation
   xs=x(idx) + dxs.*res;
end;
clear in ins dt dx ix dts dxs res idx;

% Choisir le nombre de points utilises dans le calcul
if isempty(Npt),lent=length(ts);
else            lent=Npt;end;
if flag==0,
   % Detecter si ca depasse le point finale des donnees
   % et agir en consequence (limiter le nombre de points)
   if isempty(N),N=lent;
   elseif N>lent,
      N = lent;
      if ~quiet,
	      disp(['ATTENTION: Calcul avec ',num2str(N),...
               ' points']);
      end;
   end;
else
   N=lent;
end;

% Gerer la valeur moyenne et la derive
switch(detrend),
case 1, xs(1:N)=xs(1:N)-mean(xs(1:N));
case 2; xs(1:N)=detrend(xs(1:N));
end;

% Creer le vecteur frequence
idx=[0:(N-1)/2];
fr=idx/N; % frequence reduite
f = fr/Te;

% realiser la FFT et prendre les frequences positives.
if ~quiet,disp(['Debut : ',num2str(N),' points']);end;
tic;fft_x = fft(xs,N)/N;tee=toc;
if ~quiet,disp(['Fin : ',num2str(tee),'s']);end;
fft_x = fft_x(1+idx)*2;
% Corriger les amplitudes en tenant en compte les frequences negatives
fft_x(1)=fft_x(1)*0.5;
% Obtenir l'amplitude et l'angle
ang_x=angle(fft_x);
fft_x=abs(fft_x);

% Afficher le signal temporel et son FFT
if nargout<1,
   subplot(211);
   plot(ts,xs);grid;zoom on;
   xlabel('Temps');axisbox;
   subplot(212);
   plot(f,fft_x);grid;zoom on;
   xlabel('Freq.');axisbox;
end;

