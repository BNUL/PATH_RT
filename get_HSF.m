function  Chs=get_HSF( par,SZA,SAA,VZA,VAA,Gs,Gv,uL,z,ls,lv) 
%This routine calculates the hot spot factor based on Kuusk
%input:
%par:1/SL 
%geometry in degree!
%Gs\Gv:G function at solar and viewing directions
%z:the height of a layer
%CPF:correlation function
%uL:FAVD
%z:height z
%output:
%Chs:hot spot factor at height z
%%
SZA=deg2rad(SZA);SAA=deg2rad(SAA);
VZA=deg2rad(VZA);VAA=deg2rad(VAA);
mu0=ls;
muv=lv;
f1=sqrt(Gs*Gv/(mu0*muv));
%calculate delta:
phi = VAA-SAA;
if(phi < 0 && phi > - 360)
   phi = phi + 360;
elseif(phi > 360 && phi < 720)
   phi = phi - 360;
else
   phi = phi;
end
cosgamma=cos(SZA)*cos(VZA)+sin(SZA)*sin(VZA)*cos(VAA-SAA);
delta=sqrt(1/(cos(SZA)^2)+1/(cos(VZA)^2)-2*cosgamma/(cos(SZA)*cos(VZA)));
if delta<0.00001
    delta=0.00001;
end
f2=uL* par/(delta)*(1-exp(-z*delta/par));
Chs = exp(f1 * f2);

end

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     