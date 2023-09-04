function  Chs=get_HSF_go( par,SZA,SAA,VZA,VAA,Ps_dir_go,Pv_dir_go,z)

SZA=deg2rad(SZA);SAA=deg2rad(SAA);
VZA=deg2rad(VZA);VAA=deg2rad(VAA);
mu0=cos(SZA);
muv=cos(VZA);
f1=sqrt(Ps_dir_go*Pv_dir_go*(1-Ps_dir_go)*(1-Pv_dir_go));
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
f2=exp(-delta/par*z);
Chs = f1 * f2;

end

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
