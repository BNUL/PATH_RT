function gamma = get_Gamma( SZA,SAA,VZA,VAA,ro_L,tau_L,iorien )
%This routine calculates area phase function (1/pi)Gamma(Omega'->Omega)
%extracted from SRT codes 'gind_v07.f'
%INPUT:
%SZA\SAA:incident direction (Degree !)
%VZA\VAA:scattering direction (Degree !)
%ro_L\tau_L:reflectance and transmittance of leaf
%iorien:type of leaf normal orientation
%1-planophile (a=1,b=2); 2-erectophile (a=-1,b=2);
%3-plagiophile (a=-1, b=4); 4-extremophile(a=1,b=4);
%5-uniform (a=0, b any);6-spherical

if iorien==6
    SZA=deg2rad(SZA+180);SAA=deg2rad(SAA);
    VZA=deg2rad(VZA);VAA=deg2rad(VAA);
    wleaf=ro_L+tau_L;
    cb=cos(SZA)*cos(VZA)+sin(SZA)*sin(VZA)*cos(VAA-SAA);
    b=acos(cb);
    gamma=wleaf.*(sin(b)-b*cb)/(3*pi)+tau_L.*cb/3;
    gamma=gamma./pi;
else
    a=1;b=4;
    if(iorien<=2)
       b=2; 
    end
    if(iorien==2 || iorien==3)
       a=-1; 
    end
    if(iorien==5)
       a=0; 
    end
    gamma=get_APF(a,b,ro_L,tau_L,SZA,SAA,VZA,VAA);
end

end

