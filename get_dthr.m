function dthr = get_dthr( PATH_S,PATH_H,SZA)
%This routine estimates integral G(Omega)
%(1/2*pi)int(0,pi/2)int(0,2*pi)g_L |(Omega,Omega_L)|d theta_L d fi
%Here g_L(Omega_L)=(2/pi)(1+a*cos(b*theta_L)

%type of leaf normal orientation(iorien):
%1-planophile (a=1,b=2); 2-erectophile (a=-1,b=2);
%3-plagiophile (a=-1, b=4); 4-extremophile(a=1,b=4);
%5-uniform (a=0, b any);6-spherical(g=sin(theta))
%theta,fi indicate the direction Omega in Degree !

%%
dthr = 1./sqrt(((log(PATH_S(2,1))./(log(PATH_H(2,1)))).^2-1)/(tand(SZA).^2));
end

