function get_APF = get_APF(a,b,rol,taul,cteta,phi,cteta1,phi1)
cteta1 = cteta1 + 180;
ccteta=cosd(cteta);
phi=deg2rad(phi);
ccteta1=cosd(cteta1);
phi1=deg2rad(phi1);

kpi3 = 1/(pi^3);
n = 30;
m = 4*n;
h_theta = 0.5*pi/n;
h_fi=2.*pi/m;

% steta=sqrt(abs(1.-ccteta*ccteta));
% steta1=sqrt(abs(1.-ccteta1*ccteta1));
steta = sind(cteta);
steta1 = sind(cteta1);
theta_i=0.5*h_theta;
fi_1=0.5*h_fi;
get_APF=0.;

for i = 1 : n
    fi_j=fi_1;
    xx=0.;
    c_i=cos(theta_i);
    s_i=sin(theta_i);
    for j = 1 : m
        yy=ccteta*c_i+steta*s_i*cos(phi-fi_j);
        zz=ccteta1*c_i+steta1*s_i*cos(phi1-fi_j);
        zz=zz*yy;
        if(zz <= 0.)
            xx=rol.*abs(zz);
        else
            xx=taul.*zz;
        end
        fi_j=fi_j+h_fi;    
    xx=xx*h_fi ;
    yy=1.+a*cos(b*theta_i) ;
     get_APF=get_APF+yy*xx;
       
    end
    theta_i=theta_i+h_theta;  
end
get_APF=h_theta*get_APF;
get_APF=kpi3*get_APF;
end