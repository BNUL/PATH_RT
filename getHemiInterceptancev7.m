function iD = getHemiInterceptancev7( iorien, LAI,aa)
% calculate hemeispherical interceptance
% i_D = 1/pi * int(0,2pi)int(0, pi/2){ (1-P(H, Omega)) * cos¦Èsin¦È }d¦Èdfi
% input:
% PATH_H: path distribution at the bottom layer in hemispherical space
% ZA_delta, AA_delta: angular resolution

%%

za = 0:3:90;

% za = deg2rad(za);
iv = ones(size(za));


for i = 1:size(za,2)
        G = get_G(iorien, za(i), aa);

        p = exp(-G*LAI./cosd(za(i)));
        iv(i) = (1-p) * sind(2.*za(i));
   
end

iD = trapz(deg2rad(za), iv);    % required to be checked!!!

end

