function iD = getHemiInterceptancev4( iorien, FAVD, PATH_root,aa)
% calculate hemeispherical interceptance
% i_D = 1/pi * int(0,2pi)int(0, pi/2){ (1-P(H, Omega)) * cos¦Èsin¦È }d¦Èdfi
% input:
% PATH_H: path distribution at the bottom layer in hemispherical space
% ZA_delta, AA_delta: angular resolution

%%
za = [0 : 5 : 85 89];

% za = deg2rad(za);
iv = ones(size(za));

count = 0;
for i = 1:size(za,2)
        G = get_G(iorien, za(i), aa);
        
        PATH_V = load([PATH_root, sprintf('AA%d_ZA%d_layers.mat', aa, za(i))]);
        PATH_V = PATH_V.PATH;
        PATH_V = permute(PATH_V, [2,3,1]);
        p = getDirGap(G, FAVD, PATH_V(:,:));

        
       
        iv(i) = (1-p) * sind(2.*za(i));
   
end


iD = trapz(deg2rad(za), iv);    % required to be checked!!!

end

