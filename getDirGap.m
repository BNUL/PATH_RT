function Pdir = getDirGap(G, FAVD, PATHz)
% calculate canopy directional gap at height z
% P(z, omega) = sum( exp(-G * FAVD * l) * p(l) )
% input:
% G: mean projection of leaf normal at direction omega
% FAVD: foliage area volume density
% PATH: directional path distribution at layer z

%%
Pdir = 0;
for i = 1: size(PATHz, 2)
    Pdir = Pdir + exp(-G * FAVD * PATHz(1, i)) * PATHz(2, i);
end

end
