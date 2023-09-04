function Pdir = getDirGapHomo(G, FAVD, h , m)
% calculate canopy directional gap at height z
% P(z, omega) = sum( exp(-G * FAVD * l) * p(l) )
% input:
% G: mean projection of leaf normal at direction omega
% FAVD: foliage area volume density
% PATH: directional path distribution at layer z

%%


Pdir = exp(-G * FAVD * h /m); 


end
