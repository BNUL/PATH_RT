script_full_path = mfilename('fullpath');
[script_path, ~, ~] = fileparts(script_full_path);
PATH_root = [script_path '\PLD\'];

%% Sun and observation geometry
SZA = 41;  
SAA = 147;
VAAband = [327,327,327,327,327,327,327,327,327,327,327,327,327,327,147,147,147,147,147,147,147,147,147,147,147,147,147,147,147,147,147,147];
VZAband = [70,65,60,55,50,45,40,35,30,25,20,15,10,5,0,5,10,15,20,25,30,35,39,40,41,43,45,50,55,60,65,70];


%% Canopy information
iorien = 6;            %leaf orientation type£º1-planophile; 2-erectophile; 3-plagiophile;4-extremophile;5-uniform;6-spherical
FAVD = 0.8;            %foliage area volume density
Height = 19.2;         %Canopy height
HotSpotPar = 0.1;      %hotspot parameter at leaf scale (leaf size related)
tau_L = [0.08998 0.02048 0.43134]; %leaf transmittance; each column corresponds to one band
ro_L = [0.1109 0.05928 0.45082];  %leaf reflectance
ro_s = [ones(size(VZAband,2),1).*0.045298 ...
    ones(size(VZAband,2),1).*0.027397 ...
    ones(size(VZAband,2),1).*0.191457]; %soil reflectance
step1 = 1; %BRF_S;   '1' means the process is included
step2 = 1; %BRF_L
step3 = 1; %BRF_CM
step4 = 1; %BRF_GCM
% BRF = BRF_S + BRF_L + BRF_CM + BRF_GCM
Gassume = 0; % G can be field measured or assumed based the parameter iorien
if Gassume ~=0
G= [0.489890417605104	0.492020429278576	0.494434880542201	0.497108066956753	0.500248212410529	0.504034853843981	0.508435439666468	0.513468913953419	0.518941090533616	0.524110919285445	0.528360598840716	0.531440438474483	0.533221318660591	0.533672797840561	0.532333850979834	0.533672797840561	0.533221318660591	0.531440438474483	0.528360598840716	0.524110919285445	0.518941090533616	0.513468913953419	0.509484097788791	0.508435439666468	0.507529119549790	0.505822903506383	0.504034853843981	0.500248212410529	0.497108066956753	0.494434880542201	0.492020429278576	0.489890417605104];
end
%% end of Canopy information


PATH_S = load([PATH_root, sprintf('AA%d_ZA%d_layers.mat', SAA, SZA)]);
PATH_S = PATH_S.PATH;
PATH_S = permute(PATH_S, [2,3,1]);
%path distribution at viewing direction
PATH_V = ones(size(PATH_S)); 
PATH_H = load([PATH_root, 'AA',num2str(SAA),'_ZA0_layers.mat']);
PATH_H = PATH_H.PATH;
PATH_H = permute(PATH_H, [2,3,1]);
mtm = PATH_H(1,:)*PATH_H(2,:)';
Crowndeepth = mtm./(1-PATH_H(2,1));
Znum = floor(Crowndeepth./HotSpotPar.*FAVD);   %number of layers
dthr = get_dthr(PATH_S,PATH_H,SZA);
go_par = dthr*Crowndeepth; %hotspot parameter at crown scale (crown size related)
LAI = getLAI( FAVD, PATH_H )
w = tau_L + ro_L;   %leaf albedo
iD = getHemiInterceptancev4( iorien, FAVD, PATH_root,SAA);
PATH_RT
BRF_band = BRF_band(end:-1:1,:);
clear script_path
