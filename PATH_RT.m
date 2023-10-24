%**************************************************************************
% Program Name: PATH model
% 
% Description: Modeling the canopy BRF based on path length distribution.
% 
% Author: Weihua Li, Beijing Normal University, Beijing, China.
% Email: weihuali@mail.bnu.edu.cn
% Phone: +(86) 130-3288-6110
%**************************************************************************

%% MATLAB code starts here...
BRF_band = ones(size(VZAband,2),size(w,2));
pathAngle =0;
tic
for m = 1:length(VZAband)   %1:length(VZAband)    
%% Step 1: calculate BRF_sgl
    VZA = VZAband(m);
    VAA = VAAband(m);    
    PATH_V = load([PATH_root, sprintf('AA%d_ZA%d_layers.mat', VAA, VZA)]);
    PATH_V = PATH_V.PATH;
    PATH_V = permute(PATH_V, [2,3,1]);
    if Gassume == 1
        Gs = get_G(iorien, SZA, SAA);
        Gv = get_G(iorien, VZA, VAA);
    else
        HSindex = find(VZAband == SZA);
        Gs = G(HSindex);
        Gv = G(m);
    end
    Ps_dir_inKz =  getDirGap(Gs, FAVD, PATH_S(:,2:end))./(1-PATH_S(2,1));
%     Crowndeepth_from_gap = getlequivalent(Gv, FAVD,30, Ph_dir_inKz);
    if pathAngle == 1
        ls = Crowndeepth/cosd(SZA) * (1-PATH_H(2,1));
        lv = Crowndeepth/cosd(VZA) * (1-PATH_H(2,1));
        %ls = -Crowndeepth/cosd(SZA) * log(PATH_H(2,1));
        %lv = -Crowndeepth/cosd(VZA) * log(PATH_H(2,1));        
    else
        % ls = getlequivalent(Gs, FAVD,30, Ps_dir_inKz);
        % lv = getlequivalent(Gv, FAVD,30, Pv_dir_inKz);        
        lv = PATH_V(1, :)*PATH_V(2, :)'/(1-PATH_V(2,1)); %path length;
        ls = PATH_S(1, :)*PATH_S(2, :)'/(1-PATH_S(2,1)); %path length;        
    end
    ls3 = Crowndeepth/ls; % cos
    lv3 = Crowndeepth/lv;
%% Step 1.1: calculate BRF_S
    BRFs = 0;
    BRFs_kt = 0;
    if step1 == 1
        Ps_dir =  getDirGapHomo(Gs, FAVD, Crowndeepth, ls3);
        Pv_dir =  getDirGapHomo(Gv, FAVD, Crowndeepth, lv3);
        Chs = get_HSF(HotSpotPar,SZA,SAA,VZA,VAA,Gs,Gv,FAVD,Crowndeepth,ls3,lv3);
        BRFs = getBiDirGap(Pv_dir, Ps_dir, Chs) .* ro_s(m,:);
        BRFs_kt = getBiDirGap(Pv_dir, Ps_dir, 1) .* ro_s(m,:);
    end
%% Step 1.2: calculate BRF_L
    BRF_v1 = 0;
    BRF_v1_kt = 0;
    if step2 == 1
        gamma = get_Gamma( SZA,SAA,VZA,VAA,ro_L,tau_L,iorien );
        P_bidir = ones(1, Znum);  % bi-directional gap at each layer
        P_bidir_kt = ones(1, Znum);
        for z = 1: Znum
            h = Crowndeepth / Znum * z;
            Ps_z = getDirGapHomo(Gs, FAVD, h, ls3);
            Pv_z = getDirGapHomo(Gv, FAVD, h, lv3);
            Chs = get_HSF(HotSpotPar,SZA,SAA,VZA,VAA,Gs,Gv,FAVD,h,ls3,lv3);
            P_bidir(1, z) = getBiDirGap(Pv_z, Ps_z, Chs);
            P_bidir_kt(1, z) = getBiDirGap(Pv_z, Ps_z, 1);
        end
        z = Crowndeepth / Znum : Crowndeepth / Znum: Crowndeepth;
        BRF_v1 = gamma .* FAVD ./ lv3 .* trapz(z, P_bidir).*pi./cosd(SZA);
        BRF_v1_kt = sqrt(Ps_dir_inKz).*gamma .* FAVD ./ lv3 .* trapz(z, P_bidir_kt).*pi./cosd(SZA);
    end
%% Step 2: calculate four components
    if pathAngle == 1
        Ps_dir_go =  exp(log(PATH_H(2,1))./cosd(SZA));
        Pv_dir_go =  exp(log(PATH_H(2,1))./cosd(VZA));
    else
        Ps_dir_go = PATH_S(2,1);
        Pv_dir_go = PATH_V(2,1);
    end
    
    if Height - PATH_H(1,end) < 0
        lmax = Crowndeepth;
    end
    Kg = Ps_dir_go*Pv_dir_go + get_HSF_go(go_par,SZA,SAA,VZA,VAA,Ps_dir_go,Pv_dir_go,Height - lmax + 1/2*Crowndeepth);
    Kz = Pv_dir_go - Kg;
    Kct = 1-Pv_dir_go;
    delta = cosd(SZA)*cosd(VZA)+sind(SZA)*sind(VZA)*cosd(VAA-SAA); %
    phi = acosd(delta);
    c = mtm./max(PATH_H(1,:));
    delta = cosd(phi.*(1 - sin(pi.*c./2))); 
    Kc = 0.5*(1+delta)*Kct;
    Kt = Kct - Kc;
%% Step 3: calculate BRF_CM
    BRF_vm = 0;
    if(step3)        
        i0 = 1- getDirGap(Gs, FAVD, PATH_S(:,:));
        iv = 1- getDirGap(Gv, FAVD, PATH_V(:,:));
        P_recol = 1 - iD / LAI;    % recollision probability
        tao_v = iv / (2 * LAI);    % escape probability at viewing direction
        tao_hemi = iD / (2 * LAI); % escape probability in hemispherical space
        BRF_vm = i0 .* w.^2 .* P_recol .* tao_v ./ (1 - P_recol .* w);
    end   
%% Step 4: calculate BRF_GCM
    BRF_vs = 0;
    if(step4)        
        Tdn = 1 - i0 + i0 .* w .* tao_hemi ./ (1 - P_recol .* w);
        Tup = 1 - iv + iD .* w .* tao_v ./ (1 - P_recol .* w);   % t_dif replaced by Pv_dir, required to be checked !
        Rdn = iD .* w .* tao_hemi ./ (1 - P_recol .* w);
        BRF_vs = ro_s(m,:) ./ (1 - ro_s(m,:) .* Rdn) .* Tdn .* Tup -(1 - i0) .* ro_s(m,:) .* (1 - iv);        
    end       
%% Final: calculate BRF    
    BRF_band(m,:) = Kc.*(BRFs + BRF_v1)...
        + Kt.*(BRFs_kt + BRF_v1_kt)...
        +(Kg+ Kz.*Ps_dir_inKz).*ro_s(m,:)...
        + BRF_vm...
        + BRF_vs;    
end

toc

%% clear
clear step1 step2 step3 step4 pathAngle Znum z PATH_root Height h Crowndeepth_from_gap
clear VAA VZA m lv ls ls3 lv3 HSindex delta phi SAA SZA 
clear P_bidir P_bidir_kt Ph_dir_inKz 
clear Ps_dir Ps_z Ps_dir_go Ps_dir_inKz
clear Pv_dir Pv_z Pv_dir_go Pv_dir_inKz  
clear go_par HotSpotPar Chs Crowndeepth FAVD G
clear iv i0 iD P_recol tao_v tao_hemi tau_L Tdn Tup Rdn  ro_L ro_s w  
clear Gs Gv iorien Gassume gamma mtm   
clear Kc Kct Kg Kt Kz
clear BRF BRF_v1 BRFs_kt BRF_v1_kt BRF_vs BRF_vm BRFs
