function P_bidir = getBiDirGap(Pv_dir, Ps_dir, Chs)
% calculate canopy bi-directional gap at height z
% P(z, omega_s, omega_v) = P(z, omega_s) * P(z, omega_v) * Chs
% input:
% P(z, omega_s): gap at solar direction
% P(z, omega_v): gap at viewing direction
% Chs: hotspot

%%
P_bidir = Pv_dir * Ps_dir * Chs;

