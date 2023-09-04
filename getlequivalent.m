function M = getlequivalent(Gs, FAVD,H, gap)
%GETLEQUVALENT 此处显示有关此函数的摘要
%   此处显示详细说明
m = [0.01:0.01:H];
gaps = exp(-Gs*FAVD*m);
dis = abs(gaps-gap);
index = find(dis == min(dis));

index = index(1);

M = m(index);
end

