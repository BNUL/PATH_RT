function M = getlequivalent(Gs, FAVD,H, gap)
%GETLEQUVALENT �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
m = [0.01:0.01:H];
gaps = exp(-Gs*FAVD*m);
dis = abs(gaps-gap);
index = find(dis == min(dis));

index = index(1);

M = m(index);
end

