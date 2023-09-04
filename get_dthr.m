function dthr = get_dthr( PATH_S,PATH_H,SZA)
dthr = 1./sqrt(((log(PATH_S(2,1))./(log(PATH_H(2,1)))).^2-1)/(tand(SZA).^2));
end

