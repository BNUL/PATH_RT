function LAI = getLAI( FAVD, PATH_H )
% calculate LAI
% input: PATH_H path distribution at the bottom layer in hemispherical space

%%
LAI = 0;
PATH = PATH_H(:, :, 1);  % path distribution at ZA = 0
for i = 1: size(PATH, 2)
    if PATH(2,i)<= 0.001
        LAI = LAI + 0;
    else
        LAI = LAI + FAVD * PATH(1, i) * PATH(2, i);
    end
end

end

