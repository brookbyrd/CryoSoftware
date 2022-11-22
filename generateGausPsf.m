function [h] = generateGausPsf(config, n, s1, s2)

clear h;
% Make xx yy in um space to use sigma in um space
[xx,yy] = meshgrid(linspace(-n/2*42,n/2*42,n));
rr2 = xx.^2 + yy.^2;

sigma = s1.*config.sliceThickness + s2;
h = exp(-rr2/(2.*sigma^2)); % 

end
%%
