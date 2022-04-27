%make custom colormap
%function [nx,ny,ndy]=JG_downsample_v3(y,filterby,ds) %ds is a number
function [cmap]=JG_cmap_v1(zmax,coloron,coloroff)

zmin=-zmax;
z = [-200 -coloron -mean([coloron coloroff]) -coloroff 0 coloroff mean([coloron coloroff]) coloron 200];
r = [0 0 .4 1 1 1 .9 .4 .4];
g = [0.2 0.2 .7 1 1 1 .5 0 0];
b = [.4 .4 .8 1 1 1 .4 0.1 0.1];

cmap=[interp1(z, r, interp1([1 256],[zmin zmax],1:256, 'pchip'),'pchip')',interp1(z, g, interp1([1 256],[zmin zmax],1:256, 'pchip'),'pchip')',interp1(z, b, interp1([1 256],[zmin zmax],1:256, 'pchip'),'pchip')'];
end