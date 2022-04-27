%make custom colormap

%zmax is max z-score the colormap goes up to
%coloron is the z-score that the colormap saturates at (usually set to zmax)
%coloroff is the z-score below which the colormap is white

function [cmap]=JG_cmap_v1_redder(zmax,coloron,coloroff)
zmin=-zmax;
z = [-200 -coloron -mean([coloron coloroff]) -coloroff 0 coloroff mean([coloron coloroff]) coloron 200];
r = [0 0 .4 1 1 1 .9 1 1];
g = [0.2 0.2 .7 1 1 1 .5 0 0];
b = [1 1 .8 1 1 1 .4 0 0];
cmap=[interp1(z, r, interp1([1 256],[zmin zmax],1:256, 'pchip'),'pchip')',interp1(z, g, interp1([1 256],[zmin zmax],1:256, 'pchip'),'pchip')',interp1(z, b, interp1([1 256],[zmin zmax],1:256, 'pchip'),'pchip')'];
end