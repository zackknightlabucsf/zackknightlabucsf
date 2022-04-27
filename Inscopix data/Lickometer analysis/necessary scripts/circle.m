function [c]=circle(x,y,r)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
c=[x+xp;y+yp]';
end
%square
% p=[0,0; 0,1; 0,2; 1,2; 2,2; 2,1; 2,0; 1,0; 0,0]
% plot(p(:,1),p(:,2))
