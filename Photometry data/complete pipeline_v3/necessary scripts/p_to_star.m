function [star]=p_to_star(p)
if(p>0.05)
    star='ns';
elseif(p>0.01)
    star='*';
elseif(p>0.001)
    star='**';
elseif(p<0.001)
    star='***';
end
end