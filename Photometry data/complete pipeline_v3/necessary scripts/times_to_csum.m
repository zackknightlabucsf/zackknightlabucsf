%Converts timestamps to cumulative sums
%units are the amount by which x increases to endtime
%i.e. how often you want cumulative sum measured
%timestamps before begintime are removed
function [x,y]=times_to_csum(timestamps, endtime,units,begintime)
x=(0:endtime/units)*units;
y=[];
timestamps(find(timestamps<begintime))=[];
for i=x;
    y=[y length(find(timestamps<i))];
end
end