%Converts timestamps to binary representation
%units are the amount by which x increases to endtime
%i.e. how often you want cumulative sum measured
%timestamps before begintime are removed
function [x,y]=times_to_bin(timestamps, endtime,units,begintime)
x=(0:endtime/units)*units;
y=0*x;
timestamps=reshape(timestamps,1,length(timestamps));
timestamps(find(timestamps<begintime))=[];
x=[x,timestamps];
y=[y, repelem(1, length(timestamps))];
[~,nind]=sort(x);
x=x(nind);
y=y(nind);
end