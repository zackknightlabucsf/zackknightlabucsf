function [nx,ny,ndy]=JG_downsample_v3(y,filterby,ds) %ds is a number
y=y';
windowSize=filterby;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
ny=[];
y2=flipud(y);
for(i=1:size(y,2))
    if(ds>0)
        tempy=filter(b, a, y(:,i));
        tempy2=filter(b, a, y2(:,i));
        tempy2=flipud(tempy2);
        tempy(1:filterby)=tempy2(1:filterby);
        ny=[ny;downsample(tempy',ds)];
    else
        tempy=filter(b, a, y(:,i));
        tempy2=filter(b, a, y2(:,i));
        tempy2=flipud(tempy2);
        tempy(1:filterby)=tempy2(1:filterby);
        ny=[ny;tempy'];
    end
end
ndy=std(ny,1)/sqrt(size(ny,1));ndy=ndy';
ny=mean(ny,1);ny=ny';
nx=1:size(ny,1);nx=nx';
end