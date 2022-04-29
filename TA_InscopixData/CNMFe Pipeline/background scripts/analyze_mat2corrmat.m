function [ corrmat ] = analyze_mat2corrmat( matrix )
%UNTITLED2 get matrix and then do cross correlation across the matrix rows
%   Detailed explanation goes here
siz=size(matrix);
thiscell=mat2cell(matrix,[1:siz(1)]*0+1,siz(2));

c_co={};
for i=1:siz(1);
    c_co{i}=cellfun(@(x) corr2(x,thiscell{i}),thiscell,'UniformOutput',false);
end

c_co=cellfun(@(x) cell2mat(x),c_co,'UniformOutput',false);
corrmat=cell2mat(c_co);

end

