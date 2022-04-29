for i = 1:4 
%i=1
    this = psth_all{i};
    first = this{1};
    second = this{13};
    third = this{25};
    for j = 1:size(first,1)
        temp = [first(j,:); second(j,:); third(j,:)];
        psth_16_mean(j,:) = mean(temp,1);
        a=mean(temp,1);
        averages(j) = mean(a(30:end));
    end
    mouse_Suc{i} = psth_16_mean;
    means_Suc{i} = averages; 

end

paired = [];
this = means{1}; %TA120
this_Suc = means_Suc{1};
TA120_Suclose = [1 2 5 7 8 9 11 18];
TA120_Suc = [1 2 3 4 5 6 8 14];

for i = 1:length(TA120_Suc)
    paired = [paired; this(TA120_Suclose(i)) this_Suc(TA120_Suc(i))];
end

this = means{2}; %TA121
this_Suc = means_Suc{2};
Suclose = [1 3 4 5 8 9 10 13 14 15 16 17 18 19 20 21 22];
Suc = [1 2 22 17 8 5 9 11 10 15 23 16 18 21 19 20 4];

for i = 1:length(Suc)
    paired = [paired; this(Suclose(i)) this_Suc(Suc(i))];
end

this = means{3}; %TA124
this_Suc = means_Suc{3};
Suclose = [1 2 3 4 5 6 8 10 12 13 14 15 16 17 18 19 20 21 22 23 25 26 27 28 29];
Suc = [8 4 1 7 12 14 2 19 39 21 38 20 13 9 37 16 6 24 29 30 34 32 36 33 10];

for i = 1:length(Suc)
    paired = [paired; this(Suclose(i)) this_Suc(Suc(i))];
end

this = means{4}; %TA130
this_Suc = means_Suc{4};
Suclose = [1 2 3 4 6 7 9 10 12 13 14 15];
Suc = [1 4 3 5 17 10 14 8 22 21 18 16];

for i = 1:length(Suc)
    paired = [paired; this(Suclose(i)) this_Suc(Suc(i))];
end

%% 
AA = 0; AN=0;AI=0;NA=0;NN=0;NI=0;IA=0;IN=0;II=0;
for i = 1:size(paired,1)
    if paired(i,1) >=1 
        if paired(i,2) >=1
            AA = AA+1;
        elseif paired(i,2) <=-1
            AI = AI+1;
        else
            AN = AN+1;
        end
    elseif paired(i,1) <1 && paired(i,1) >-1
        if paired(i,2) >=1
            NA = NA+1;
        elseif paired(i,2) <=-1
            NI = NI+1;
        else
            NN = NN+1;
        end
    elseif paired(i,1) <=-1
        if paired(i,2) >=1
            IA = IA+1;
        elseif paired(i,2) <=-1
            II = II+1;
        else
            IN = IN+1;
        end
    end
end