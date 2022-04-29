%% Davis_LickRate
Davis_file = { %output file from Davis converted to csv
    'TA120_fastedSucrose2_1213.csv'
    'TA121_fastedSucrose2_1213.csv'
    'TA124_fastedSucrose2_1213.csv'
    'TA130_fastedSucrose2_1213.csv'
    
    %'TA120_fastedSucralose2_1215.csv'
    %'TA121_fastedSucralose2_1215.csv'
    %'TA124_fastedSucralose2_1215.csv'
    %'TA130_fastedSucralose2_1215.csv'
};
for i = 1:length(Davis_file)
    %% Read in licks per presentation
    davis = csvread(Davis_file{i}, 11, 8, [11 8 46 8]);
    time = csvread(Davis_file{i},48,1);
    I = find(time==0); time(I) =NaN;
    %% Sort licks
    if i==1
        Suc32 = []; Suc16 = []; Suc8 = []; Suc4 = []; Suc2=[]; Suc1=[]; Suc0_5=[];
        Suc0_25=[]; Suc0_125=[]; Water=[];
    end
    
    for j=1:36
        if sum(isnan(time(j))) == 0 %ignore trials with only one lick
            if j==9 || j==21 || j==33
                Suc32 = [Suc32; davis(j)/nansum(time(j,:))];
            elseif j==1 || j==13 || j==25
                Suc16 = [Suc16; davis(j)/nansum(time(j,:))];
            elseif j==2 || j==14 || j==26
                Suc8 = [Suc8; davis(j)/nansum(time(j,:))];
            elseif j==4 || j==16 ||j==28
                Suc4 = [Suc4; davis(j)/nansum(time(j,:))];
            elseif j==11 || j==23 || j==25
                Suc2 = [Suc2; davis(j)/nansum(time(j,:))];
            elseif j==5 || j==17 || j==29
                Suc1 = [Suc1; davis(j)/nansum(time(j,:))];
            elseif j==12 ||j==24 ||j==36
                Suc0_5 = [Suc0_5; davis(j)/nansum(time(j,:))];
            elseif j==8 || j==20 || j==32
                Suc0_25 = [Suc0_25; davis(j)/nansum(time(j,:))];
            elseif j==6||j==18||j==30
                Suc0_125 = [Suc0_125; davis(j)/nansum(time(j,:))];
            else
                Water = [Water; davis(j)/nansum(time(j,:))];
            end
        end
    end
end

%prism_input = [Suc32 Suc16 Suc8 Suc4 Suc2 Suc1 Suc0_5 Suc0_25 Suc0_125 Water];