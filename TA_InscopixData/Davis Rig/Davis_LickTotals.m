Davis_file = { %output file from Davis converted to csv
    %'TA120_fastedSucrose2_1213.csv'
    %'TA121_fastedSucrose2_1213.csv'
    %'TA124_fastedSucrose2_1213.csv'
    %'TA130_fastedSucrose2_1213.csv'
    
    %'TA120_fastedSucralose2_1215.csv'
    %'TA121_fastedSucralose2_1215.csv'
    %'TA124_fastedSucralose2_1215.csv'
    %'TA130_fastedSucralose2_1215.csv'
    
    %'0331TA121_LipidCurve.csv'
    %'0331TA124_LipidCurve.csv'
    %'0401TA175_LipidCurve.csv'
    %'0401TA176_LipidCurve.csv'
    %'0401TA186_LipidCurve.csv'
    
    '0405TA121_tastePanel.csv'
    '0405TA124_tastePanel.csv'
    '0406TA175_tastePanel.csv'
    '0406TA176_tastePanel.csv'
    
};
Curve = 0;
if Curve == 1
    for i = 1:length(Davis_file)
        davis = csvread(Davis_file{i}, 11, 8, [11 8 46 8]);
        
        if i==1
            Suc32 = []; Suc16 = []; Suc8 = []; Suc4 = []; Suc2=[]; Suc1=[]; Suc0_5=[];
            Suc0_25=[]; Suc0_125=[]; Water=[];
        end
        for j=1:36
            if j==9 || j==21 || j==33
                Suc32 = [Suc32; davis(j)];
            elseif j==1 || j==13 || j==25
                Suc16 = [Suc16; davis(j)];
            elseif j==2 || j==14 || j==26
                Suc8 = [Suc8; davis(j)];
            elseif j==4 || j==16 ||j==28
                Suc4 = [Suc4; davis(j)];
            elseif j==11 || j==23 || j==35
                Suc2 = [Suc2; davis(j)];
            elseif j==5 || j==17 || j==29
                Suc1 = [Suc1; davis(j)];
            elseif j==12 ||j==24 ||j==36
                Suc0_5 = [Suc0_5; davis(j)];
            elseif j==8 || j==20 || j==32
                Suc0_25 = [Suc0_25; davis(j)];
            elseif j==6||j==18||j==30
                Suc0_125 = [Suc0_125; davis(j)];
            else
                Water = [Water; davis(j)];
            end
        end
        Davis.Suc32(i) = mean(Suc32(end-2:end));
        Davis.Suc16(i) = mean(Suc16(end-2:end));
        Davis.Suc8(i) = mean(Suc8(end-2:end));
        Davis.Suc4(i) = mean(Suc4(end-2:end));
        Davis.Suc2(i) = mean(Suc2(end-2:end));
        Davis.Suc1(i) = mean(Suc1(end-2:end));
        Davis.Suc0_5(i) = mean(Suc0_5(end-2:end));
        Davis.Suc0_25(i) = mean(Suc0_25(end-2:end));
        Davis.Suc0_125(i) = mean(Suc0_125(end-2:end));
        Davis.Water(i) = mean(Water(end-8:end));
        
    end
    clearvars Suc32 Suc16 Suc8 Suc4 Suc2 Suc1 Suc0_5 Suc0_25...
        Suc0_125 Water
    Davis.prism = [Davis.Suc32' Davis.Suc16' Davis.Suc8'...
        Davis.Suc4' Davis.Suc2' Davis.Suc1' Davis.Suc0_5'...
        Davis.Suc0_25' Davis.Suc0_125' Davis.Water'];
else
    for i = 1:length(Davis_file)
        davis = csvread(Davis_file{i}, 11, 8, [11 8 46 8]);
        
        if i==1
            Ensure = []; Polycal = []; Intralipid = []; Sucralose = []; 
            Salt=[]; CitricAcid=[]; SiliconeOil=[];
            Sucrose=[]; Quinine=[]; Water=[];
        end
        for j=1:36
            if j==2 || j==14 || j==26
                Ensure = [Ensure; davis(j)];
            elseif j==3 || j==15 || j==27
                Polycal = [Polycal; davis(j)];
            elseif j==4 || j==16 || j==28
                Intralipid = [Intralipid; davis(j)];
            elseif j==5 || j==17 || j==29
                Sucralose = [Sucralose; davis(j)];
            elseif j==6 || j==18 || j==30
                Salt = [Salt; davis(j)];
            elseif j==7 || j==19 || j==31
                CitricAcid = [CitricAcid; davis(j)];
            elseif j==8 || j==20 || j==32
                SiliconeOil = [SiliconeOil; davis(j)];
            elseif j==10 || j==22 || j==34
                Sucrose = [Sucrose; davis(j)];
            elseif j==11||j==23||j==35
                Quinine = [Quinine; davis(j)];
            else
                Water = [Water; davis(j)];
            end
        end
        Davis.Ensure(i) = mean(Ensure(end-2:end));
        Davis.Polycal(i) = mean(Polycal(end-2:end));
        Davis.Intralipid(i) = mean(Intralipid(end-2:end));
        Davis.Sucralose(i) = mean(Sucralose(end-2:end));
        Davis.Salt(i) = mean(Salt(end-2:end));
        Davis.CitricAcid(i) = mean(CitricAcid(end-2:end));
        Davis.SiliconeOil(i) = mean(SiliconeOil(end-2:end));
        Davis.Sucrose(i) = mean(Sucrose(end-2:end));
        Davis.Quinine(i) = mean(Quinine(end-2:end));
        Davis.Water2(i) = mean(Water(end-8:end));
    end
    clearvars Ensure Polycal Intralipid Sucralose Salt CitricAcid...
        SiliconeOil Sucrose Quinine Water
    Davis.Prism = [Davis.Ensure' Davis.Polycal' Davis.Intralipid'...
        Davis.Sucralose' Davis.Salt' Davis.CitricAcid' ...
        Davis.SiliconeOil' Davis.Sucrose' Davis.Quinine' Davis.Water2'];
    
end
    