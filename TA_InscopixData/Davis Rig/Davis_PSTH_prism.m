
Suc32 = [9 21 33];
Suc16 = [1 13 25];
Suc8 = [2 14 26];
Suc4 = [4 16 28];
Suc2 = [11 23 35];
Suc1 = [5 17 29];
Suc0_5 = [12 24 36];
Suc0_25 = [8 22 32];
Suc0_125 = [6 18 30];
Water = [3 7 10 15 19 22 27 31 34];

mice = 4;

psth.Suc32 = []; psth.Suc16=[]; psth.Suc8=[];psth.Suc4=[];psth.Suc2=[];psth.Suc1=[];
psth.Suc05=[];psth.Suc025=[];psth.Suc0125=[];psth.Water=[];
for i = 1:mice
    mouse = psth_all_act{i};
    for j = 1:length(Suc32)
        psth.Suc32 = [psth.Suc32; mouse{Suc32(j)}];
    end
    for j = 1:length(Suc16)
        psth.Suc16 = [psth.Suc16; mouse{Suc16(j)}];
    end
    for j = 1:length(Suc8)
        psth.Suc8 = [psth.Suc8; mouse{Suc8(j)}];
    end
    for j = 1:length(Suc4)
        psth.Suc4 = [psth.Suc4; mouse{Suc4(j)}];
    end
    for j = 1:length(Suc2)
        psth.Suc2 = [psth.Suc2; mouse{Suc2(j)}];
    end
    for j = 1:length(Suc1)
        psth.Suc1 = [psth.Suc1; mouse{Suc1(j)}];
    end
    for j = 1:length(Suc0_5)
        psth.Suc05 = [psth.Suc05; mouse{Suc0_5(j)}];
    end
    for j = 1:length(Suc0_25)
        psth.Suc025 = [psth.Suc025; mouse{Suc0_25(j)}];
    end
    for j = 1:length(Suc0_125)
        psth.Suc0125 = [psth.Suc0125; mouse{Suc0_125(j)}];
    end
    for j = 1:length(Water)
        psth.Water = [psth.Water; mouse{Water(j)}];
    end
end

% %% for prism
% psth.Suc32_m = mean(psth.Suc32,1); psth.Suc32_s = std(psth.Suc32,0,1)./sqrt(size(psth.Suc32,1));
% psth.Suc16_m = mean(psth.Suc16,1); psth.Suc16_s = std(psth.Suc16,0,1)./sqrt(size(psth.Suc16,1));
% psth.Suc8_m = mean(psth.Suc8,1); psth.Suc8_s = std(psth.Suc8,0,1)./sqrt(size(psth.Suc8,1));
% psth.Suc4_m = mean(psth.Suc4,1); psth.Suc4_s = std(psth.Suc4,0,1)./sqrt(size(psth.Suc4,1));
% psth.Suc2_m = mean(psth.Suc2,1); psth.Suc2_s = std(psth.Suc2,0,1)./sqrt(size(psth.Suc2,1));
% psth.Suc1_m = mean(psth.Suc1,1); psth.Suc1_s = std(psth.Suc1,0,1)./sqrt(size(psth.Suc1,1));
% psth.Suc05_m = mean(psth.Suc05,1); psth.Suc05_s = std(psth.Suc05,0,1)./sqrt(size(psth.Suc05,1));
% psth.Suc025_m = mean(psth.Suc025,1); psth.Suc025_s = std(psth.Suc025,0,1)./sqrt(size(psth.Suc025,1));
% psth.Suc0125_m = mean(psth.Suc0125,1); psth.Suc0125_s = std(psth.Suc0125,0,1)./sqrt(size(psth.Suc0125,1));
% psth.Water_m = mean(psth.Water,1); psth.Water_s = std(psth.Water,0,1)./sqrt(size(psth.Water,1));