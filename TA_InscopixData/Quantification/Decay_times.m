% Finding the time to half max for each neuron
% Run Heatmap_quant beforehand

%If don't have neurons_act, put a line here to make it so.
Neurons_act = mean(all_mrerank,1);
t_halfmax = []; %storage vector for the output
for i = 1:size(Neurons_act)
   avg_trace = movmean(Neurons_act(i,:), 50); %make trace less noisy
   halfmax = max(avg_trace(600:end))/2;
   index = find(avg_trace(600:end)>=halfmax, 1, 'last');
   t_halfmax = [t_halfmax; index];
end
%% For Davis Curves

Decay.Suc32 = []; Decay.Suc16 = []; Decay.Suc8=[];Decay.Suc4=[];
Decay.Suc2=[];Decay.Suc1=[];Decay.Suc05=[];Decay.Suc025=[];Decay.Suc0125=[];
Decay.Water=[];

for i = 1:size(psth.Suc32,1)
    halfmax = max(psth.Suc32(i,30:end))/2;
    index = find(psth.Suc32(i,30:end)>=halfmax,1, 'last');
    Decay.Suc32 = [Decay.Suc32; index];
end
for i = 1:size(psth.Suc16,1)
    halfmax = max(psth.Suc16(i,30:end))/2;
    index = find(psth.Suc16(i,30:end)>=halfmax,1, 'last');
    Decay.Suc16 = [Decay.Suc16; index];
end
for i = 1:size(psth.Suc8,1)
    halfmax = max(psth.Suc8(i,30:end))/2;
    index = find(psth.Suc8(i,30:end)>=halfmax,1, 'last');
    Decay.Suc8 = [Decay.Suc8; index];
end    
for i = 1:size(psth.Suc4,1)
    halfmax = max(psth.Suc4(i,30:end))/2;
    index = find(psth.Suc4(i,30:end)>=halfmax,1, 'last');
    Decay.Suc4 = [Decay.Suc4; index];
end
for i = 1:size(psth.Suc2,1)
    halfmax = max(psth.Suc2(i,30:end))/2;
    index = find(psth.Suc2(i,30:end)>=halfmax,1, 'last');
    Decay.Suc2 = [Decay.Suc2; index];
end
for i = 1:size(psth.Suc1,1)
    halfmax = max(psth.Suc1(i,30:end))/2;
    index = find(psth.Suc1(i,30:end)>=halfmax,1, 'last');
    Decay.Suc1 = [Decay.Suc1; index];
end
for i = 1:size(psth.Suc05,1)
    halfmax = max(psth.Suc05(i,30:end))/2;
    index = find(psth.Suc05(i,30:end)>=halfmax,1, 'last');
    Decay.Suc05 = [Decay.Suc05; index];
end
for i = 1:size(psth.Suc025,1)
    halfmax = max(psth.Suc025(i,30:end))/2;
    index = find(psth.Suc025(i,30:end)>=halfmax,1, 'last');
    Decay.Suc025 = [Decay.Suc025; index];
end
for i = 1:size(psth.Suc0125,1)
    halfmax = max(psth.Suc0125(i,30:end))/2;
    index = find(psth.Suc0125(i,30:end)>=halfmax,1, 'last');
    Decay.Suc0125 = [Decay.Suc0125; index];
end
for i = 1:size(psth.Water,1)
    halfmax = max(psth.Water(i,30:end))/2;
    index = find(psth.Water(i,30:end)>=halfmax,1, 'last');
    Decay.Water = [Decay.Water; index];
end
    