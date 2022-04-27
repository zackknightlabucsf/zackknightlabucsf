function [A_combined_psth_start, A_combined_psth_end, B_combined_psth_start, B_combined_psth_end] = rerank_for_inscopix_analysis_JG(all_psth_start, all_psth_end, sorted_by)

%make structures
A_combined_psth_start=struct('avgpsth',[],'dff',[], 'dff_total',[], 'zC', [], 'experiment', []);
B_combined_psth_start=struct('avgpsth',[],'dff',[], 'dff_total',[], 'zC', [], 'experiment', []);
A_combined_psth_end=struct('avgpsth',[],'dff',[], 'dff_total',[], 'zC', [], 'experiment', []);
B_combined_psth_end=struct('avgpsth',[], 'dff',[], 'dff_total',[], 'zC', [], 'experiment', []);

%make A structure
whichref=find(sum(strcmp(struct2cell(all_psth_start), 'A'),1)==1);
f = fieldnames(A_combined_psth_start);
for i = 1:length(f)
    for r = 1:length(whichref)
        A_combined_psth_start.(f{i}) = [A_combined_psth_start.(f{i});all_psth_start(whichref(r)).(f{i})];
        A_combined_psth_end.(f{i}) = [A_combined_psth_end.(f{i}); all_psth_end(whichref(r)).(f{i})];
    end
end

%make B structure
whichref=find(sum(strcmp(struct2cell(all_psth_start), 'B'),1)==1);
f = fieldnames(B_combined_psth_start);
for i = 1:length(f)
    for r = 1:length(whichref)
        B_combined_psth_start.(f{i}) = [B_combined_psth_start.(f{i});all_psth_start(whichref(r)).(f{i})];
        B_combined_psth_end.(f{i}) = [B_combined_psth_end.(f{i}); all_psth_end(whichref(r)).(f{i})];
    end
end

%sorting id
if(sorted_by== "A start")
    [~,rerank]=sort(A_combined_psth_start.dff);
end

if(sorted_by== "B start")
    [~,rerank]=sort(B_combined_psth_start.dff);
end

%sorted by total dff
if(sorted_by== "A total")
    [~,rerank]=sort(A_combined_psth_start.dff_total);
end
if(sorted_by== "B total")
    [~,rerank]=sort(B_combined_psth_start.dff_total);
end



%sort
A_combined_psth_start.dff=A_combined_psth_start.dff(rerank);
if(size(A_combined_psth_start.avgpsth,1)>0)
    A_combined_psth_start.avgpsth=A_combined_psth_start.avgpsth(rerank,:);
end
if(size(A_combined_psth_start.zC,1)>0)
    A_combined_psth_start.zC=A_combined_psth_start.zC(rerank,:);
end

A_combined_psth_end.dff=A_combined_psth_end.dff(rerank);
if(size(A_combined_psth_end.avgpsth,1)>0)
    A_combined_psth_end.avgpsth=A_combined_psth_end.avgpsth(rerank,:);
end
if(size(A_combined_psth_end.zC,1)>0)
    A_combined_psth_end.zC=A_combined_psth_end.zC(rerank,:);
end

B_combined_psth_start.dff=B_combined_psth_start.dff(rerank);
if(size(B_combined_psth_start.avgpsth,1)>0)
    B_combined_psth_start.avgpsth=B_combined_psth_start.avgpsth(rerank,:);
end
if(size(B_combined_psth_start.zC,1)>0)
    B_combined_psth_start.zC=B_combined_psth_start.zC(rerank,:);
end

B_combined_psth_end.dff=B_combined_psth_end.dff(rerank);
if(size(B_combined_psth_end.avgpsth,1)>0)
    B_combined_psth_end.avgpsth=B_combined_psth_end.avgpsth(rerank,:);
end
if(size(B_combined_psth_end.zC,1)>0)
    B_combined_psth_end.zC=B_combined_psth_end.zC(rerank,:);
end

end