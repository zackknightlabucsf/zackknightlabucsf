function [A_combined_psth_start, A_combined_psth_end, B_combined_psth_start, B_combined_psth_end] = rerank_for_inscopix_analysis_JG(all_psth_start, all_psth_end, sorted_by)

%make structures
A_combined_psth_start=struct('psth',{[]},'avgpsth',[],'ttl_bout',[],'time_ttl_bout',[],'dff',[],'boutsize',[], 'crossreg',[], 'mouseid',[], 'zC', [], 'experiment', []);
B_combined_psth_start=struct('psth',{[]},'avgpsth',[],'ttl_bout',[],'time_ttl_bout',[],'dff',[],'boutsize',[], 'crossreg', [], 'mouseid',[], 'zC', [], 'experiment', []);
A_combined_psth_end=struct('psth',{[]},'avgpsth',[],'ttl_bout',[],'time_ttl_bout',[],'dff',[],'boutsize',[], 'crossreg',[], 'mouseid',[], 'zC', [], 'experiment', []);
B_combined_psth_end=struct('psth',{[]},'avgpsth',[],'ttl_bout',[],'time_ttl_bout',[],'dff',[],'boutsize',[], 'crossreg', [], 'mouseid',[], 'zC', [], 'experiment', []);

%make A structure
whichref=find(sum(strcmp(struct2cell(all_psth_start), 'A'),1)==1);
f = fieldnames(A_combined_psth_start);
for i = 1:length(f)
    for r = 1:length(whichref)
        A_combined_psth_start(1).(f{i}) = all_psth_start(whichref).(f{i});
        A_combined_psth_end(1).(f{i}) = all_psth_end(whichref).(f{i});
    end
end

%make B structure
whichref=find(sum(strcmp(struct2cell(all_psth_start), 'B'),1)==1);
f = fieldnames(B_combined_psth_start);
for i = 1:length(f)
    for r = 1:length(whichref)
        B_combined_psth_start(1).(f{i}) = all_psth_start(whichref).(f{i});
        B_combined_psth_end(1).(f{i}) = all_psth_end(whichref).(f{i});
    end
end

%sorting id
if(sorted_by== "A start")
[~,rerank]=sort(A_combined_psth_start.dff);
end




%sort
A_combined_psth_start.dff=A_combined_psth_start.dff(rerank);
A_combined_psth_start.avgpsth=A_combined_psth_start.avgpsth(rerank,:);
A_combined_psth_start.zC=A_combined_psth_start.zC(rerank,:);

A_combined_psth_end.dff=A_combined_psth_end.dff(rerank);
A_combined_psth_end.avgpsth=A_combined_psth_end.avgpsth(rerank,:);
A_combined_psth_end.zC=A_combined_psth_end.zC(rerank,:);

B_combined_psth_start.dff=B_combined_psth_start.dff(rerank);
B_combined_psth_start.avgpsth=B_combined_psth_start.avgpsth(rerank,:);
B_combined_psth_start.zC=B_combined_psth_start.zC(rerank,:);

B_combined_psth_end.dff=B_combined_psth_end.dff(rerank);
B_combined_psth_end.avgpsth=B_combined_psth_end.avgpsth(rerank,:);
B_combined_psth_end.zC=B_combined_psth_end.zC(rerank,:);

end