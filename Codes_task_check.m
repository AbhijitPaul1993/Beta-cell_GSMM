%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Metabolic task check %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

% Import "Tasks.xlsx" as string array
Tasks(:,7:12) = [];

load('reaction_description.mat')
[a,b,c]=intersect(reaction_description(:,1),Tasks(:,3));

Check=[];
for i=1:10
    load(['model_d_' num2str(i) '.mat']);
    
    for j=1:length(a)
        eval(sprintf('model=model_d_%d;', i))
        model.c(b(j))=1;
        sol1=optimizeCbModel(model);
        sol2=optimizeCbModel(model,'min');
        Check(c(j),i)=max(abs(sol1.f),abs(sol2.f));
    end
end

for i=11:20
    load(['model_nd_' num2str(i-10) '.mat']);
    
    for j=1:length(a)
        eval(sprintf('model=model_nd_%d;', i-10))
        model.c(b(j))=1;
        sol1=optimizeCbModel(model);
        sol2=optimizeCbModel(model,'min');
        Check(c(j),i)=max(abs(sol1.f),abs(sol2.f));
    end
end

meta_task_check= [Tasks Check];
meta_task_check(2,7:16)="Diabetes";
meta_task_check(2,17:end)="Control";
% save meta_task_check meta_task_check