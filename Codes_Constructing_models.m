%%%%%%%%%%%%%%%%%%%%%%% Bounds of metabolic reactions %%%%%%%%%%%%%%%%%%%%%
%% Bound of reactions using the gene expression data

clear
initCobraToolbox()

% load generic metabolic model
load('Recon2.v04.mat')

% load gene expresssion value of metabolic genes
load('meta_gene_exp.mat')
% first ten columns respresent the expresssion data of control subjects
% andlater ten columns respresent the expresssion data of diabetic subjects

%% Calculating reactions flux bounf by using E-flux method

rxn_exp=zeros(length(modelR204.rxns),size(meta_gene_exp,2));
for i=1:length(modelR204.rxns)
    if ~isempty(modelR204.grRules{i})
        ans1=strsplit(modelR204.grRules{i},'or');
        A=zeros(length(ans1),size(meta_gene_exp,2));
        for j=1:length(ans1)
            ans2=strfind(ans1{j},'and');
            if isempty(ans2)
                ans3=ans1{j}(setdiff(1:length(ans1{j}),union(union(strfind(ans1{j},' '),strfind(ans1{j},')')),strfind(ans1{j},'('))));
                [a,b]=intersect(modelR204.genes,ans3);
                if length(b)==0
                    A(j,:)=0;
                else
                    A(j,:)=meta_gene_exp(b,:);
                end
            else
                ans4=strsplit(ans1{j},'and');
                B=zeros(length(ans4),size(meta_gene_exp,2));
                for k=1:length(ans4)
                    ans5=ans4{k}(setdiff(1:length(ans4{k}),union(union(strfind(ans4{k},'('),strfind(ans4{k},')')),strfind(ans4{k},' '))));
                    [a,b]=intersect(modelR204.genes,ans5);
                    B(k,:)=meta_gene_exp(b,:);
                end
                for m=1:size(B,2)
                    if length(find(B(:,m)))==0
                        A(j,m)=0;
                    else
                        A(j,m)=min(B(find(B(:,m)),m));
                    end
                end
            end
        end
        rxn_exp(i,:)=sum(A,1);
    end
end
% save rxn_exp rxn_exp


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GSMMs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constructing the subject-specific models by fixing the bounds of reactions and with no objective function


load('Recon2.v04.mat')
load('rxn_exp.mat')

A=rxn_exp/max(max(rxn_exp));

for i=1:10
    model=modelR204;
    model.c(find(model.c))=0;
    model.ub(find(A(:,i)))=A(find(A(:,i)),i);
    model.ub(find(A(:,i)==0))=1;
    
    model.lb(intersect(find(model.lb),find(A(:,i))))=-A(intersect(find(model.lb),find(A(:,i))),i);
    model.lb(setdiff(find(model.lb),find(A(:,i))))=-1;
    
    eval(sprintf('model_nd_%d = model;', i))
    eval(['save ' ['model_nd_',num2str(i)] ' ' ['model_nd_',num2str(i)] ''])
end

for i=11:20
    model=modelR204;
    model.c(find(model.c))=0;
    model.ub(find(A(:,i)))=A(find(A(:,i)),i);
    model.ub(find(A(:,i)==0))=1;
    
    model.lb(intersect(find(model.lb),find(A(:,i))))=-A(intersect(find(model.lb),find(A(:,i))),i);
    model.lb(setdiff(find(model.lb),find(A(:,i))))=-1;
    
    eval(sprintf('model_d_%d = model;', i-10))
    eval(['save ' ['model_d_',num2str(i-10)] ' ' ['model_d_',num2str(i-10)] ''])
    
end
