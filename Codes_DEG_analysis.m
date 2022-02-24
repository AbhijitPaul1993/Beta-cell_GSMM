%% Differentially expressed metabolic genes

clear

% load generic metabolic model
load('Recon2.v04.mat')

% load gene expresssion value of metabolic genes
load('meta_gene_exp.mat')
% first ten columns respresent the expresssion data of control subjects
% andlater ten columns respresent the expresssion data of diabetic subjects

A=meta_gene_exp(:,1:10);
B=meta_gene_exp(:,11:20);
meta_gene_fold=median(B,2)./median(A,2);

rng('default')
meta_gene_p_val= mattest(B, A, 'permute', true);
% save meta_gene_fold meta_gene_fold
% save meta_gene_p_val meta_gene_p_val

p=find(meta_gene_p_val <=0.05 & meta_gene_p_val > 0);

p1=p(find(meta_gene_fold(p)>1.2));
[a1,b1]=sort(meta_gene_fold(p1),'descend');
meta_up_gene{1,1}=p1(b1);
meta_up_gene{1,2}=a1;
meta_up_gene{1,3}=meta_gene_p_val(p1(b1));
meta_up_gene{1,4}=modelR204.genes(p1(b1));

p2=p(find(meta_gene_fold(p)<1/1.2));
[a2,b2]=sort(meta_gene_fold(p2));
meta_down_gene{1,1}=p2(b2);
meta_down_gene{1,2}=a2;
meta_down_gene{1,3}=meta_gene_p_val(p2(b2));
meta_down_gene{1,4}=modelR204.genes(p2(b2));


