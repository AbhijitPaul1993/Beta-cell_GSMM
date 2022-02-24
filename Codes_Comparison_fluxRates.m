%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Comparison between the flux distributions

% load('rxn_exp.mat')
% A=rxn_exp;
% min(min(A))/max(max(rxn_exp))=1.5584e-07;

clear
load('Recon2.v04.mat')
load('reaction_description')
load('flux_dia.mat')
load('flux_non_dia.mat')


A=flux_non_dia;
B=flux_dia;
A(find(abs(A)<10^-8))=0;
B(find(abs(B)<10^-8))=0;

flux_fold=median(B,2)./median(A,2);
rng('default')
flux_p_val = mattest(B, A, 'permute', true);

% p-value cut-off
p=find(flux_p_val <=0.05 & flux_p_val > 0);
p_r=p(union(find(flux_fold(p)<0),find(flux_fold(p)==Inf)));
p=setdiff(p,p_r);

% up-regulated altered reactions
p1=p(find(flux_fold(p)>1.2)); %248
[a1,b1]=sort(modelR204.subSystems(p1));
% [flux_dia(p1(b1),:) zeros(length(p1),1) flux_non_dia(p1(b1),:)];
up_reaction{1,1}=p1(b1);
up_reaction{1,2}=flux_fold(p1(b1));
up_reaction{1,3}=flux_p_val(p1(b1));
up_reaction{1,4}=reaction_description(p1(b1),:);

% down-regulated altered reactions
p2=p(find(flux_fold(p)<1/1.2)); %187
[a2,b2]=sort(modelR204.subSystems(p2));
% [flux_dia(p2(b2),:) zeros(length(p2),1) flux_non_dia(p2(b2),:)];
down_reaction{1,1}=p2(b2);
down_reaction{1,2}=flux_fold(p2(b2));
down_reaction{1,3}=flux_p_val(p2(b2));
down_reaction{1,4}=reaction_description(p2(b2),:);

% save up_reaction up_reaction
% save down_reaction down_reaction