%%
%%%%%%%%%%%%%%%%%%%%%%%%% Reporter pathway analysis %%%%%%%%%%%%%%%%%%%%%%%

%%

clear
load('reaction_description.mat')
load('Recon2.v04.mat')
load('flux_p_val.mat')
load('up_reaction.mat')
load('down_reaction.mat')

Z=norminv(1-flux_p_val);
Z(Z==Inf)=0;
Pathways=unique(reaction_description(:,5));

h=tabulate(reaction_description(:,5));
up=up_reaction{1,1};
h_u=tabulate(reaction_description(up,5));
down=down_reaction{1,1};
h_d=tabulate(reaction_description(down,5));

Z_score_pathway=[];
pathway_meta={};
X=strings(length(Pathways),5);
for i=1:length(Pathways)
    r=find(reaction_description(:,5)==Pathways(i));
    Z_score_pathway(i,1)=sum(Z(r))/sqrt(length(r));
    
    M=full(modelR204.S(:,r));
    M(find(M))=1;
    M=sum(M,2);
    pathway_meta{i,1}=Pathways{i};
    pathway_meta{i,2}=modelR204.mets(find(M));
    pathway_meta{i,3}=modelR204.metNames(find(M));
    
    P_h=find(string(h(:,1))==Pathways(i));
    P_hu=find(string(h_u(:,1))==Pathways(i));
    P_hd=find(string(h_d(:,1))==Pathways(i));
    X(i,1)=Pathways(i);
    X(i,2)=h{P_h,2};
    if ~isempty(P_hu)
    X(i,3)=h_u{P_hu,2};
    else
      X(i,3)=0;
    end
     if ~isempty(P_hd)
    X(i,4)=h_d{P_hd,2};
    else
      X(i,4)=0;
    end      
    X(i,5)=double(X(i,3))+double(X(i,4));
end
% save pathway_meta pathway_meta
% save Z_score_pathway Z_score_pathway

p1=1-normcdf(Z_score_pathway);
pathway_p_val=p1;
pathway_p_val(p1==0)= min(p1(p1~=0));
Pathway_enrich=[X pathway_p_val];
% save Pathway_enrich Pathway_enrich

% 1st column ---> Name of the pathways
% 2nd column ---> Number of reactions invloved in that pathway
% 3rd column ---> Number of up-regulated reactions in that pathway
% 4th column ---> Number of down-regulated reactions in that pathway
% 5th column ---> Number of altered reactions in that pathway
% 6th column ---> Obtained p-value

