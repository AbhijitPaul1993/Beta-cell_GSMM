%%%%%%%%%%%%%%%%%%%%%%%%%%%% Flux distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Flux distribution of 10 diabeteic models

initCobraToolbox()

load('model_d_1.mat')
storedata=gpSampler(model_d_1,130000);
load('flux_dia.mat')
flux_dia(:,1)=mean(storedata.points,2);
% save flux_dia flux_dia

load('model_d_2.mat')
storedata=gpSampler(model_d_2,130000);
load('flux_dia.mat')
flux_dia(:,2)=mean(storedata.points,2);
% save flux_dia flux_dia

load('model_d_3.mat')
storedata=gpSampler(model_d_3,130000);
load('flux_dia.mat')
flux_dia(:,3)=mean(storedata.points,2);
% save flux_dia flux_dia

load('model_d_4.mat')
storedata=gpSampler(model_d_4,130000);
load('flux_dia.mat')
flux_dia(:,4)=mean(storedata.points,2);
% save flux_dia flux_dia

load('model_d_5.mat')
storedata=gpSampler(model_d_5,130000);
load('flux_dia.mat')
flux_dia(:,5)=mean(storedata.points,2);
% save flux_dia flux_dia

load('model_d_6.mat')
storedata=gpSampler(model_d_6,130000);
load('flux_dia.mat')
flux_dia(:,6)=mean(storedata.points,2);
% save flux_dia flux_dia

load('model_d_7.mat')
storedata=gpSampler(model_d_7,130000);
load('flux_dia.mat')
flux_dia(:,7)=mean(storedata.points,2);
% save flux_dia flux_dia

load('D:\Work\Diabetes_Beta_cell\Matlab\saturation_mean_flux\f_d_8_check.mat')
flux_dia(:,8)=f_d_8_check(:,2);
% save flux_dia flux_dia

load('model_d_9.mat')
storedata=gpSampler(model_d_9,130000);
load('flux_dia.mat')
flux_dia(:,9)=mean(storedata.points,2);
% save flux_dia flux_dia

load('model_d_10.mat')
storedata=gpSampler(model_d_10,130000);
load('flux_dia.mat')
flux_dia(:,10)=mean(storedata.points,2);
% save flux_dia flux_dia


%% Flux distribution of 10 non-diabeteic models

initCobraToolbox()

load('model_nd_1.mat')
storedata=gpSampler(model_nd_1,130000);
load('flux_non_dia.mat')
flux_non_dia(:,1)=mean(storedata.points,2);
% save flux_non_dia flux_non_dia

load('D:\Work\Diabetes_Beta_cell\Matlab\saturation_mean_flux\f_nd_2_check.mat')
flux_non_dia(:,2)=f_nd_2_check(:,2);
% save flux_non_dia flux_non_dia

load('model_nd_3.mat')
storedata=gpSampler(model_nd_3,130000);
load('flux_non_dia.mat')
flux_non_dia(:,3)=mean(storedata.points,2);
% save flux_non_dia flux_non_dia

load('model_nd_4.mat')
storedata=gpSampler(model_nd_4,130000);
load('flux_non_dia.mat')
flux_non_dia(:,4)=mean(storedata.points,2);
% save flux_non_dia flux_non_dia

load('model_nd_5.mat')
storedata=gpSampler(model_nd_5,130000);
load('flux_non_dia.mat')
flux_non_dia(:,5)=mean(storedata.points,2);
% save flux_non_dia flux_non_dia

load('model_nd_6.mat')
storedata=gpSampler(model_nd_6,130000);
load('flux_non_dia.mat')
flux_non_dia(:,6)=mean(storedata.points,2);
% save flux_non_dia flux_non_dia

load('model_nd_7.mat')
storedata=gpSampler(model_nd_7,130000);
load('flux_non_dia.mat')
flux_non_dia(:,7)=mean(storedata.points,2);
% save flux_non_dia flux_non_dia

load('model_nd_8.mat')
storedata=gpSampler(model_nd_8,130000);
load('flux_non_dia.mat')
flux_non_dia(:,8)=mean(storedata.points,2);
% save flux_non_dia flux_non_dia

load('model_nd_9.mat')
storedata=gpSampler(model_nd_9,130000);
load('flux_non_dia.mat')
flux_non_dia(:,9)=mean(storedata.points,2);
% save flux_non_dia flux_non_dia

load('model_nd_10.mat')
storedata=gpSampler(model_nd_10,130000);
load('flux_non_dia.mat')
flux_non_dia(:,10)=mean(storedata.points,2);
% save flux_non_dia flux_non_dia
