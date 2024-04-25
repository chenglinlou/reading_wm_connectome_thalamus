function [info_target,info_source,cc,cp] = cost_read_corr(lamda)
load("bh.mat")
load('fa_subnet.mat')
% gender_label = cova(:,1);
% gender_male = gender_label == 0;
% gender_female = gender_label == 1;
load("network_pre_nos.mat","allh")
load('routing_results_main_thresed.mat','Cinfo','Ctrans',"visits")
ROI_label =  [77;11;13;17;29;55;63;65;79;81;83;85;89];
hublabel = allh{20};
group_l = cell(96,1);
for p = 1:96
    group_l{p} = [ROI_label;setdiff(hublabel{p}',ROI_label)];
end
%% cost of routing through the left thalamus (77th parcel in AAL)
a = Ctrans(:,:,:,lamda);
info_source = zeros(96,1);
info_target = zeros(96,1);
visit = visits(:,:,:,lamda);
visit_thala = zeros(96,1);
for i = 1:96
     label = group_l{i};
info = a(label,label,i);
iss = sum(info,2)';
info_source(i) = iss(:,1);
itt = sum(info,1);
info_target(i) = itt(:,1);
thala = visit(label,label,i);
b = sum(thala,1);
c = sum(thala,2);
visit_thala(i,1) = b(1)+c(1);
end
% par_exclude = [26,29,30,34,53,54,55,63,66,67,84,96,93,106,109,4,51,86];
% par_include = setdiff([1:112],par_exclude);
% [c, p] = partialcorr([info_source,visit_thala,bhnew],[cova(:,[1,3]),fa_subnet]);
[cc, cp] = partialcorr([info_target+info_source,bh],[bh(:,[6:7]),fa_subnet]);



% [c,p] = partialcorr([info_target(gender_male == 1,:) + info_source(gender_male == 1,:),bhnew(gender_male == 1,:)],[cova(gender_male == 1,3),str_fa_sub(gender_male == 1,1)]);
% [cc,cp] = partialcorr([info_target(gender_female == 1,:) + info_source(gender_female == 1,:),bhnew(gender_female == 1,:)],[cova(gender_female == 1,3),str_fa_sub(gender_female == 1,1)]);
%  [c,p] = partialcorr([sum(thala_trans_target+thala_trans_source,2)(gender_male == 1,:), sum(thala_info_target+thala_info_source,2)(gender_male == 1,:),bhnew(gender_male == 1,:)],cova(gender_male == 1,3));
