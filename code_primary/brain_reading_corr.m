
% load('network.mat','nos','fa')
load("network_pre_nos.mat", "label", "allh",'all_thres','all_bin')
% load('group_net.mat','group_l')
load('bh.mat') %% sharing behavioral and demographic data has not been approved by ethics at Western University
load('sub.mat')
load("data\FA_mtx_raw.mat");

%load('roi.mat')

ROI_label =  [77;11;13;17;29;55;63;65;79;81;83;85;89];
%ROI_label =  [11;13;29;55;63;65;77;79;81;85];

all_subnet = cell(1,length(sub));
fa_subnet = zeros(length(sub),1);
for i = 1:length(sub)
   group_l = allh{20}{i};
   ROI = [ROI_label;setdiff(group_l',ROI_label)];
   all_subnet{1,i} = all_thres{i}(ROI,ROI);
   f_sub = FA_mtx_raw(:,:,i);

   f_sub = weight_conversion(f_sub,'autofix');%map fa net
   f_sub = f_sub .* all_bin{i};
   AA = triu(f_sub(ROI,ROI));
   fa_subnet(i,1) = mean(AA(AA~=0));
% fa_subnet(i,1) = sum(AA(:));
%     fa_subnet(i,1) = sum(triu(f_sub),'all')/sum(1:length(f_sub)-1);
end
save("fa_subnet.mat","fa_subnet")
%%
%topological properties
thala_mean = zeros(length(sub),1);
cc_all = cell(length(sub),1);
cc_nrm_all = cell(length(sub),1);
den_all = zeros(length(sub),1);
D_all = cell(1,length(sub));
hops_all = cell(1,length(sub));
Pmat_all = cell(1,length(sub));
cpl_all = zeros(length(sub),1);
cpl_nrm_all = zeros(length(sub),1);
dgr_all = cell(length(sub),1);
str_all = cell(length(sub),1);
pth_tran_all = cell(1,length(sub));
% btc_all = zeros(length(sub),length(ROI));
% Ebtc_all = cell(length(sub),length(ROI));
Q1_all = zeros(length(sub),1);
% M_all = zeros(length(sub),length(ROI));
% Z_all = zeros(length(sub),length(ROI));
% P_all = zeros(length(sub),length(ROI));
% Pnorm_all = zeros(length(sub),length(ROI));
% Presidual_all = zeros(length(sub),length(ROI));
% betweenmodk_all = zeros(length(sub),length(ROI));
Eglob_all = zeros(length(sub),1);
Eloc_all =cell(length(sub),1);
L_all = cell(1,length(sub));

for i = 1:length(sub)
    a = all_subnet{i};   %thresholded network
    th = a(1,:);
    thala_mean(i,1) = sum(th(th ~= 0));
    b = weight_conversion(a,'binarize');
    a_rand = randmio_und(a,1);
    b_rand = weight_conversion(a_rand,'binarize');
    a_nrm = weight_conversion(a,'normalize');  % normalized network
    a_rand_nrm = weight_conversion(a_rand,'normalize');
    cc = clustering_coef_wu(a); %clustering coefficient
    cc = cc';
    cc_rand = clustering_coef_wu(a_rand);
    cc_rand = cc_rand';
    cc_nrm = zeros(1,length(a));
    for j = 1:length(a)
        cc_nrm(j) = cc(j)/cc_rand(j);
    end
    cc_all{i} = cc;
    cc_nrm_all{i} = cc_nrm;
    %cc_nrm_all(isinf(cc_nrm_all)|isnan(cc_nrm_all)) = 0;
    
    den = density_und(a);  %density
    den_all(i) = den;
    
    L = weight_conversion(a,'lengths'); %characteristic path length
    L_all{i} = L;
    L_rand = randmio_und(L,1);
    %D = distance_wei(L);
    %load('~/Documents/DTI_RD/network_1length(ROI)3/AALdti/length/network_pre_len.mat','all_thres_l')
    %LL = all_thres_l{length(ROI)};
    [D,hops,Pmat] = distance_wei_floyd(L); %Pmat for retrieve shortest path
    D_all{i} = D;
    hops_all{i} = hops;
    Pmat_all{i} = Pmat;
    [D_rand,hops_rand,Pmat_rand] = distance_wei_floyd(L_rand);
    % cpl = charpath(D,0,0);
    % cpl_rand = charpath(D_rand,0,0);
    % cpl_nrm = cpl/cpl_rand;
    % cpl_all(i) = cpl;
    % cpl_nrm_all(i) = cpl_nrm;
    
    dgr = degrees_und(a); %degree
    dgr_all{i} = dgr;
    str = strengths_und(a); %strength
    str_all{i}= str;
    
    %Cs = subgraph_centrality(a); %subgraph_centrality
    
    %pth_tran = path_transitivity(L); %path_transitivity
    %pth_tran_all{i} = pth_tran;
    
    %binn = all_bin{i};
    btc = betweenness_wei(L); %betweeness centrality
%     btc = btc ./ (dgr' .* dgr');
%     btc = sqrt(btc);
    btc = btc';
%     btc_all(i,:) = btc;
    Ebtc = edge_betweenness_wei(L); %edge betweeness centrality
%     Ebtc_all{i} = Ebtc;
    
%     n = size(a,1);   %optimal community structure (modularity with iterations)
%     M = 1:n;
%     Q0 = -1; Q1 = 0;
%     while Q1-Q0 > 1e-5
%         Q0 = Q1;
%         [M,Q1] = community_louvain(a,[],M);
%     end
%     Q1_all(i) = Q1;
%     Z = module_degree_zscore(a,M,0); %degree centrality z-score
%     Z = Z';
%     Z_all(i,:) = Z;
%     [PC_norm,PC_residual,PC,between_mod_k] = participation_coef_norm(a,M,1000); %participation coefficient
%     PC_norm = PC_norm';
%     PC_residual = PC_residual';
%     PC = PC';
%     between_mod_k = between_mod_k';
%     P_all(i,:) = PC;
%     Pnorm_all(i,:) = PC_norm;
%     Presidual_all(i,:) = PC_residual;
%     betweenmodk_all(i,:) = between_mod_k;
%     MM = M;
%     M = M';
%     M_all(i,:) = M;
    
    
    
    Eglob = efficiency_wei(a); %global efficiency
    Eglob_all(i) = Eglob;
    Eloc = efficiency_wei(a,2); %local efficiency
    Eloc = Eloc';
    Eloc_all{i} = Eloc;
end

cc_thala = zeros(length(sub),1);
Eloc_thala = zeros(length(sub),1);
for p = 1:length(sub)
    cc_thala(p,1) = cc_all{p}(1);
    Eloc_thala(p,1)= Eloc_all{p}(1);
end

% gender_label = cova(:,1);
% gender_male = gender_label == 0;
% gender_female = gender_label == 1;
% load("ov.mat")
aaa = [cc_thala,Eloc_thala,bhnew,cova_0(:,[1,2,4]),fa_subnet];


[cc, cp] = partialcorr(aaa(:,1:6),aaa(:,[7,8,10]))
% p = polyfit(cc_thala,bh(:,3),1);
% yfit = polyval(p,cc_thala);
% 
% plot(cc_thala,bh(:,3),'bo',cc_thala,yfit,'k-');

%scatter(cc_thala,bh(:,3))
% [ccm,cpm] = partialcorr([thala_count(gender_male == 1,:),cc_all(gender_male == 1,18),thala_mean(gender_male == 1,1),Eloc_all(gender_male == 1,18),bhnew(gender_male == 1,:)],[cova(gender_male == 1,3),str_fa_sub(gender_male == 1,1)])
% [ccf,cpf] = partialcorr([thala_count(gender_female == 1,:),cc_all(gender_female == 1,18),thala_mean(gender_female == 1,1),Eloc_all(gender_female == 1,18),bhnew(gender_female == 1,:)],[cova(gender_female == 1,3),str_fa_sub(gender_female == 1,1)])


