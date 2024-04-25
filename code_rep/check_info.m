function [imm_info, imm_tran_cost] = check_info(lmd)
routing_strg = routing(lmd);
P = routing_strg{1,1};
P_t = P(:,:,15);
load("../network_pre_nos.mat","all_thres")
load("../graph_para.mat","L_all")
L = L_all{1,1};
SC = all_thres{1,1};
%     SC(eye(90)>0) = 0;
%     kkk = sum(SC,1); kkk(kkk==0) = 0.000001;
%     kkk = kkk(ones(90,1),:).*eye(90,90);
%     P_ref = kkk\SC;
kkk = sum(SC,2);
P_ref = SC./kkk;
    
kk = P_t.*(log2(P_t) - log2(P_ref));
        kk(isnan(kk)) = 0; kk(isinf(kk)) = 0;
        imm_info = sum(kk,2);
        d = log2(P_t) - log2(P_ref);
        d(isnan(d)) = 0; d(isinf(d)) = 0;
        imm_tran_cost = sum(d(:));
%  imm_tran_cost_j = P_t .* L;
%         imm_tran_cost = sum(imm_tran_cost_j,2);