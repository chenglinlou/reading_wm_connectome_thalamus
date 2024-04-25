clear all
close all
clc

path2code = fullfile(pwd);
%addpath('~/mat_tools/BCT_2017'); % add BCT to path
load("network_pre_nos.mat","all_thres") % load SC from .mat file; SC is NxNxM - N nodes, M subjects
load('fa_subnet.mat')
load('sub.mat')
% load("LAU2_scale3.mat","SC")
SC = zeros(90,90,64);
for i = 1:64
    SC(:,:,i) = all_thres{1,i};
    % SC(:,:,i) = randmio_und(all_thres{1,i},500);
end
SCall = SC; %clear SC;

[N,~,M] = size(SCall);

lambda_vals = logspace(1.45,-2.2,30);
lambda_vals = [0,lambda_vals(end:-1:1)];

model_ops = {'SPL_W_log','SPL_W_inv'}; % weight transforms 

model = model_ops{1};
rho = [.15]; % optional density parameter to ensure same density across all subjects

path2results = fullfile(pwd,'..','results',num2str(rho),model);
if ~exist(path2results,'dir')
    mkdir(path2results);
end

mask = find(~eye(N));
KLref = zeros(N,N,M,length(lambda_vals));
Cinfo = zeros(N,N,M,length(lambda_vals));
Ctrans = zeros(N,N,M,length(lambda_vals));
visits = zeros(N,N,M,length(lambda_vals));
for s = 1:M
    
    SC = SCall(:,:,s);
    SC = weight_conversion(SC,'normalize');
%     SC = randmio_und(SC,500);
    SC(eye(N)>0) = 0;
    
    thr = prctile(SC(mask),(1-rho)*100); % threshold SC to desired dens
    SC(SC <= thr) = 0;
    if length(unique(get_components(SC))) > 1
        error('\n Network %d is disconnected \n',s)
    end   
    
    k = sum(SC,1);
    k = k(ones(N,1),:).*eye(N,N);
    Pref = k\SC;
    
    % calculate Local and Global matrices
    [G,D,B,steps] = f_local_global(SC,[],model);
%     Wn = weight_conversion(SC,'normalize');
%     D = weight_conversion(Wn,'lengths');
%     [G,~,~] = distance_wei_floyd(D);
    
%     KLref = zeros(N,N,length(lambda_vals));
%     Cinfo = zeros(N,N,length(lambda_vals));
%     Ctrans = zeros(N,N,length(lambda_vals));
%     visits = zeros(N,N,length(lambda_vals));
    
    tic;
    for lind = 1:length(lambda_vals)
        fprintf('lambda = %f \n',lambda_vals(lind));
        for t = 1:N
            % define G(t) - global term
            Gt = G(:,t)';
            Gt = Gt(ones(N,1),:) + D;
            
            % prob of going from i to j, given a target t
            FG = exp(-((lambda_vals(lind)*Gt) + D)).*(SC>0);
            Z_lamb_t = sum(FG,2);
            Z_lamb_t = Z_lamb_t(:,ones(1,N)); % Normalize to ensure sum(P,2) = ones(N,1)
            P = FG./Z_lamb_t;
            P(isnan(P))=0;
%             if sum(sum(P,2),1) ~= N  % sanity check 
%                 error(sprintf('\n ERROR -- subj %d \t lind %d \t target %d',s,lind,t));
%             end
            % KL divergence
            if lind > 1
                kk = P.*(log2(P) - log2(Pref));
                kk(isnan(kk)) = 0; kk(isinf(kk)) = 0;
                KLref(:,t,s,lind) = sum(kk,2);
            end
            % calculate mean absorbtion times
            Q = P;
            Q(t,:) = [];  % make target an absorbing state
            Q(:,t) = [];
            FM  = inv(eye(N-1) - Q);    % Fundamental Matrix
            ind = true(N,1); ind(t) = false;
            Cinfo(ind,t,s,lind) = (bsxfun(@(a,b) a./b,FM,sum(FM,2)))*KLref(ind,t,s,lind);  %line 83 sum of KLref
            Ctrans(ind,t,s,lind) = FM*(sum(P(ind,:).*(D(ind,:)),2));   %sum
            visits(ind,t,s,lind) = mean(FM,1)';
            
        end
    end
    %fname = fullfile(path2results,sprintf('BRW_%d.mat',s));   
    toc
end
save('routing_results_main_thresed.mat','lambda_vals','KLref','Cinfo','Ctrans','visits')
% lambda_vals_rand = lambda_vals; KLref_rand = KLref; Cinfo_rand = Cinfo; Ctrans_rand = Ctrans;visits_rand = visits;
% save('routing_results_main_thresed_rand.mat','lambda_vals_rand','KLref_rand','Cinfo_rand','Ctrans_rand','visits_rand')

% tran_t = zeros(64,90,31);
% tran_s = zeros(64,90,31);
% info_t = zeros(64,90,31);
% info_s = zeros(64,90,31);
% 
% for i = 1:31
%     atran = Ctrans(:,:,:,i);
%     ainfo = Cinfo(:,:,:,i);
%     for j = 1:64
%         aatran = atran(:,:,j);
%         aainfo = ainfo(:,:,j);
%         tran_t(j,:,i) = sum(aatran,1)';
%         tran_s(j,:,i) = sum(aatran,2);
%         info_t(j,:,i) = sum(aainfo,1)';
%         info_s(j,:,i) = sum(aainfo,2);
%     end
% end
