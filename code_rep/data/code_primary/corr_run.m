%efficiency corr run
net_pre_nos(3); %network preprocessing
brain_reading_corr % get brain-reading correlations

%% cost corr run
main_BRW_LAU2_glob_lambda; %extract routing strategy matrices

trans_all = zeros(64,30); 
for j = 2:31
    [a,b,cc,cp] = cost_read_corr(j);
    trans_all(:,j-1) = a+b;
end
trans = mean(trans_all,2); %mean cost across all lambda
[c, p] = partialcorr([info,bhnew],[cova(:,[1,3]),fa_subnet]);