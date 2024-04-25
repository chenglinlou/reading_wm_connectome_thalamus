%efficiency corr run
net_pre_nos(3);
brain_reading_corr

%% cost corr run
main_BRW_LAU2_glob_lambda;
%%
load("bh.mat")
load("routing_results_main_thresed.mat","Ctrans")
load("fa_subnet.mat","fa_subnet")
info_all = zeros(96,30);
for j = 2:31
    [a,b,cc,cp] = cost_read_corr(j);
    info_all(:,j-1) = a+b;
end
%%
info = mean(info_all,2);
[c, p] = partialcorr([info,bh(:,1:5)],[bh(:,[7]),fa_subnet])
