function net_pre_nos(thre)
%load("sub.mat")
all_orin = {};
all_fix = {};
all_bin = {};

for m = 1:64%length(sub)
    load("data\NOS_mtx_raw.mat");
    all_orin{m} = NOS_mtx_raw(:,:,m);
    W = weight_conversion(all_orin{m}, "autofix");
    all_fix{m} = W;
end

all_thres = {};

for k = 1:64%length(sub)
    W = all_fix{k};
    W_thr = threshold_absolute(W, thre);
     %for i = 1:90
     %   for j = 1:90
      %      v = W(i,j);
       %     if(v <= thre)
        %        W(i,j) = 0;
         %   end
      %  end
    %end
    Wb = weight_conversion(W_thr,"binarize");
    all_bin{k} = Wb;
    all_thres{k} = W_thr;
end

%rich-club coefficient and normalized rich-club coefficient
all_rand = {};
all_rnorm = {};
for i = 1:64%length(sub)
    Wb = all_thres{i};
    Wrand = randmio_und(Wb,5000);
    all_rand{i} = Wrand;
    r = rich_club_wu(Wb,50);
    rr = rich_club_wu(Wrand,50);
    rnorm = [];
    for j = 1:length(r)
        rnorm(j) = r(j)/rr(j);
    end
    all_rnorm{i} = rnorm;
end

for i = 1:64%length(sub)
    a = all_rnorm{i};
    a(isnan(a))=0;
    all_rnorm{i} = a;
end


 k = zeros(50,64);
for j = 1:50
kk = [];
for i = 1:64%length(sub)
a = all_rnorm{i};
b = a(j);
kk = [kk,b];
end
k(j,:) = kk;
end

K = k';
m_K = mean(K);

H = [];
P = [];
%CI = {};
for i = 1:50
x = k(i,:);
[hx, px] = ttest(x,1,0.05,1);
H = [H,hx];
P = [P,px];
%CI{i} = CIx;
end

degree_nos = {};
rrank = {};
for m = 1:64%length(sub)
    a = all_thres{m};
    s = degrees_und(a);
    [~,p] = sort(s,'descend');
    r = 1:length(s);
    r(p) = r;
    rrank{m} = r;
    degree_nos{m} = s;
end

allh = {};
allnoh = {};
for m = 1:50
h = {};
nh = {};
for i = 1:64%length(sub)
    htem = [];
    nhtem = [];
    b = degree_nos{i};
    for j = 1:90
        d = b(j);
        if (d >= m)
            htem = [htem, j];
        else
            nhtem = [nhtem, j];
        end
    end
    h{i} = htem;
    nh{i} = nhtem;
end
allh{m} = h;
allnoh{m} = nh;
end

save("network_pre_nos.mat", 'm_K', 'k', "label", "all_thres", "all_orin", "all_fix", "all_bin", "all_rand", "all_rnorm", "H", "P", "degree_nos", "allh","allnoh")
end