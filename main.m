%% load the real input data
getdata;
%% preprocess the real input data
A = sparse(abs(A')); %WSI&Meth
B = sparse(abs(B')); %MSI&CNV
%% SVD iniialization
X = [X1 X2 X3];
K = 15;
[U,S,V]=svd(X);
s_sum = sum(S); 
W(:,1) = sqrt(S(1,1))*U(:,1);
H(1,:) = sqrt(S(1,1))*V(:,1)';
for j = 2 : K
    x = U(:,j);
    y = V(:,j);
    xp = (x>=0).*x;
    xn = (x<0).*(-x);
    yp = (y>=0).*y;
    yn = (y<0).*(-y);
    xpnrm = norm(xp);
    ypnrm = norm(yp);
    mp = xpnrm * ypnrm;
    xnnrm = norm(xn);
    ynnrm = norm(yn);
    mn = xnnrm * ynnrm;
    if mp > mn
        u = xp/xpnrm;
        v = yp/ypnrm;
        sigma = mp;
    else
        u = xn/xnnrm;
        v = yn/ynnrm;
        sigma = mn;
    end
    W(:,j) = sqrt(S(j,j)*sigma)*u;
    H(j,:) = sqrt(S(j,j)*sigma)*v';
end
W = abs(W);
save W_original.mat W;
H = abs(H);
H1 = H(:,1:150);
H2 = H(:,151:18660);
H3 = H(:,18661:24068);
save H1_original.mat H1;
save H2_original.mat H2;
save H3_original.mat H3;
%% applying MDJNMF
L1 = 0.001; L2 = 0.01; r1 = 1; r2 =1; K = 15; a = 0.001;
    tic
    [W,H1,H2,H3] = MCJNMF_comodule(X1,X2,X3,A,B,a,r1,r2,L1,L2,K);
    toc
 %% get comodules(index)
tt = 2;
[Co_module] =Comodule_selection(W, H1, H2, H3, tt);
save Comodule_tt_2.mat Co_module;
%% module elements extraction
B1=zeros(15,1500);
B2=zeros(15,1500);
B3=zeros(15,1500);
B1(B1==0)=[];
B2(B2==0)=[];
B3(B3==0)=[];
for i=1:15
    A1=Co_module(i,1);
    A2=Co_module(i,2);
    A3=Co_module(i,3);
    A11=cell2mat(A1);
    A22=cell2mat(A2);
    A33=cell2mat(A3);
    k1=length(A11);
    k2=length(A22);
    k3=length(A33);
    for j1=1:k1
        B1(i,j1)=A11(1,j1);
    end
    for j2=1:k2
        B2(i,j2)=A22(1,j2);
    end
    for j3=1:k3
        B3(i,j3)=A33(1,j3);
    end
end
B1(B1==0)=NaN;
B2(B2==0)=NaN;
B3(B3==0)=NaN;

xlswrite('Co_module_tt_2_K_15_svd_WSI_2.xlsx',B1');
xlswrite('Co_module_tt_2_K_15_svd_METH_2.xlsx',B2');
xlswrite('Co_module_tt_2_K_15_svd_CNV_2.xlsx',B3');