function [Co_module] = Comodule_selection(W, H1, H2, H3, tt)
%
% INPUT
% W         : common basis matrix
% H1        : coefficient matrix
% H2        : coefficient matrix
% H3        : coefficient matrix
% tt        : a given threshold for z-score.
% OUTPUT
% Co_module : the index list of WSI, DNA methylation, CNV
%
% Compute the mean(meadia) and std in columns of W and rows in H1, H2, H3 to determine
% the module member and output the Co-module based on W and H1, H2, H3
% matrices.
%
m1 = size(H1,2);
m2 = size(H2,2);
m3 = size(H3,2);
n = size(W,1);
K = size(W,2);

MW = mean(W,1);     MH1 = mean(H1,2);   MH2 = mean(H2,2);  MH3 = mean(H3,2);
VW = std(W,0,1);    VH1 = std(H1,0,2);  VH2 = std(H2,0,2); VH3 = std(H3,0,2);

% Co-Module
    for i = 1:K
        c1 = find(H1(i,:) > MH1(i) + tt*VH1(i));
        c2 = find(H2(i,:) > MH2(i) + tt*VH2(i));
        c3 = find(H3(i,:) > MH3(i) + tt*VH3(i));
        Co_module{i,1}=c1; Co_module{i,2}=c2; Co_module{i,3}=c3;
    end