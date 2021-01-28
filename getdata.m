% Whole solid image
X1 = csvread('WSI_150_261.csv');
X1=X1';
X1=log2(X1+1);
% DNA methylation
X2 = csvread('meth_18510_261.csv');
X2=X2';
X2=log2(X2+1);
% Copy number variation
X3 = csvread('CNV_5408_261.csv');
X3=X3';
X3=log2(X3+1);
%Correlation WSI&Meth MSI&CNV
A = csvread('meth_WSI_18510_150.csv');
B = csvread('CNV_WSI_5408_150.csv');