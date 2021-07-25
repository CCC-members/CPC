clc;clear;close all;
%% load IRIS data
load iris_dataset;
%% set data,n,lmax,pmax
S1 = irisInputs(:,1:50)';
S2 = irisInputs(:,51:100)';
S3 = irisInputs(:,101:150)';

X(:,:,1)=S1;
X(:,:,2)=S2;
X(:,:,3)=S3;

n = [50 50 50];
lmax = 50;
pmax = 4;
%% compute stepwiseCFCPC
[Lambda_f,Q_f] = stepwiseCFCPC(X,n,pmax,lmax);