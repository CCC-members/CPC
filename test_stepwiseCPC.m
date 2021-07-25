clc;clear;close all;
%% load IRIS data
load iris_dataset;
%% set data,n,lmax,pmax
S1 = irisInputs(:,1:50)';
S2 = irisInputs(:,51:100)';
S3 = irisInputs(:,101:150)';
n = [50 50 50];
lmax = 50;
pmax = 4;
%% Compute Covariance matrices
S1 = S1-mean(S1);
S2 = S2-mean(S2);
S3 = S3-mean(S3);
S(:,:,1)=(S1'*S1)./n(1);
S(:,:,2)=(S2'*S2)./n(2);
S(:,:,3)=(S3'*S3)./n(3);
%% compute stepwiseCPC
[Lambda,Q,qtilde] = stepwiseCPC(S,n,pmax,lmax);