function [theMatrix]=demeaner(varargin)
%%% first input is inputMatrix, second is dimention
%%%% if dimention is 1 or nothing it takes the average of all rows and
%%%% subtracts that to each row, if it is 2 it does the same with the
%%%% columns
inputMat=varargin{1};
if nargin==1
    meanMat=mean(inputMat,1);
    theMatrix=inputMat-(meanMat'*ones(1,size(inputMat,1)))';%takes the realization matrix and subtracts the mean to each line
elseif varargin{2}==1
    meanMat=mean(inputMat,1);
    theMatrix=inputMat-(meanMat'*ones(1,size(inputMat,1)))';%takes the realization matrix and subtracts the mean to each line
else
inputMat=inputMat';
meanMat=mean(inputMat,1);
theMatrix=inputMat-(meanMat'*ones(1,size(inputMat,1)))';%takes the realization matrix and subtracts the mean to each line
theMatrix=theMatrix';
end
end