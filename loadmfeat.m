function [X , W , ind , label , viewNum , latentDim] = loadmfeat(percentDel)
load('mfeat.mat')
viewNum = 5;
latentDim = 10;
X = cell(viewNum,1);
W = cell(viewNum,1);
label = [];
for i = 1:latentDim
   label = [label;ones(200,1)*i];
end
X{1} = mfeat_fac';
X{2} = mfeat_fou';
X{3} = mfeat_kar';
X{4} = mfeat_pix';
X{5} = mfeat_zer';

ind = ones(size(X{1},1), viewNum);
ind = splitDigitData( ind, percentDel, 1 );
for i = 1:viewNum
    item = ind(:,i);
    W{i} = diag(item);
    temp = find(item == 0); 
    X{i}(temp,:)= 0;
    X{i} = mapminmax(X{i},0,1);
    X{i} = NormalizeFea(X{i},1);
end
end