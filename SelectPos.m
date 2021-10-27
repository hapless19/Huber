function [posSample, p] = SelectPos(I,tmplsize, param, posNum,posSig)
posSample = zeros(tmplsize(1), tmplsize(2), posNum);
p = zeros(6, posNum);
posSample(:,:,posNum) = warpimg(I, param, tmplsize);
p(:,posNum) = param; 

param = repmat(affparam2geom(param(:)), [1,posNum-1]);
param = param + randn(6,posNum-1).*repmat(posSig(:),[1,posNum-1]); 
wimgs = warpimg(I, affparam2mat(param), tmplsize);
posSample(:,:, 1:posNum-1) = wimgs;
p(:,1:posNum-1) = affparam2mat(param);


