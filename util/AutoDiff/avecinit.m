function x = avecinit(v)

%AVECINIT Initialize audi vector.
%   x = AVECINIT(v) is a column vector of n=numel(v) independent audi
%   variables. It is a shortcut for [x1,x2,...,xn] = ainit(v1,v2,...,vn).
%
%   CAUTION: Subsequent computations are likely to be slow if n is large.
%
%   See also: ainit

for i = numel(v):-1:1
  x(i,1) = audi(v(i),i,numel(v),1);
end