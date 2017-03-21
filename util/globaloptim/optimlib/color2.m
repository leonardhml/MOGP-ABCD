function [group] = color2(J,p)
%

%color2	Column partition for sparse finite differences.
%
%	 GROUP = color2(J,P) returns a partition of the 
%   column corresponding to a coloring of the column-
%   intersection graph. GROUP(I) = J means column I is 
%   colored J.
%
%	 All columns belonging to a color2 can be estimated 
%   in a single finite difference. 
%

%   Copyright 1990-2006 The MathWorks, Inc.

%
[m,n] = size(J);
if nargin < 2, 
   p = 1:n; 
end
J = J(:,p);
group = zeros(n,1);
ncol = 0; 
J = spones(J);
while any(group==0)  
   % Build group for ncol
   ncol = ncol + 1;
   rows = zeros(m,1);
   index = find(group == 0);
   lenindex = length(index);
   for i = 1:lenindex
      k = index(i);
      inner = J(:,k)'*rows;
      if inner == 0
         group(k) = ncol;
         rows = rows + J(:,k);
      end
   end
end
group(p)= group;
