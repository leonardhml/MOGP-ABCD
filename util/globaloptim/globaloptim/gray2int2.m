function int = gray2int2(gray)
%GRAY2INT2 Convert a gray code array to an integer.
%   INT = GRAY2INT2(GRAY) Returns GRAY's decimal integer value.
%   Sometimes people use gray codes to interpret their binary strings because 
%   gray codes have the property that adjacent numerical values differ by only 
%   one bit, whereas straight binary can differ by many. for example 7 is four
%   bits away from 8 in binary.

%   Copyright 2003-2004 The MathWorks, Inc.

n = length(gray);

bin = zeros(1,n);
bin(n) = gray(n);
for i = (n-1):-1:1
    bin(i) = xor(bin(i+1),gray(i));
end

map = 2 .^ (0:n);
int = abs(sum(map(find(bin))));