function [ output ] = preprocessTrafficData(csv, r, c, resolution)
% r: row number to start with (0-indexed)
% c: col number to start with (0-indexed)
% resolution: Set how finely detailed the output will be. If 1, resolution
% will be 5 minutes. If 6, 30 minutes, and so on.
raw = csvread(csv, r, c);
[l, w] = size(raw);
output = [];
for i = 1:resolution:length(raw)
    day = raw(i,1);
    mon = raw(i,2);
    hh = raw(i,3);
    mm = raw(i,4);
    ss = raw(i,5);
    time = (day - 1) * 1440 + hh * 60 + mm;
    spd = mean(raw(i:i+resolution-1, 6));
    vol = mean(raw(i:i+resolution-1, 7));
    occ = mean(raw(i:i+resolution-1, 8));
    output = [output; time spd vol occ];
end

% Normalise output

normcols = normc(output(:, 2:end));
output = [output(:,1) normcols];

end

