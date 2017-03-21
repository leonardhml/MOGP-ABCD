disp('hello world');

disp('testing parallel computation');

cluster = parcluster('local');
numworkers = cluster.NumWorkers;
parpool(numworkers);

res = {};
for i = 1:20
    pause(5);
    res(i).val = i;
    disp(i);
end

save('res.mat', 'res');
