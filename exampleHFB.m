
fn = 'R:\MSS\Johnson_Lab\DATA\NMH\NM05\Tasks\Rest\BPR\NM05_rest.mat';
data = load(fn).data;



dat = cat(3, data.trial{:});

highfrex = linspace(70, 150, 81); 
tim = -600:2399;

allpow = arrayfun(@(x) tempFunc(squeeze(dat(x,:,:)),highfrex,...
                        1000, tim, 5), [1:2], ...   %edit here to increase from just channels 1:2 and instead do all channels
    'uniformoutput', false);
timOut = [-600:5:2399];



