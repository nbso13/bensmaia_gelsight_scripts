function s = stim_sinesweep(freq,amp,len,loc,samp_freq,ramp_len)
% s = stim_sinesweep(freq,amp,len,loc,samp_freq,ramp_len)

assert(length(freq)==2);
assert(length(amp)==2);

if nargin<6 || isempty(ramp_len)
    ramp_len = 0.05;
end

if nargin<5 || isempty(samp_freq)
    samp_freq = 5000;
end

if nargin<4 || isempty(loc)
    loc = [0 0];
end

if nargin<3 || isempty(len)
    len = 1;
end

ref = sin(linspace(0,2*pi*1000,samp_freq*1000));

ampsweep = linspace(amp(1),amp(2),samp_freq*len);
freqsweep = linspace(freq(1),freq(2),samp_freq*len);
trace = ampsweep.*ref(cumsum(round(freqsweep)));
trace = apply_ramp(trace,ramp_len,samp_freq);

s = Stimulus(reshape(trace,[],1),loc,samp_freq);
