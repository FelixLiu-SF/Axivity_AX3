function counts = agfilt_felix(data,t,integN_10hz,B,A)
% counts = agfilt_felix(data,t,integN_10hz,B,A);
% 
% data is AX3 data vector in g units
% t is AX3 time vector in matlab time units with same length as data
% integN_10hz is the summing window length for activity counts


% Actigraph characteristics
deadband = 0.068;
sf = 30;
peakThreshold = 2.13;
adcResolution = 0.0164;
gain = 0.965;

% get the mean average sampling frequency
filesf = 1/(mean(diff(t))*(24*60*60));

%% offset/scale the time vector into units of seconds elapsed 
t2 = (24*60*60)*[t - t(1)];

if (filesf>sf)
    %resample data down to 30Hz
    dataf = resample(data,t2,sf);
else
    %Aliasing Filter
    [B2,A2] = butter(4,[0.01 7]./(sf/2));
    dataf = filtfilt(B2,A2,data);
end

S = size(dataf);

B = B * gain;

for n=1:S(2)
    
    %apply Actigraph filter
    fx8up = filter(B,A,dataf(:,n));
    
    %downsample to 10 Hz
    fx8up = resample(fx8up,10,sf);
    
    %truncate max signal magnitude
    fx8 = pptrunc(fx8up,peakThreshold);
    
    %sum the activity counts
    counts(:,n) = runsum(floor(trunc(abs(fx8),deadband)./adcResolution),integN_10hz,0);
    
end
    