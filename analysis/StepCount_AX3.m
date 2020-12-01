function [t1,S_lp1,m1,steps1,pk_locs]=StepCount_AX3(data,cadence,pk_window,matdate_start,matdate_stop)

%% path to other functions
addpath('..\data_io')
addpath('..\activitycounts')

%% pre-defined signal parameters
Fs = 30;
T = 1/Fs;

%% read in accel data and convert data to units of g
x0 = double([data.x, data.y, data.z]);
x0 = x0*data.AccScale;

%% read in time data
[t0] = AX3_interpolatetime(data);

%% Segment and sync the data
a1 = find(t0>=matdate_start,1,'first');
b1 = find(t0<=matdate_stop,1,'last');

if(isempty(a1))
    a1 = 1;
end
if(isempty(b1))
    b1 = size(t0,1);
end

x0 = x0(a1:b1,:);
t0 = t0(a1:b1);

%% resample and filter data with J. Brond's filters
load('agcoefficients.mat');

deadband = 0.068;
peakThreshold = 2.13;
adcResolution = 0.0164;
gain = 0.965;

B = B * gain;

% get the mean average sampling frequency
filesf = 1/(mean(diff(t0))*(24*60*60));

% offset/scale the time vector into units of seconds elapsed 
time_elapsed = (24*60*60)*[t0 - t0(1)];

% resample data to 30 Hz
x1 = resample(x0,time_elapsed,Fs);
t1 = linspace(t0(1),t0(end),(t0(end)-t0(1))*(24*60*60*Fs));
t1 = [t1, t1(end) + (1/(24*60*60))/Fs];


%% construct magnitude vector2
m1 = sqrt(sum(x1(:,1).^2 + x1(:,2).^2 + x1(:,3).^2,2));
m1 = real(m1);

%% count steps (cumulative by pk_window)
step_abs_thresh = 0.3;
step_width = round(cadence/T);
step_window = round(pk_window/T);

pk_locs = [];

matwindow = datenum(0,0,0,0,0,pk_window);

[dtst_y, dtst_m, dtst_d, dtst_H, dtst_M, dtst_S] = datevec(t1(1));
start0 = datenum(dtst_y, dtst_m, dtst_d, dtst_H, dtst_M+1, 0);
% start0 = datenum(t1(1));
lim0 = floor((datenum(t1(end)) - start0)/(matwindow));
steps0 = zeros(lim0,5);
last0 = 1;
lastpks0 = 0;

hw = waitbar(0,'Counting Steps...');
for ix=1:lim0

    tmptime0 = (start0 + (ix*matwindow));
    
    seg0 = find(t1>(tmptime0),1);
    if(isempty(seg0))
        seg0 = length(t1);
    end
    
    window0 = (seg0-step_window);
    if(window0<1)
        window0 = 1;
    end
    tmpseg0 = m1(window0:seg0);
    
    step_height0 = (max(tmpseg0)-min(tmpseg0))/2 + min(tmpseg0);

    if(step_height0<step_abs_thresh)
        step_height0 = step_abs_thresh;
    end
    
    steps0(ix,1) = tmptime0;
    steps0(ix,5) = step_height0;
    
    if(max(m1(last0:seg0))<=step_height0)
        steps0(ix,2) = lastpks0;
        steps0(ix,3) = 0;
        steps0(ix,4) = 0;
    else
%         [pks0, locs0] = findpeaks(m1(last0:seg0), 'MINPEAKHEIGHT', step_height0, 'MINPEAKDISTANCE', step_width);
        [pks0, locs0] = findpeaks(m1(last0:seg0), 'MINPEAKHEIGHT', step_height0, 'MINPEAKDISTANCE', step_width,'MINPEAKPROMINENCE',0.1);
        tdiff0 = diff(locs0)*T;
        steps0(ix,2) = lastpks0 + size(pks0,1);
        steps0(ix,3) = size(pks0,1);
        if(size(locs0,1)>1)
            steps0(ix,4) = mean(tdiff0);
        else
            steps0(ix,4) = 0;
        end
        
        pk_locs = [pk_locs; ((locs0-1)+last0)];
        
    end
    
    last0 = seg0+1;
    lastpks0 = steps0(ix,2);
    
    waitbar((ix/lim0),hw,datestr(tmptime0));
    
end
delete(hw);

%% collect steps by minute
matminute = datenum(0,0,0,0,1,0);
start1 = datenum(t1(1));
lim1 = floor((datenum(t1(end)) - datenum(t1(1)))/(matminute));
steps1 = zeros(lim1,4);
last1 = 1;
lastpks1 = 0;

hw = waitbar(0,'Counting Steps...');

for ix=1:lim1

    tmptime1 = (start1 + (ix*matminute));
    
    seg1 = find(steps0(:,1)>=tmptime1,1,'first');
    if(isempty(seg1))
        seg1 = size(steps0,1);
    end
    
    tmpseg1 = steps0(last1:seg1,:);
    
    steps1(ix,1) = tmptime1 - matminute;
    steps1(ix,2) = tmpseg1(end,2);
    steps1(ix,3) = sum(tmpseg1(:,3));
    steps1(ix,5) = mean(tmpseg1(:,5));
    
    if(steps1(ix,3)>0)
        steps1(ix,4) = sum(tmpseg1(:,3).*tmpseg1(:,4))/sum(tmpseg1(:,3));
    else
        steps1(ix,4) = 0;
    end

    last1 = seg1+1;
  
    waitbar((ix/lim1),hw,datestr(tmptime1));
    
end

delete(hw);





