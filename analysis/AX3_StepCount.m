function [t1,x1,m1,stepcounts,pk_locs]=AX3_StepCount(data,cadence,pk_window,step_abs_thresh,pk_prominence,matdate_start,matdate_stop,filter_style)
% [t1,x1,m1,stepcounts,pk_locs]=AX3_StepCount(data,cadence,pk_window,step_abs_thresh,pk_prominence,matdate_start,matdate_stop,filter_style);
% 
% INPUTS
% data: AX3/AX6 data from AX3_quickdata.m
% cadence: cadence in steps/sec for step peak width (leave empty for default)
% pk_window: epoch window of seconds to check for peaks (leave empty for default)
% step_abs_thresh: peak threshold for counting step (leave empty for default)
% pk_prominence: prominence of step peaks (leave empty for default)
% matdate_start: matlab date/time for beginning of analysis
% matdate_stop: matlab date/time for end of analysis
% filter_style: int 1 for Jan Brond's filters, or 2 for rudimentary filter
% 
% OUTPUTs
% t1: resampled time vector of analysis
% x1: resampled and filtered acceleration signal
% m1: resampled and filtered magnitude signal
% stepcounts: matrix of steps/minute, cols: date/time, cumulative steps,
%   steps in minute, avg cadence in minute, avg peak half-max height
% pk_locs: locations of step peaks on t1/x1/m1
% 
% this functions applies a step counting algorithm to AX3/AX6 acceleration data for 
% time length specified. There are 2 options for filters with default values for peak
% finding. This code is NOT VALIDATED for accurate step counting yet, but
% follows many recommendations from literature. 


%% path to other functions
addpath('..\activitycounts')

%% pre-defined signal parameters
warning('off','signal:findpeaks:largeMinPeakHeight')

%% check inputs
if isempty(cadence)
    cadence = 0.300; %default 0.3 steps per second
end
if isempty(pk_window)
    pk_window = 10; %default 10 seconds
end
if isempty(step_abs_thresh)
    if filter_style==1
        step_abs_thresh = 0.2; %default is 0.3 g
        
    elseif filter_style==2
        step_abs_thresh = 0.3; %default is 0.3 g
        
    else
        step_abs_thresh = 0.3; %default is 0.3 g
    end
end
if isempty(pk_prominence)
    pk_prominence = 0.1; %default 0.1
end

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

%% resample and filter data
[t1, x1, Fs] = filter_accel(t0,x0,filter_style);

% parameters from sampling rate
T = 1/Fs;

%% construct magnitude vector
m1 = sqrt(sum(x1(:,1).^2 + x1(:,2).^2 + x1(:,3).^2,2));
m1 = real(m1);

%% count steps (cumulative by pk_window)
step_width = round(cadence/T);
step_window = round(pk_window/T);

pk_locs = [];

matwindow = datenum(0,0,0,0,0,pk_window);

[dtst_y, dtst_m, dtst_d, dtst_H, dtst_M, dtst_S] = datevec(t1(1));
start0 = datenum(dtst_y, dtst_m, dtst_d, dtst_H, dtst_M+1, 0);

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
        [pks0, locs0] = findpeaks(m1(last0:seg0), 'MINPEAKHEIGHT', step_height0, 'MINPEAKDISTANCE', step_width,'MINPEAKPROMINENCE',pk_prominence);
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
stepcounts = zeros(lim1,4);
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
    
    stepcounts(ix,1) = tmptime1 - matminute;
    stepcounts(ix,2) = tmpseg1(end,2);
    stepcounts(ix,3) = sum(tmpseg1(:,3));
    stepcounts(ix,5) = mean(tmpseg1(:,5));
    
    if(stepcounts(ix,3)>0)
        stepcounts(ix,4) = sum(tmpseg1(:,3).*tmpseg1(:,4))/sum(tmpseg1(:,3));
    else
        stepcounts(ix,4) = 0;
    end

    last1 = seg1+1;
  
    waitbar((ix/lim1),hw,datestr(tmptime1));
    
end

delete(hw);

function [t1, x1, Fs] = filter_accel(t0,x0,filter_style)

switch filter_style
    
    case 1
        %% use Jan Brond's filters from Activity Count code
        
        load('agcoefficients.mat');

        Fs = 30;
        deadband = 0.068;
        peakThreshold = 2.13;
        adcResolution = 0.0164;
        gain = 0.965;
        B = B * gain;
        
        % offset/scale the time vector into units of seconds elapsed 
        time_elapsed = (24*60*60)*[t0 - t0(1)];

        % resample data to 30 Hz
        x1 = resample(x0,time_elapsed,Fs);
        t1 = linspace(t0(1),t0(end),(t0(end)-t0(1))*(24*60*60*Fs));
        t1 = [t1, t1(end) + (1/(24*60*60))/Fs];

        % apply acceleration signal filters
        S = size(x1);
        for n=1:S(2)
            x1(:,n) = filter(B,A,x1(:,n));
        end
        
    case 2
        %% use original rudimentary filters
        
        % get the mean average sampling frequency
        filesf = 1/(mean(diff(t0))*(24*60*60));
        
        Fs = round(filesf);
        
        % Resample data
        t1 = linspace(t0(1),t0(end), (t0(end)-t0(1))*(24*60*60*Fs) )';
        S = size(x0);
        
        accel1 = zeros(length(t1),S(2));
        for n=1:S(2)
            accel1(:,n) = interp1(t0,x0(:,n),t1,'pchip',0);
        end

        L1 = size(t1,1);
        
        %% create low pass frequency filter: ramp up to 1, ramp down to 0.1, and zero after Hz>30

        LPc = 2; % Hz low pass to 1 (ramp down)
        LPc2 = 5; % Hz low pass to 0.1
        LPc0 = 0.13; % Hz high pass to 1 (ramp up)

        %frequency axis SIGNAL 1
        L = L1;
        YHf = (Fs/L)*[0:((L/2)-1)];

        %ramp/zero cutoffs
        fc = find(YHf>=LPc,1);
        fc2 = find(YHf>=LPc2,1);
        fc0 = find(YHf>=LPc0,1);
        fc30 = find(YHf>=30,1);

        %start with a baseline filter of 0.1 at all frequencies
        r = 0.1*ones(L,1);

        %low frequency pass
        r(1:(fc+1)) = 1;
        r(L-(fc+1):L) = 1;

        %ramp down between LPc to LPc2, this reduces Gibbs phenomenon
        r((fc+1):fc2) = linspace(1,0.1,length([(fc+1):fc2]));
        r((L-fc2):(L-(fc+1))) = linspace(0.1,1,length((L-fc2):(L-(fc+1))));

        %ramp up between 0 Hz and LPc0, this removes offsets
        r(1:fc0) = linspace(0,0.5,length([1:fc0]));
        r(L-fc0:L) = linspace(0.5,0,length([L-fc0:L]));

        %remove high frequencies
        r(fc30:(L-fc30)) = 0;

        YHf1 = YHf;
        r1 = r;
        
        %% filter the signals
        x1 = [];

        for n = 1:S(2)
            % accel signal
            tmp_accel = accel1(:,n);
            % fourier transform to frequency domain
            Y = fft(tmp_accel);
            % filter FT and inverse-fourier transform back to time domain
            x1(:,n) = ifft(Y.*r1);
        end

end


