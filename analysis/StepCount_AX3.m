function [t1,S_lp1,m1,steps1,pk_locs]=StepCount_AX3(data,cadence,pk_window,matdate_start,matdate_stop)


%% pre-defined signal parameters
Fs = 100;
T = 1/Fs;

%% read in date files
% x1 = dlmread(csv_in);
x0 = [data.x, data.y, data.z];
x0 = double(x0);
x0 = x0/256;

[t0] = AX3_interpolatetime(data);
t0 = t0';

%% Sync the data
% sync_start = matdate_start;
% sync_stop = matdate_stop;

% t1 = x1(1:end,1);
a1 = find(t0>=matdate_start,1,'first');
b1 = find(t0<=matdate_stop,1,'last');
if(isempty(a1))
    a1 = 1;
end
if(isempty(b1))
    b1 = size(t0,1);
end

x1 = x0(a1:b1,:);
a2 = t0(a1);
b2 = t0(b1);

% t1 = t0(a1:b1,1);
% accel1 = x1;

% Resample data
t1 = linspace(a2,b2, (b2-a2)*(24*60*60*Fs) )';
accel1 = zeros(length(t1),3);
for ax=1:3
    accel1(:,ax) = interp1(t0(a1:b1,1),x1(:,ax),t1,'pchip',0);
end

L1 = size(t1,1);

clear t0 x1;

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
clear YHf r L x1;

%% filter the signals
S_lp1 = [];

for dx = 1:3
    % accel signal 1
    S = accel1(:,dx);
    % fourier transform to frequency domain
    Y = fft(S);
    % filter FT and inverse-fourier transform back to time domain
    S_lp1(:,dx) = ifft(Y.*r1);
end


%% construct magnitude vector2
m1 = sqrt(sum(S_lp1(:,1).^2 + S_lp1(:,2).^2 + S_lp1(:,3).^2,2));
m1 = real(m1);

clear accel1;

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





