function [wtv]=AX3_weartime(data,epoch_m,bigsearch_n,epoch_min,active_threshold)
% [wtv]=AX3_weartime(data,epoch_m,bigsearch_n);
% 
% INPUTS
% data: AX3/AX6 data from AX3_quickdata.m
% epoch_m: epoch length in minutes for output (1 minute to 1 day range)
% bigsearch_n: number of epochs to search at once. larger n is faster. 
% epoch_min: minimum size of epoch to check for wear time, leave empty for
%   default, 30 seconds.
% active_threshold: threshold of mininum epoch seconds per epoch_m for epoch to be
%   classified as wear-time. Leave empty for default, 60 seconds. 
% 
% OUTPUTS
% wtv: output cell with wear time classification by epoch, with cols
%   1) date/time of epoch start
%   2) [0:nonwear-time detected, 1:wear-time detected]
%   3) standard deviation of x,y,z axes in epoch
%   4) range of x,y,z axes in epoch
%   5) active period: sum of consecutive mininum epochs with wear-time, in seconds
%   6) binary indicator vector of wear time by minimum epoch
% 
% Calculate AX3 monitor wear time. The wear time is originally based on
% algorithms from Axivity OpenMovement, but modified for speed and is also 
% slightly heuristic; it takes into account previous and future epochs 
% when determining wear time. 

% parse inputs
if(isnumeric(epoch_m))
    if(epoch_m<=1) %invalid epoch time, default to 1 minute
        epoch_m = 1; 
        
    elseif(epoch_m>=(60*24)) %invalid epoch time, default to 1 day
        epoch_m=(60*24);
    end
else %invalid epoch time, default to 1 minute
    epoch_m = 1;
end

if(bigsearch_n<1)
    bigsearch_n = 1;
end

if(isempty(epoch_min))
    epoch_min = 30;
end

if(isempty(active_threshold))
    active_threshold = 60;
end

% set epoch parameters
s_min = epoch_min; % minimum epoch to check, in seconds

% get epoch parameters
s_total = epoch_m*60; %epoch in seconds
n_wtv = ceil(s_total/s_min); % number of minimum chunks in output epoch        


%declare variables
wtv = []; %wear-time output
epoch_std = zeros(1,3); %array for epoch st. devs.
epoch_r = zeros(1,3); %array of epoch ranges

% interpolate time data
[t] = AX3_interpolatetime(data);

%get start & stop times
[yr1,mo1,day1,hr1,mn1,~] = datevec(t(1));

[yr2,mo2,day2,hr2,mn2,~] = datevec(t(end));

%get times in number of epochs
min_to_start = (hr1*60) + mn1;
min_to_stop = (hr2*60) + mn2;

epoch_before_start = floor(min_to_start/epoch_m);
epoch_after_stop = ceil(min_to_stop/epoch_m);

start_time = datenum(yr1,mo1,day1,0,epoch_before_start*epoch_m,0);
stop_time = datenum(yr2,mo2,day2,0,epoch_after_stop*epoch_m,0);

recorded_days = stop_time - start_time;
recorded_min = recorded_days*24*60;
last_epoch = round(recorded_min/epoch_m);

%calculate wear-time in bug chunks and recursively for any detections
bx=1;
while(bx<=last_epoch)
    
    % big search
    epoch1 = (bx-1)*datenum(0,0,0,0,epoch_m,0) + start_time;
    epoch2 = (bx+bigsearch_n-1)*datenum(0,0,0,0,epoch_m,0) + start_time;
    
    % get indices of big epoch
    t1 = find(t>=epoch1,1,'first');
    t2 = find(t<epoch2,1,'last');
    
    % sanitize indices
    t1 = min([t1,length(t)]);
    t2 = max([t2,1]);

    % get acceleration data of big epoch
    if(t2-t1>1)
        x1 = double(data.x(t1:t2))*data.AccScale;
        y1 = double(data.y(t1:t2))*data.AccScale;
        z1 = double(data.z(t1:t2))*data.AccScale;
    else
        x1=0;
        y1=0;
        z1=0;
    end

    % get standard deviation and range of epoch
    epoch_std(1,1) = std(x1);
    epoch_std(1,2) = std(y1);
    epoch_std(1,3) = std(z1);

    epoch_r(1,1) = max(x1) - min(x1);
    epoch_r(1,2) = max(y1) - min(y1);
    epoch_r(1,3) = max(z1) - min(z1);

    % check st dev and range for wear time
    n_std = size(find(epoch_std<(3/1000)),2);   %number of axes below std threshold
    n_r = size(find(epoch_r<(50/1000)),2);      %number of axes below range threshold
    
    if(n_std>1 || n_r>1 || sum(isnan(epoch_std)>1) || sum(isnan(epoch_r))>1) 
        % too many axes below threshold -> nonwear
        
        % set all epochs in big epoch range to nonwear
        epoch_cell = [epoch1:datenum(0,0,0,0,epoch_m,0):(epoch2-datenum(0,0,0,0,epoch_m,0))]';
        epoch_cell = mat2cell(epoch_cell,ones(bigsearch_n,1),[1]);
        
        wtv = [wtv; [epoch_cell, ...
                repcell([bigsearch_n,1],0), ...
                repcell([bigsearch_n,1],epoch_std), ...
                repcell([bigsearch_n,1],epoch_r), ...
                repcell([bigsearch_n,1],0), ...
                repcell([bigsearch_n,1],0)]];
        
        bx=bx+bigsearch_n;
        
    else
        % wear-time detected in big epoch, check individual epochs for wear-time
        
        for ix=bx:(bx+bigsearch_n-1)
            
            % search a single epoch
            epoch1 = (ix-1)*datenum(0,0,0,0,epoch_m,0) + start_time;
            epoch2 = (ix)*datenum(0,0,0,0,epoch_m,0) + start_time;

            % get indices of one epoch
            t1 = find(t>=epoch1,1,'first');
            t2 = find(t<epoch2,1,'last');

            % sanitize indices
            t1 = min([t1,length(t)]);
            t2 = max([t2,1]);

            % get acceleration data of epoch
            if(t2-t1>1)
                x1 = double(data.x(t1:t2))*data.AccScale;
                y1 = double(data.y(t1:t2))*data.AccScale;
                z1 = double(data.z(t1:t2))*data.AccScale;
            else
                x1=0;
                y1=0;
                z1=0;
            end

            % get standard deviation and range of epoch
            epoch_std(1,1) = std(x1);
            epoch_std(1,2) = std(y1);
            epoch_std(1,3) = std(z1);

            epoch_r(1,1) = max(x1) - min(x1);
            epoch_r(1,2) = max(y1) - min(y1);
            epoch_r(1,3) = max(z1) - min(z1);

            % check st dev and range for wear time
            n_std = size(find(epoch_std<(3/1000)),2);   %number of axes below std threshold
            n_r = size(find(epoch_r<(50/1000)),2);      %number of axes below range threshold

            if(n_std>1 || n_r>1 || sum(isnan(epoch_std)>1) || sum(isnan(epoch_r))>1) 
                % too many axes below threshold -> nonwear
                
                wtv = [wtv; {epoch1, 0, epoch_std, epoch_r, 0, 0}];

            else % check in 30 second periods for how much activity detected

                wtv_periods = zeros(n_wtv,1);
                wtv_periods = AX3_recursive_weartime(wtv_periods,[x1,y1,z1],t(t1:t2),s_total,1,s_total,s_total/2,s_min);

                
                %% total active periods of wear time in minimum epochs
                
                % active_periods = min([size(find(wtv_periods==1),1)*s_min, s_total]); %deprecated
        
                % get total wear time from consecutive periods only
                
                % get indices at which wear time binary switches 
                ap_indc = [1; find(diff(wtv_periods)) + 1];
                % get lengths of consecutive epochs of wear or non-wear
                ap_len = diff([ap_indc; length(wtv_periods) + 1]);
                % multiple consecutive length with wear-time binary
                ap_wlen = wtv_periods(ap_indc).*ap_len;
                
                % get sum of consecutive wear time epochs only
                active_periods = min([sum(ap_wlen(ap_wlen>1))*s_min,s_total]);

                if(active_periods>=active_threshold)
                    wtv = [wtv; {epoch1, 1, epoch_std, epoch_r, active_periods, wtv_periods}];
                else
                    wtv = [wtv; {epoch1, 0, epoch_std, epoch_r, active_periods, wtv_periods}];
                end
                
            end
        end
        bx = ix+1;
    end
end

% Close "holes" if output epochs are 30 minutes or less. Doesn't make sense
% to close holes if epochs are much larger than that, in my opinion. 
if(epoch_m<=30)
    
    % close 'holes' on non-wear between wear-time if 2 or less epochs wide
    wtv_cell = cellfun(@num2str,wtv(:,2),'UniformOutput',0)';
    wtv_str = cellfun(@horzcat,wtv_cell);

    holes = findstr(wtv_str,'1001');
    while(~isempty(holes))
        wtv_str(holes(1):(holes(1)+3)) = '1111';
        holes = findstr(wtv_str,'1001');
    end

    holes = findstr(wtv_str,'101');
    while(~isempty(holes))
        wtv_str(holes(1):(holes(1)+2)) = '111';
        holes = findstr(wtv_str,'101');
    end

    fill_str = mat2cell(wtv_str',ones(size(wtv(:,2))),[1]);
    fill_cell = cellfun(@str2num,fill_str,'UniformOutput',0);

    wtv(:,2) = fill_cell;


    % close 'holes' on wear-time between non-wear,
    % if 1-2 epoch wide and isolated (no wear-time nearby)
    wtv_cell = cellfun(@num2str,wtv(:,2),'UniformOutput',0)';
    wtv_str = cellfun(@horzcat,wtv_cell);

    holes = findstr(wtv_str,'0001000');
    while(~isempty(holes))
        wtv_str(holes(1):(holes(1)+6)) = '0000000';
        holes = findstr(wtv_str,'0001000');
    end

    holes = findstr(wtv_str,'00011000');
    while(~isempty(holes))
        wtv_str(holes(1):(holes(1)+7)) = '00000000';
        holes = findstr(wtv_str,'00011000');
    end

    fill_str = mat2cell(wtv_str',ones(size(wtv(:,2))),[1]);
    fill_cell = cellfun(@str2num,fill_str,'UniformOutput',0);

    wtv(:,2) = fill_cell;

end

function wtv_periods = AX3_recursive_weartime(wtv_periods,period_in,t_in,s_total,wtv_ix,s_in,s_now,s_min)

% round s_now to nearest factor of s_min seconds
s_double = round((s_now*2)/s_min)*s_min;
s_now = round(s_now/s_min)*s_min;
if((s_double-s_now)<s_min)
    s_now = s_now - s_min;
end
if(s_now<s_min)
    s_now = s_min;
end

% calculate indices
last_period = ceil(s_in/s_now);
n_wtv = size(wtv_periods,1);

for jx=1:last_period
    
    ix_1 = ((jx-1)*s_now)/s_min + wtv_ix;
    ix_2 = (jx*s_now)/s_min  + wtv_ix -1;
    if(ix_2>n_wtv)
        ix_2 = n_wtv;
    end
    
    try
    period1 = (jx-1)*datenum(0,0,0,0,0,s_now) + t_in(1);
    period2 = (jx)*datenum(0,0,0,0,0,s_now) + t_in(1);

    p1 = find(t_in>=period1,1,'first');
    p2 = find(t_in<period2,1,'last');
    if(isempty(p1))
        p1=1;
    end
    if(isempty(p2))
        p2=size(t_in,2);
    end

    x1 = period_in(p1:p2,1);
    y1 = period_in(p1:p2,2);
    z1 = period_in(p1:p2,3);
    
    period_std(1,1) = std(x1);
    period_std(1,2) = std(y1);
    period_std(1,3) = std(z1);

    period_r(1,1) = max(x1) - min(x1);
    period_r(1,2) = max(y1) - min(y1);
    period_r(1,3) = max(z1) - min(z1);

    p_std = size(find(period_std<(3/1000)),2);  %number of axes below std threshold
    p_r = size(find(period_r<(50/1000)),2);     %number of axes below range threshold
    catch err
        disp(err.message);
        p_std = 0;
        p_r = 0;
    end
    
    if(p_std>1 || p_r>1) 
        %too many axes below threshold -> whole period is nonwear
        wtv_periods(ix_1:ix_2,1) = 0;
        
    elseif(s_now<=s_min)
        %reached minimal period size, count this period as wear time
        wtv_periods(ix_1:ix_2,1) = 1; 

    else
        %call recursive function to further subdivide this period
        wtv_periods = AX3_recursive_weartime(wtv_periods,[x1,y1,z1],t_in(p1:p2),s_total,ix_1,s_now,s_now/2,s_min);
        
    end
end