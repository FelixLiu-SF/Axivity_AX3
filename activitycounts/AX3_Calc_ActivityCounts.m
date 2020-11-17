function [ActivityCounts, t_start, metadata] = AX3_Calc_ActivityCounts(cwa_file,AC_filterA,AC_filterB)
% This function calculates Actigraph-equivalent Activity Counts from AX3
% using the filters from Jan Brond
% 
% returns Activity Counts per minute in x,y,z axes


% Acitivty Count filter parameters
hz_count = 10;
sec_buffer = 60;
sec_drift = 0;

% load the AX3 cwa file
[data,metadata]=AX3_quickdata(cwa_file);
[t]=AX3_interpolatetime(data);

% get time vectors for the AX3 recording
[Y1,M1,D1,hr1,min1,sec1] = datevec(t(1));
[Y2,M2,D2,hr2,min2,sec2] = datevec(t(end));

% throw away at least 1 minute from beginning and end
t_start = datenum(Y1,M1,D1,hr1,min1+3,0);
t_end = datenum(Y2,M2,D2,hr2,min2-3,0);

% AX3 buffer and drift times for fine-tuning
ax3_buffer = datenum(0,0,0,0,0,sec_buffer);
ax3_drift = datenum(0,0,0,0,0,sec_drift);

% get length of AX3 recording in terms of minutes
min_length = round((t_end - t_start)*(24*60));

% get length of AX3 recording in terms of hours
hours_length = floor(min_length/60);
hours_modulus = mod(min_length,60);

min_rem_start = min_length - (hours_modulus-1);

% calculate AX3 activity counts for every minute
ActivityCounts = zeros(min_length,3);
mx=1;

for hx=1:hours_length
    
    mx1 = (hx-1)*60 + 1;
    mx2 = (hx)*60;
    
    mx = [mx1:mx2];

    try
        % start/stop times for 1 minute time length
        start_time = t_start + datenum(0,0,0,(hx-1),0,0);
        stop_time = start_time + datenum(0,0,0,1,0,0); 
        
        % start/stop times with time buffer/drift included
        t1 = find(t>=(start_time - ax3_buffer + ax3_drift),1,'first');
        t2 = find(t<(stop_time + ax3_buffer + ax3_drift),1,'last');

        T1 = t(t1:t2);
        
        if(~isempty(T1))

            % convert data to units of g
            x1 = double(data.x(t1:t2))*data.AccScale;
            y1 = double(data.y(t1:t2))*data.AccScale;
            z1 = double(data.z(t1:t2))*data.AccScale;

            acc_data = [];
            acc_data = [x1,y1,z1];

            % calculate activity counts per second at 10 Hz
            counts_buffered = [];
            counts_buffered = agfilt_felix(acc_data,T1,hz_count,AC_filterB,AC_filterA);

            % ignore activity counts from buffer data
            cut_start = sec_buffer + 1;
            cut_stop = size(counts_buffered,1) - sec_buffer;

            counts = counts_buffered(cut_start:cut_stop,:);

            % sum activity counts for each minute
            for hr_mx=1:60
                
                nx1 = (hr_mx-1)*60 + 1;
                nx2 = (hr_mx)*60;
                
                act_counts = [0,0,0];
                act_counts = [sum(counts(nx1:nx2,1)), sum(counts(nx1:nx2,2)), sum(counts(nx1:nx2,3))];
            
                ActivityCounts(mx(hr_mx),:) = act_counts;
            end
            
        else
            % impute 0 activity counts
            ActivityCounts(mx,:) = 0;
        end
        
    catch err
        % impute 0 activity counts
        ActivityCounts(mx,:) = 0;
    end

end %mx


for mx=min_rem_start:min_length
    
    try
        % start/stop times for 1 minute time length
        start_time = t_start + datenum(0,0,0,0,(mx-1),0);
        stop_time = start_time + datenum(0,0,0,0,1,0); 
        
        % start/stop times with time buffer/drift included
        t1 = find(t>=(start_time - ax3_buffer + ax3_drift),1,'first');
        t2 = find(t<(stop_time + ax3_buffer + ax3_drift),1,'last');

        T1 = t(t1:t2);
        
        if(~isempty(T1))

            % convert data to units of g
            x1 = double(data.x(t1:t2))/256;
            y1 = double(data.y(t1:t2))/256;
            z1 = double(data.z(t1:t2))/256;

            acc_data = [];
            acc_data = [x1,y1,z1];

            % calculate activity counts per second at 10 Hz
            counts_buffered = [];
            counts_buffered = agfilt_felix(acc_data,T1,hz_count,B,A);

            % ignore activity counts from buffer data
            cut_start = sec_buffer + 1;
            cut_stop = size(counts_buffered,1) - sec_buffer;

            counts = counts_buffered(cut_start:cut_stop,:);

            % sum activity counts for 1 minute
            act_counts = [0,0,0];
            act_counts = [sum(counts(:,1)), sum(counts(:,2)), sum(counts(:,3))];
            
            ActivityCounts(mx,:) = act_counts; 
            
        else
            % impute 0 activity counts
            ActivityCounts(mx,:) = 0;
        end
        
    catch err
        % impute 0 activity counts
        ActivityCounts(mx,:) = 0;
    end
    
end

