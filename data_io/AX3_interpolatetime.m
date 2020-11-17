function [t]=AX3_interpolatetime(data)
% Interpolate the AX3 time & offset data to the recorded accelerometry
% [t]=AX3_interpolatetime(data);

%align AX3 offsets
ix_offset = (data.SampleCount*([1:size(data.offset,1)]' - 1)) + data.offset; %index of offsets are the true TimeStamps


%get rid of duplicate offset index stamps
first_duplicate_ix = find(diff(ix_offset)<1); 
ix_offset(first_duplicate_ix,:) = [];
data.Time(first_duplicate_ix,:) = [];

%interpolate time
t = interp1(ix_offset,data.Time,[1:size(data.x,1)],'linear'); %interpolate inbetween true times

%extrapolate both ends using SampleRate
mat_sec = 1/(24*60*60); 
if(ix_offset(1)>1)
    t(1:ix_offset(1)) = (t(ix_offset(1)) - (mat_sec)*((ix_offset(1)-1)/data.SampleRate)):(mat_sec/data.SampleRate):t(ix_offset(1));
end
if(ix_offset(end)<size(t,2))
    t(ix_offset(end):end) = (t(ix_offset(end))):(mat_sec/data.SampleRate):(t(ix_offset(end)) + (mat_sec)*(data.SampleCount-data.offset(end))/data.SampleRate);
end