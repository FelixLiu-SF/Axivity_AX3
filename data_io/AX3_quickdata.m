function [data,metadata]=AX3_quickdata(filename)
% [data,metadata]=AX3_quickdata(filename)
% simple reading of AX3 .cwa files, time is not accurate

%check filename exists
if(~exist(filename,'file'))
    disp(horzcat('Input file not found: ',filename));
    return;
end

%create data structure
data = struct(...
    'Packet1',[],...
    'DeviceID',[],...
    'SessionID',[],...
    'StartTime',[],...
    'StopTime',[],...
    'Annotation',[],...
    'Packet2',[],...
    'Battery',[],...
    'TimeStamp',[],...
    'Temperature',[],... % centigrade
    'SampleRate',[],...
    'SampleCount',[],...
    'PackingFormat',[],...
    'e',[],...
    'x',[],...
    'y',[],...
    'z',[],...
    'Time',...
    'offset'); 

metadata = struct(...
    'Packet1',[],...
    'DeviceID',[],...
    'SessionID',[],...
    'StartTime',[],...
    'StopTime',[],...
    'Annotation',[],...
    'Packet2',[],...
    'Battery',[],...
    'TimeStamp',[],...
    'SampleRate',[],...
    'SampleCount',[],...
    'PackingFormat',[]); 

%open file for reading
fid = fopen(filename,'r');

%read header packet

data.Packet1 = char(fread(fid,2,'char')');      %@0 +2

fseek(fid,5,'bof');
data.DeviceID = fread(fid,1,'uint16');          %@5 +2
data.SessionID = fread(fid,1,'uint32');         %@7 +4

fseek(fid,13,'bof');
tmpdt = fread(fid,1,'uint32');                  %@13 +4
data.StartTime = parseDT(tmpdt);
tmpdt = [];

tmpdt = fread(fid,1,'uint32');                  %@17 +4
data.StopTime = parseDT(tmpdt);
tmpdt=[];

fseek(fid,64,'bof');
data.Annotation = char(fread(fid,448,'char')'); %64 +448

fseek(fid,1024,'bof');
data.Packet2 = char(fread(fid,2,'char')');      %@1024+0 +2

fseek(fid,1038,'bof');
tmpdt = fread(fid,1,'uint32');                  %1024+14 +4
data.TimeStamp = parseDT(tmpdt);
tmpdt=[];

fseek(fid,1047,'bof');
data.Battery = fread(fid,1,'uint8');            %1024+23 +1

fseek(fid,1048,'bof');
RateCode = fread(fid,1,'uint8');                %1024+24 +1
data.SampleRate = 3200/(bitshift(1,(15-bitand(15,RateCode))));

fseek(fid,1049,'bof');
data.PackingFormat = fread(fid, 1, 'uint8');    %1024+25 +1

fseek(fid,1052,'bof');
data.SampleCount = fread(fid,1,'uint16');       %1024+28 +2

%get length of file in bytes
fseek(fid,0,'eof');                             % go to end of file
lengthBytes = ftell(fid);                       % get position
numPackets = floor(lengthBytes / 512);          % number of packets

%read timestamp data
fseek(fid,1038,'bof');
TIME = fread(fid,numPackets,'uint32',508);
data.Time = parseDT(TIME);
clear TIME;

%read temperature data
fseek(fid,1044,'bof');                          %1024+20 +2
RAWTEMP = fread(fid,numPackets,'uint16',510);
data.Temperature = (double(RAWTEMP)*150.0 - 20500)/1000;
clear RAWTEMP;

%read sample offset data
fseek(fid,1050,'bof');
data.offset = fread(fid,numPackets,'int16',510);

%read accelerometer data
fseek(fid,1054,'bof');

if data.PackingFormat == (3*16 + 0) % packed data, 120 samples, 32-bit
    
    ACC = fread(fid,[data.SampleCount,Inf],horzcat(num2str(data.SampleCount),'*uint32'),32);
    ACCsz = size(ACC);
    ACC = reshape(ACC,(ACCsz(1)*ACCsz(2)),1);

    e = int32(2 .^ bitshift(ACC, -30));
    x = int32(bitand(         ACC      , 1023)); 
    y = int32(bitand(bitshift(ACC, -10), 1023)); 
    z = int32(bitand(bitshift(ACC, -20), 1023)); 

    data.e = e;
    data.x = e .* (x - (int32(1024) .* int32(x >= 512)));
    data.y = e .* (y - (int32(1024) .* int32(y >= 512)));
    data.z = e .* (z - (int32(1024) .* int32(z >= 512)));

elseif data.PackingFormat == (3*16 + 2) % unpacked data, 80 samples, 48-bit
    
    ACC = fread(fid,[data.SampleCount,Inf],horzcat(num2str(3*data.SampleCount),'*int16'),32);
    ACCsz = size(ACC);
    ACC = reshape(ACC,(ACCsz(1)*ACCsz(2)),1);

    data.e = 1; % no modulation in unpacked data format
    data.x = ACC(1:3:(end-2));
    data.y = ACC(2:3:(end-1));
    data.z = ACC(3:3:end);
    
else % unknown data format
    
    data.e = [];
    data.x = [];
    data.y = [];
    data.z = [];
    
end

fclose(fid);

metadata.Packet1 = data.Packet1;
metadata.DeviceID = data.DeviceID;
metadata.SessionID = data.SessionID;
metadata.StartTime = data.StartTime;
metadata.StopTime = data.StopTime;
metadata.Annotation = data.Annotation;
metadata.Packet2 = data.Packet2;
metadata.Battery = data.Battery;
metadata.TimeStamp = data.TimeStamp;
metadata.SampleRate = data.SampleRate;
metadata.SampleCount = data.SampleCount;
metadata.PackingFormat = data.PackingFormat;


function datetime = parseDT(t)
%parse packed date/time

datearr = zeros(size(t,1),6);

datearr(:,1) = bitshift(t,-26) + 2000;
datearr(:,2) = bitshift(bitand(t,67108863),-22);
datearr(:,3) = bitshift(bitand(t,4194303),-17);
datearr(:,4) = bitshift(bitand(t,131071),-12);
datearr(:,5) = bitshift(bitand(t,4095),-6);
datearr(:,6) = bitand(t,63);

datetime = datenum(datearr);
