function [metadata]=AX3_metadata_only(filename)
% [metadata]=AX3_metadata_only(filename)
% simple metadata reading of AX3 .cwa files, time is not accurate

%check filename exists
if(~exist(filename,'file'))
    disp(horzcat('Input file not found: ',filename));
    return;
end

%create data structure

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
    '',[],...
    'Time',...
    'offset'); 

%open file for reading
fid = fopen(filename,'r');

%read header packet

metadata.Packet1 = char(fread(fid,2,'char')');      %@0 +2

fseek(fid,5,'bof');
metadata.DeviceID = fread(fid,1,'uint16');          %@5 +2
metadata.SessionID = fread(fid,1,'uint32');         %@7 +4

fseek(fid,13,'bof');
tmpdt = fread(fid,1,'uint32');                  %@13 +4
metadata.StartTime = parseDT(tmpdt);
tmpdt = [];

tmpdt = fread(fid,1,'uint32');                  %@17 +4
metadata.StopTime = parseDT(tmpdt);
tmpdt=[];

fseek(fid,64,'bof');
metadata.Annotation = char(fread(fid,448,'char')'); %64 +448

fseek(fid,1024,'bof');
metadata.Packet2 = char(fread(fid,2,'char')');      %@1024+0 +2

fseek(fid,1038,'bof');
tmpdt = fread(fid,1,'uint32');                  %1024+14 +4
metadata.TimeStamp = parseDT(tmpdt);
tmpdt=[];

fseek(fid,1047,'bof');
metadata.Battery = fread(fid,1,'uint8');            %1024+23 +1

fseek(fid,1048,'bof');
RateCode = fread(fid,1,'uint8');                %1024+24 +1
metadata.SampleRate = 3200/(bitshift(1,(15-bitand(15,RateCode))));

fseek(fid,1049,'bof');
metadata.PackingFormat = fread(fid, 1, 'uint8');    %1024+25 +1

fseek(fid,1052,'bof');
metadata.SampleCount = fread(fid,1,'uint16');       %1024+28 +2

%get length of file in bytes
fseek(fid,0,'eof');                             % go to end of file
lengthBytes = ftell(fid);                       % get position
numPackets = floor(lengthBytes / 512);          % number of packets

%read timestamp data
fseek(fid,1038,'bof');
TIME = fread(fid,numPackets,'uint32',508);
metadata.Time = parseDT(TIME);
clear TIME;

%read sample offset data
fseek(fid,1050,'bof');
metadata.offset = fread(fid,numPackets,'int16',510);

%DO NOT read accelerometer data

fclose(fid);


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
