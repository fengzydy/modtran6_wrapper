%%%%%%% Jing Feng (jing.feng@noaa.gov), Oct 10, 2022
%%%%%%% read output from .chn
function [wavenum,rad]=read_channels(filename,version)
fid=fopen(filename);
if ~exist('version','var')
      version=5;
end
if version==5
data=textscan(fid,'%12.6f%5d%15f%15f%10.5f%15f%11.5f%11.5f%11.5f%11.5f%15.7f%15.7f%15.7f%15.7f%12.8f%12.8f%9s%12f%5s%8s%8f%5s',...
              'Headerlines',5);
        wavenum=data{1};
        rad=data{3};
else
      data=textscan(fid,'%12.6f%5d%5d%15f%15f%10.5f%15f%11.5f%11.5f%11.5f%11.5f%15.7f%15.7f%15.7f%15.7f%12.8f%12.8f%9s%12f%5s%8s%8f%5s',...
              'Headerlines',5);
        wavenum=data{1};
        rad=data{4};
end


fclose(fid)