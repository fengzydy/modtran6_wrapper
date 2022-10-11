%%%%%%% Jing Feng (jing.feng@noaa.gov), Oct 10, 2022
%%%%%%% read output from .plt
function [wavnum, rad] = read_pltout(filename)
% [wavnum, rad] = read_pltout(filename)

if (nargin < 1)
  filename = './pltout';
end

fid = fopen(filename,'r');
data = fscanf(fid,'%g', [2, inf]);
fclose(fid);

wavnum = data(1,:);
rad = data(2,:);

return;
