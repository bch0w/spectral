% convert a grd file into .mat format for mesh creation
% SRTM30_PLUS has 30 second resolution which is roughly 1 km
% e.g. (6371 * 2pi)/360 * 1/60 * 1/60 * 30 = 0.927
% e140s10.nc is a section of SRTM30_PLUS containing parts of oceania
% LLHC 140E60S / URHC 180E10S

[lon_temp,lat_temp,Zt] = grdread2('e140s10.nc');

Zt = double(Zt); 

[ny, nx]=size(Zt);
Xlon = repmat(lon_temp,ny,1);
Ylat = repmat(lat_temp',1,nx);

% taken from Carl's scripts, not sure what the latter values are but 'd'
% is not used so leaving it how it is
d = [-180,180,-90,90,-10898,8271,0,0.0166666666666667,0.0166666666666667];

save('srtm30_e140s10','Xlon','Ylat','Zt','nx','ny','d');

