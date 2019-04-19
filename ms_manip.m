
% --------------------------------------------------------------------
% subroutine to perform manipulations on loaded variables
% note: preferably run "main"
% --------------------------------------------------------------------



% --------------------------------------------------------------------
% manipulations: general
% --------------------------------------------------------------------


% get grid dimensions
[nlat, nlon] = size(lat_mod);


% get number of srex regions
nreg = length(srex_reg);


% get number of probabilities(=percentages)
npercentages = length(percentages);


% get land pixel indices - excluding Antarctica
island                = pct_land > 50; 
island(lat_mod < -60) = 0;             % remove antarctica


% get land pixel indices - excluding Antarctica
isirr                = pct_irr > 10;
isirr(lat_mod < -60) = 0;            % remove antarctica


% get PINB (Pakistan, India Bangladesh) indices
countries      = m_shaperead('ne_10m_admin_0_countries'); % load country borders
countries.p    = find(ismember(countries.NAME,'Pakistan'));
countries.i    = find(ismember(countries.NAME,'India'));
countries.n    = find(ismember(countries.NAME,'Nepal'));
countries.b    = find(ismember(countries.NAME,'Bangladesh'));
countries.p_in = inpolygon(lon_mod, lat_mod, countries.ncst{countries.p}(:,1), countries.ncst{countries.p}(:,2) );
countries.i_in = inpolygon(lon_mod, lat_mod, countries.ncst{countries.i}(:,1), countries.ncst{countries.i}(:,2) );
countries.n_in = inpolygon(lon_mod, lat_mod, countries.ncst{countries.n}(:,1), countries.ncst{countries.n}(:,2) );
countries.b_in = inpolygon(lon_mod, lat_mod, countries.ncst{countries.b}(:,1), countries.ncst{countries.b}(:,2) );
ispinb         = island & ( countries.p_in | countries.i_in | countries.n_in | countries.b_in );
save('mw_ispinb','ispinb');

% get boundary around PINB for plotting polygon on map
lon_ispinb = double(lon_mod(ispinb));
lat_ispinb = double(lat_mod(ispinb));
ispinb_boundary = boundary(lon_ispinb, lat_ispinb, 1);


% create date vectors
date_vec = datevec(datenum(time_begin):1:datenum(time_end));


% get corners of SREX regions used in this study and use them to generate srex masks
for i=1:nreg
    [~, pol_lat(i), pol_lon(i)] = mf_srex(lat_mod, lon_mod, TXx_ctl, area, island, [], srex_reg{i}); %#ok<*SAGROW>
    issrex(:,:,i)               = inpolygon(lon_mod, lat_mod, pol_lon{i}, pol_lat{i}) & island;
end


