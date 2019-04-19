

% --------------------------------------------------------------------
% function to load 2D model data
% --------------------------------------------------------------------


function [lat, lon, varargout] = mf_load(file_name, varargin)


% check if file is there
file_exist = exist(file_name,'file');


% 1. if file is there:
if file_exist == 2

    
% load grid data
% try
    lat = flipud(ncread(file_name, 'lat')); % latitude [° N] (axis: Y)
    lon =        ncread(file_name, 'lon') ; % longitude [° E] (axis: X)
% catch errmessage
%     disp(errmessage.message)
%     disp('loading -LATIXY- and -LONGXY- instead')
%     lat = flipud(ncread(file_name, 'LATIXY')); % latitude [° N] (axis: Y)
%     lon =        ncread(file_name, 'LONGXY') ; % longitude [° E] (axis: X)
% end


% circshift lon data so that Greenwich is in the center of the matrix
shiftsize       = length(lon)./2;
lon             = circshift(lon, shiftsize);
lon(lon >= 180) = lon(lon >= 180) - 360;     % set range tp [-180 177.5]


% get 2D grid
[lon, lat] = meshgrid(lon, lat);


% loop over variable names
for i=1:length(varargin)

    
    % load variable data
    VAR = varargin{i};            % get variable name
    VAR = ncread(file_name, VAR); % load the variable


    % check variable dimensions and tread accordingly
    [nx ny nz nt] = size(VAR);
    if     numel(size(VAR)) == 2
        VAR = rot90(VAR); 
    elseif numel(size(VAR)) == 3
        VARr = NaN(ny,nx,nz);
        for j=1:size(VAR,3)
            VARr(:,:,j) = rot90(VAR(:,:,j));
        end
        VAR = VARr;
    elseif numel(size(VAR)) == 4
        VAR = permute(VAR,[1 2 4 3]);
        % check dimensions again after permute command
        if     numel(size(VAR)) == 3    % it was a 'fake' fourth dimension (e.g. U10,V10)
            VARr = NaN(ny,nx,nz);
            for j=1:size(VAR,3)
                VARr(:,:,j) = rot90(VAR(:,:,j));
            end
            VAR = VARr;
        elseif numel(size(VAR)) == 4    % it was a real fourth dimension (e.g. QV_hm)

            % undo permute (result: vertical: 3th, time: 4th)
            VAR = permute(VAR,[1 2 4 3]);
            VARr = NaN(ny,nx,nz,nt);
            for j=1:size(VAR,3)
                for k=1:size(VAR,4)
                    VARr(:,:,j,k) = rot90(VAR(:,:,j,k));
                end
            end
            VAR = VARr;
        end
    end

    
    % circshift data so that Greenwich is in the center of the matrix
    VAR = circshift(VAR, [0 shiftsize]);

    
    % store data
    varargout{i} = VAR;
    
end


% 2. if file is not there: send out warning
else
disp(['file ', file_name,' does not exist'])
end


end

