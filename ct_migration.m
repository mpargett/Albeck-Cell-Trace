%CT_MIGRATION
%   Calculates cell motion metrics from X and Y coordinates in CellTrace
%   data.
%
%   z = ct_migration(x, ...)
%       returns motion metrics in output structure z, for the data in input
%       x (either the path/filename of the processed data file, or its
%       contents in a structure).  Process parameters are passed as
%       Name/Value pairs.
%
%   Output:
%   z - structure containing fields:
%       vel     - nCell x nTime-1 x 2 array of cell velocity 3rd dimension
%                   holds x, then y, velocities 
%       spd     - nCell x nTime-1 array of cell speed
%       ang     - nCell x nTime-1 array of cell angular direction 
%       turn    - nCell x nTime-2 array of cell turning angles
%       dnet    - vector of net distances traveled per cell
%       snet    - vector of net speed
%       anet    - vector of net angular direction
%       dtot    - vector of total distance traveled per cell
%       stot    - vector of average speed per cell
%
%   Parameters:
%       xysc    - scaling for the X and Y coordinates (e.g. to absolute
%                   units). xysc may be scalar (scaling both X and Y by the
%                   same value) or a two element array ([xscale, yscale]). 
%       xf      - (for post data-proc'ed input), fieldname of X coordinates
%       yf      - (for post data-proc'ed input), fieldname of Y coordinates
%       plot    - Logical, TRUE to also plot traces via ct_migrationplot.
%                   All parameters will be forwarded, so plotting
%                   parameters may also be included.
%       ploth   - Handle to axes on which to plot the migration traces.
%
%
%   Example usages:
%       --- Calling filename, without scaling
%   z = ct_migration('L:\Processed Data\FolderName\Filename_xy01.mat');
%   
%       --- Pre-loading data, with scaling from MetaData
%   gbl = load('L:\Processed Data\FolderName\Filename_Global.mat');
%   xysc = [gbl.GMD.cam.ip.bkmd.cam.PixSizeX, ...
%           gbl.GMD.cam.ip.bkmd.cam.PixSizeY];
%   in = load('L:\Processed Data\FolderName\Filename_xy01.mat');
%   z = ct_migration(in, 'xysc', xysc);
%
%       --- Also plot the migration tracks (in desire axes)
%   figure; ah = axes;
%   z = ct_migration('L:\Processed Data\FolderName\Filename_xy01.mat', ...
%       'plot', true, 'ploth', ah);
%


function z = ct_migration(x, varargin)
%% Parse parameters 
%Initialize as default
p.xysc = 1;     %Spatial scaling parameter (default no scaling)
p.xf = [];  p.yf = [];  %Fieldnames for data structure X and Y
p.plot = false;         %Logical to also plot tracks
p.ploth = [];           %Plotting handle handle

%Input option pair parsing:
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end
%Split pairs to the parameter structure
for s = 1:2:nin;    p.(lower(varargin{s})) = varargin{s+1};     end

%% Parse input data
%If numeric, assume nCells x nTime x [X,Y]
if isnumeric(x);    xi = 1; yi = 2;
else
    %If filename provided, load data
    if ischar(x);   x = load(x, 'vcorder', 'valcube'); 	end
    if isstruct(x)
        if isfield(x,'valcube') %IF loaded from IMAN processed data
            xi = find(strcmpi('xcoord', x.vcorder));    %X indices
            yi = find(strcmpi('ycoord', x.vcorder));    %Y indices
            x = x.valcube;      %Data matrix
        elseif isfield(x,'data') %IF from post-proc'ed data
            x = cat(3, x.data.(p.xf), x.data.(p.yf));    %Keep coordinates
            xi = 1; yi = 2;                             %X and Y indices
        end
    end
end


%% Process input data
%   Get data size
nc = size(x,1);
%Get number of time points per cell
ntime = sum(~isnan(x(:,:,xi)),2);

%Scale XY coordinates if indicated
if numel(p.xysc) == 1     %Same scale for each
    x(:,:,[xi,yi]) = x(:,:,[xi,yi]).*p.xysc;
elseif numel(p.xysc) == 2 %Individual scales
    x(:,:,xi) = x(:,:,xi).*p.xysc(1);
    x(:,:,yi) = x(:,:,yi).*p.xysc(2);
end

%Loop over cells to get first and last positions (skipping NaNs)
[xyst, xynd]  = deal(nan(nc,2));
for s = 1:nc
    %Get valid time points (skip if no valid points)
    vt = ~isnan(x(s,:,xi));    if ~any(vt); continue; end
    vst = find(vt, 1, 'first');     vnd = find(vt, 1, 'last');
    
    %Get x,y start and ending values
    xyst(s,1:2) = x(s,vst,[xi,yi]);
    xynd(s,1:2) = x(s,vnd,[xi,yi]);
end


%% Calculate and assign output metrics
%Velocity, speed and direction (as functions of time)
%   Velocity vector
z.vel = [diff(x(:,:,[xi,yi]), 1, 2), NaN(nc,1,2)];
%   Speed (magnitude of movement)
z.spd = sqrt(nansum(z.vel.^2, 3));
%   Angle of movement
z.ang = atan2(z.vel(:,:,2), z.vel(:,:,1));
%   Turning rate (angular)
z.turn = [NaN(nc,1), diff(z.ang,1,2)];

%Final motion metrics (cumulative per cell)
%   Total distance moved (integrated over time)
z.dtot = nansum(z.spd, 2);
%   Total (mean) speed
z.stot = z.dtot./ntime;

%   Net velocity
vnet = (xynd - xyst);
%   Net distance moved
z.dnet = sqrt( nansum( vnet.^2, 2 ) );
%   Net speed
z.snet = z.dnet./ntime;
%   Net angle
z.anet = atan2(vnet(:,2), vnet(:,1));

%% Run Migration Plot
if p.plot
    vpass = [{x(:,:,xi), x(:,:,yi)}, varargin];
    if ishandle(p.ploth);   vpass = [{p.ploth},vpass];    end
    if isfield(p, 'framesize'); p.framesize = p.framesize(:).*p.xysc(:); 
    vpass = [vpass, {'framesize', p.framesize}];     end
    ct_migrationplot(vpass{:});
end    
    
end
