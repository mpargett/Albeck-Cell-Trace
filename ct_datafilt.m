%CT_DATAFILT
%   Filter CellTrace data based on value, derivative, time frame, or custom
%   filters.  
%
%   [d, f] = ct_datafilt(d, fType, fChan, fPar, ...)
%       filters data d (a CellTrace data structure) using the filter type,
%   fType, on Channel fChan, with filter parameter fPar. Filter definitions
%   are given below; fChan may be a vector of channel indices or a cell
%   array of channel names.  The filtered indices, f, are returned as
%   logical arrays
%
%   [d, f] = ct_datafilt(d, ft, ...)
%       filters d for a filter definition structure, with fields t, c, and
%   p for fType, fChan, and fPar as above.
%
%   [d, f] = ct_datafilt(d,..., pName, pValue)
%       additonally sets operation parameter pName to pValue.  See below
%   for available paramters.  Parameters may also be provided in a
%   structure, with pName as a fieldname containing pValue.
%
%   Filter and Parameter inputs may be entered in any order.
%
%   Filters
%       Max     - Enforce a maximum value. fPar must be a scalar.  All data
%               points greater than fPar will be rejected.
%       Min     - Enforce a minimum value. fPar must be a scalar.  All data
%               points less than fPar will be rejected.
%       dMax    - Enforce a maximum derivative. fPar must be a scalar.  All
%               data derivatives greater than fPar will be rejected.
%       dMin    - Enforce a minimum derivative. fPar must be a scalar.  All
%               data derivatives less than fPar will be rejected. dMax and
%               dMin may be used together to make an absolute value filter.
%       tWin    - Require that traces have no invalid indices within the
%               desired time range. fPar must be a 2 element array defining
%               the time range (i.e. [tStart, tEnd]).  Not for use with
%               LeaveGaps on.
%       custom  - Define a custom filter to apply.
%
%   Parameters
%       LeaveGaps   - Logical. If TRUE, instead of removing all tracks with
%                   any bad index, set bad indices to NaN and return all
%                   tracks. Run ct_trimdata afterward to remove tracks with
%                   oversize gaps and interpolate over acceptable gaps.
%
%
%   Example: A minimum value filter on the channel "RFP", and maximum value
%       filter on "DAPI"
%   [d,f] = ct_datafilt(d, 'Min', 'RFP', 0, 'Max', 'DAPI', 5000);
%
%   Example: A minimum value filter on the channel "RFP", and a required
%       time window
%   [d,f] = ct_datafilt(d, 'Min', 'RFP', 0, 'tWin', 'RFP', [20, 80]);
%
%   Example: A multi-filter call predefining the filter structure, with the
%       LeaveGaps parameter set on.
%   ft(1).t = 'Min';  ft(1).c = 'RFP';  ft(1).p = 0;
%   ft(2).t = 'Max';  ft(2).c = 'DAPI'; ft(2).p = 5000;
%   ft(3).t = 'dMax'; ft(3).c = 'DAPI'; ft(3).p = 1000;
%   [d, f] = ct_datafilt(d, ft, 'LeaveGaps', 1);


function [d, f] = ct_datafilt(d, varargin)
%List possible filter types
fm = {'Max', 'Min', 'dMax', 'dMin', 'tWin', 'custom'};
%Set default parameters
p.LeaveGaps = false;
%   Get possible paramter names
pm = fieldnames(p);    s = 1;
%   Initialize filter structure
ftn = {'t', 'c', 'p'};  ft = cell2struct(cell(1,3), ftn, 2);

%% Parse filter and parameter inputs
sf = 1;  nin = numel(varargin); 
while s < nin    
    %Allow filter and parameter structures
    if isstruct(varargin{s})	
        if all(isfield(varargin{s}, ftn)) %IF matching a filter, assign
            ft(s + 1:numel(varargin{s})) = varargin{s};  
        elseif all(isfield( varargin{s}, pm )) %IF matching parameters
            for vm = fieldnames(varargin{s})'   %Apply each parameter
                p.(vm{1}) = varargin{s}.(vm{1});
            end
        else    %Error if a structure input not matching proper fields
            error('CT:DATAFILT:BadInput', ['An input structure, either ',...
                'filter or parameter, contains an invalid field']);
        end; s = s+1; %Increment to next input
        
    %Allow input of Name/Value pairs (or triplets, for filters)
    elseif ischar(varargin{s})  %IF a string, use as Name/Value
        if any(strcmpi(varargin{s}, fm));   %Check name against filters
            [ft(sf).t, ft(sf).c, ft(sf).p] = deal(varargin{s:s+2}); %Add
            sf = sf + 1;    s = s + 3;      %Increment counters
        else pnm = strcmpi(varargin{s}, pm);   %Check against parameters
            if any(pnm);  p.(pm{pnm}) = varargin{s+1};  s = s + 2; %Add
            else warning('CT:DATAFILT:BadInput', [num2str(varargin{s}), ...
                    'is not a recognized filter type or parameter.']);
            end
        end
    end
end
%   Get number of filters
nf = numel(ft);


%% Input data type management
[d, iv, dp, ci] = ct_compform(d, [], [], {ft.c});
[ft.c] = deal(ci{:});   %Apply mapped channel indices
%Data are NaN-padded arrays, in a cell array
nx = numel(d);


%% Perform filter calculations
%   Identify any derivitive filters
dfs = ~cellfun(@isempty, regexpi({ft.t}, '^dm'));

%Loop over XYs
f = cell(nx,1);
for s = 1:nx
    [nc, nt, nchan] = size(d{s});   %Get size of data array
    
    if any( dfs )   %IF any filters are derivative
        sc = unique([ft(dfs).c]);       %Get relevant channels
        dtmp = nan(nc, nt, nchan);  	%Initialize derivative matrix
        dtmp(:,2:end,sc) = diff( d{s}(:,:,sc), 1, 2);  %Fill channels
    end
    
    %Loop over filters requested to derive each
    f{s} = false(nc, nt, nchan);  %Initialize filter cell array
    for sf = 1:nf
        %Determine filtration indices based on type, parameters
        switch lower(ft(sf).t)
            case 'max'
                f{s}(:,:,sf) = any( d{s}(:,:,ft(sf).c) > ft(sf).p, 3);
            case 'min'
                f{s}(:,:,sf) = any( d{s}(:,:,ft(sf).c) < ft(sf).p, 3);
            case 'dmax'
                f{s}(:,:,sf) = any( dtmp(:,:,ft(sf).c) > ft(sf).p, 3);
            case 'dmin'
                f{s}(:,:,sf) = any( dtmp(:,:,ft(sf).c) < ft(sf).p, 3);
            case 'twin'
                    f{s}(:,ft(sf).p(1):ft(sf).p(2),sf) = any( isnan(...
                        d{s}(:, ft(sf).p(1):ft(sf).p(2), ft(sf).c) ), 3);
            case 'custom'
                f{s}(:,:,sf) = any( ft(sf).p(d{s}(:,:,ft(sf).c)), 3);
        end
    end
    %Compile filtration indices (across filters), inverts to good indices
    f{s} = ~any(f{s},3);
    
    %Eliminate bad track or only bad data if gaps allowed
    if p.LeaveGaps
        d{s}( repmat(~f{s},1,1,nchan) ) = NaN;	%Convert bad indices to NaN
    else
        %Get index for tracks with violations and filter
        f{s} = all(f{s},2);       d{s} = d{s}(f{s},:,:);
        %Filter index vector (linfo)
        if ~isempty(iv.lin); 	%IF linfo was present
            iv.lin{s} = iv.lin{s}(f{s});	%Remove bad indices
            %Map to new indices
            imap = nan(nc,1);   imap(f{s}) = 1:nnz(f{s}); 
            %Adjust new index vector references
            for sl = unique(iv.lin{s}(~isnan(iv.lin{s})))';   
                iv.lin{s}(iv.lin{s} == sl) = imap(sl);  end
        end
        %Filter cell index to remove bad indices
        if ~isempty(iv.cel); iv.cel{s} = iv.cel{s}(f{s}); end 
    end
    
end     %END Loop over XYs 

%Reframe data to structure provided at input
[d] = ct_compform(d, iv, dp);


end


