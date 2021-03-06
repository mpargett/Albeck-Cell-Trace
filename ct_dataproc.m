%CT_DATAPROC
%   Processing for Cell Trace data, from raw values (nCells x nTime x
%   nChannels) to named outputs.  Processing functions used are
%   user-defined, or selected from commonly used reporter models.
%   
%   [d,p] = ct_dataproc(fname, ...)
%       delivers outputs d, the prepared data (a structure with fields
%   'cellindex', 'data', and 'linfo' if present), and p, the parameters
%   associated. 
%
%   Parameters may be provided as name, value pairs, as in
%   [d,p] = ct_dataproc(fname, name1, value1, name2, value2, ...)
%   
%   savename    - cell array of names under which to save processed data
%   name        - cell array of names of channels ({'c1','c2'})
%   ind         - cell array in indices from valcube to use for each
%                channel ({[1,2], 4})
%   fun         - cell array of functions to apply to each channel
%   tsamp       - sampling period (not actually used)
%   pwrat       - power ratio for FRET channel conversion
%   gapmax      - maximum number of bad time points to bridge
%   dlegnthmin  - minimum length of an acceptable trace
%   startonly   - TRUE to only keep tracks that begin at the start of imaging
%   packascell  - TRUE to return data with NaNs removed as a cell array
%                (depricated), FALSE will return each channel as a
%                NaN-padded  cells x time array 
%   keeplinfo   - Set to FALSE to discard Lineage info (default: TRUE)
%   verbose     - TRUE to show all warnings

function [d,p] = ct_dataproc(fname, varargin)

%% Pre-processing
%Default parameters
p = struct( 'savename',{{[]}},  'name',{{'ekar','fra','x','y'}}, ...
    'ind',{{[4,5],3,7,8}}, 'fun',{{[]}},  'tsamp',1,  'pwrat',[],  ...
    'stripidx', false, 'verbose', true, ...
	'gapmax',2,  'dlengthmin',100,  'startonly',false, ...
    'keeplinfo', true); %<- For trimming

%Input option pair parsing:
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end
%Split pairs to the parameter structure
for s = 1:2:nin;    p.(lower(varargin{s})) = varargin{s+1};     end

%Block warnings if not running verbose
if p.verbose;   warning('ON', 'CT:DATAPROC');
else            warning('OFF','CT:DATAPROC');   end

%Check input file name
if ischar(fname); fname = {fname}; end
if ischar(p.savename); p.savename = {p.savename}; end
nf = length(fname);
if nf ~= numel(p.savename)
    warning('CT:DATAPROC', ['Number of file names should match number',...
        ' of save names.  Any file names with no corresponding save ',...
        'name will not be saved.']);
    [p.savename{end:numel(fname)}] = deal([]);
end

%Check parameter structure
if ischar(p.name); p.name = {p.name}; end;      nn = numel(p.name);
if length(p.fun) ~= nn; [p.fun{end+(1:(nn-length(p.fun)))}] = deal([]); end
if isnumeric(p.ind);    
    if nn == 1; p.ind = {p.ind}; 
    elseif length(p.ind) == nn; p.ind = num2cell(p.ind); 
    else error('Number of Names (name) and Indices (ind) must match'); 
    end
end


%% Refine data
dc = cell(nf,1);
for s = 1:length(fname)
    if isempty(fname{s});	continue;   %Return empty cell for bad input
    elseif ~exist(fname{s},'file')
        warning('CT:DATAPROC',['Invalid filename encountered.  ',...
            'Resulting cell (',num2str(s),') will be empty.']); continue;
    end
    %Load data - check first if lineage info is present
    vchk = whos('-file',fname{s});  
    haslinfo = any(strcmp({vchk.name},'linfo')) && p.keeplinfo;
    if haslinfo     %Collect linfo, if available
            td = load(fname{s}, 'vcorder', 'valcube', 'bkg', 'linfo');
    else    td = load(fname{s}, 'vcorder', 'valcube', 'bkg');
    end;    rmbkg = false;  %Initialize bkg removal as FALSE
    
    %If note indicates individual channels are raw, allow bkg removal
    %   Note, this is only included for legacy -raw channels are depricated
    if isfield(td,'vcorder') && ...
       ~isempty(regexpi(td.vcorder{end}, 'indiv.*are[^(no)]*raw'))
   
        indivs = ~cellfun('isempty', ... %ID indiv channels
            regexpi(td.vcorder, '^.fp_(nuc|cyt)[^/]*$'));
        if any(indivs([p.ind{:}])); rmbkg = true;  nch = nnz(indivs)/2;
        if isempty(td.bkg)
            tn = regexp(fname{s}, ['(?<d>^.*\\)(?<f>[^\\]*)',...
                '_xy(?<n>\d*)(?<a>.*\.mat)'], 'names');
            gmd = load([tn.d,tn.f,'_Global.mat']);
            axy = gmd.p.bk(str2double(tn.n)).altxy;
            if axy < 10; dh = '0'; else dh = []; end
            tmpd = load([tn.d,tn.f,'_xy',[dh,num2str(axy)],tn.a], 'bkg');
            td.bkg = tmpd.bkg;   %Replace bkg with alternate
        end
        end
    else
    end
    
    %Extract desired channels prior to cleanup
    for sc = 1:nn
        if isempty(p.fun{sc});  
            p.fun{sc} = subf_deffuns(p.name{sc}, p.ind{sc}, p.pwrat); 
        end
        intemp = arrayfun(@(x)td.valcube(:,:,x), p.ind{sc}, ...
            'UniformOutput', false);
        %Remove background from individidual channels (as needed)
        if rmbkg && any( indivs(p.ind{sc}) )
            rbi = indivs(p.ind{sc});
            %Get relevant channel indices and remove background
            chi = mod(p.ind{sc}-1,nch)+1;
            intemp(rbi) = cellfun(@(x,b)bsxfun(@minus,x,b), intemp(rbi), ...
                num2cell(td.bkg(:,chi),1), 'UniformOutput', false);
        end
        %Apply processing function and retain result
        pd.(p.name{sc}) = p.fun{sc}(intemp{:}); 
    end
    %Run data trimming procedure
    d = ct_trimdata(pd, p);
    
    %IF linfo is present, apply trimming index to linfo, append to output
    if haslinfo && ~isempty(d.cellindex) && size(td.linfo,1) >= max(d.cellindex)
        %Filter linfo
        linfo = td.linfo(d.cellindex,:);
        %   Build map to new indices
        imap = nan(d.cellindex(end),1);          %Initialize 
        imap(d.cellindex) = 1:nnz(d.cellindex);     %Place new index values
      
        %Adjust linfo references
        for sl = reduce_empties(unique(linfo(~isnan(linfo)))')
            linfo( linfo == sl ) = imap(sl);
        end
        %Append linfo to the output
        d.linfo = linfo;
    end
    
    %Save data IF filename provided
    if ~isempty(p.savename{s});    save(p.savename{s},'d','p'); end
    
    %Store data from this XY in a cell
    dc{s} = d;     %Store data (and indices)
    %   Index stripping now disallowed
%     if p.stripidx; dc{s} = d.data;  %Store data
%         if haslinfo; dc{s}.linfo = d.data.linfo; end    %Include linfo
%     else dc{s} = d;     %Store data (and indices)
%     end
end
d = dc; %Assign final output

end


%Subfunction for default function definitions
function sfun = subf_deffuns(n, ind, pwrat)

%Search for pre-defined names
rn = {'kar', 'ktr'};    %Recognized names
rmat = find( ~cellfun('isempty', regexpi(n, rn, 'start')) );

%Check for inconsistent matches
nr = length(rmat);  dtxt = 'Applying no processing function.';
if nr < 1;     warning('CT:DATAPROC',['Name "',n,'" not recognized.  ',dtxt]); rmat = 0;
elseif nr > 2; warning('CT:DATAPROC',['Unparsable name, "',n,'".  ',dtxt]); rmat = 0;
end

%Define desired processing function
switch rmat
    case 1  %For Kinase Activity Reporters (KAR, FRET-based)
        if ~isempty(pwrat); 
            if numel(ind) == 2;     sfun = @(x1,x2) 1 - (x1./x2)./pwrat;
            elseif numel(ind) == 1; sfun = @(x1)    1 - (x1)    ./pwrat; 
            else warning('CT:DATAPROC',['FRET-based data must passed ',... 
                'with 1 or 2 channel indices only.  (1 is for ',...
                'pre-computed ratios, 2 for individual color channels)']);   
                sfun = @(x)x;
            end
        else warning('CT:DATAPROC',['Power ratio must be supplied for ',...
        	'FRET-based Kinase Activity Reporters ("',n,'").  ',dtxt]);  
            sfun = @(x)x;
        end
    case 2  %For Kinase Translocation Reporters (KTR)
        if numel(ind) == 2;     sfun = @(x1,x2)x1./x2;  %Take ratio
        elseif numel(ind) == 1; sfun = @(x1)   x1;  %Do nothing, if 1 index
            warning('CT:DATAPROC',['KTR "',n,'" provided only one channel',...
                ' index (ok if input channel is a precomputed ratio.', dtxt]);
            else warning('CT:DATAPROC',...
                    ['No Channel index found for ("',n,'").  ',dtxt]);  
                sfun = @(x)x;
        end
    otherwise;  sfun = @(x)x;
end

end

