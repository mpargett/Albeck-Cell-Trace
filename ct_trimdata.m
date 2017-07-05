%CT_TRIMDATA
%   Routine to trim Cell Trace data to useable segments, and optionally
%   fill small gaps via interpolation.
%
%   out = ct_trimdata(in, param)
%   
%   Returns a structure with fields 'data', the trimmed data, presented as
%   provided (i.e. in a structure, or cell array per channel), 'hasdata',
%   and index of the cells that contained usable trace, and if requesting
%   outputs as a cell array, 'timei' the indices of time [start, end]
%   provided in the output.
%
%   Parameters:
%   gapmax
%   dlengthmin
%   startonly
%   packascell

function out = ct_trimdata(in, param)

%Manage input type, structure or cell
if isstruct(in); infnames = fieldnames(in);  in = struct2cell(in); 
    skp = cellfun(@(x)~(isnumeric(x) || iscell(x)), in); 
    inpass = in(skp); in = in(~skp);
end;   if isnumeric(in); in = {in};   end;   nchan = numel(in);

%Manage operational parameters
%   Default values
p = struct('gapmax', 10, 'dlengthmin', 100, 'startonly', true, ...
    'packascell', false);
%   Assign provided parameters
if exist('param','var') && ~isempty(param);  pnames = fieldnames(param);
    for s = pnames(:)'; p.(lower(s{1})) = param.(s{1}); end
end

%FOR each channel that is a cell array, reconstitute as NaN padded array
isc = cellfun(@iscell,in);
in(isc) = cellfun(@(x)sub_nanpad(x), in(isc), 'UniformOutput', false);

%Define indices of 'good' data (nonzero), in all channels simultaneously
goodi = cellfun(@(x)~isnan(x), in, 'UniformOutput', false);
if all(cellfun(@(x)all(x(:)), goodi)); %FIXME. When other methods nailed down better
    goodi = cellfun(@(x)x>0, in, 'UniformOutput', false);
end

for s = 2:numel(goodi); goodi{s} = goodi{s-1} & goodi{s}; end
goodi = goodi{end};
[nCells, nTime] = size(goodi);

%Indicate data range as region between first and last data in each Cell
drng = false(nCells,nTime);
dst = goodi & ~[false(nCells,1), goodi(:,1:end-1)];
dnd = goodi& ~[goodi(:,2:end), false(nCells,1)];
for s = 1:nCells
    if p.startonly  %IF only wants starts, consider start always in range
        drng(s, 1 : find(dnd(s,:),1,'last') ) = true;
    else  %ELSE take bad start time out of range
        drng(s, find(dst(s,:),1,'first') : find(dnd(s,:),1,'last') ) = true;
    end
end


%% TRIM DATA AND IDENTIFY SMALL GAPS TO FILL
gst = cell(nCells,1); gnd = gst;
% in = cellfun(@(x)num2cell(x,2), in, 'UniformOutput', false);
for s = 1:nCells
    %Identify gaps as being in the data range, but not a good index
    gaps = ~goodi(s,:) & drng(s,:);
    if any(gaps)
        %Find indices demarking beginning and end of each gap
        gst{s} = find( gaps & [false, ~gaps(1:end-1)] );
        gnd{s} = find( gaps & [~gaps(2:end), false] );
        
        %Check for starting gaps, if startonly
        if p.startonly && length(gnd{s}) > length(gst{s})
            if gnd{s}(1) > p.gapmax; drng(s,:) = false; continue;
                %IF starting gap is more than max, neglect trace, else...
            else  drng(s,1:gnd{s}(1)) = true;   %Set range as ok
                for sc = 1:nchan    %Fill start with first data value
                    in{sc}(s, 1:gnd{s}(1)) = ...
                        in{sc}(s, gnd{s}(1)+1);
                end
                gnd{s}(1) = [];     %Remove gap end flag
            end
        end
        
        %Proceed with gap evaluation
        gapsz = gnd{s}-gst{s} + 1;    %Size of gaps
        if any(gapsz > p.gapmax)  %IF any gaps are too big
            %Neglect data separated by too large of gaps
            biggap = gapsz > p.gapmax;    %Which gaps are 'big'
            bgnd = [find(drng(s,:),1,'first')-1, gnd{s}(biggap)]; %Add start
            bgst = [gst{s}(biggap), find(drng(s,:),1,'last')+1];  %Add end
            if p.startonly
                %Retain only starting segment
                drng(s, bgst(1):end) = false;
            else
                %Consider size of data segments, and retain only largest
                [mx, mi] = max(bgst-bgnd);  %Index of largest segment
                drng(s,:) = false;  drng(s, bgnd(mi)+1:bgst(mi)-1) = true;
            end
            %Neglect gap starts/ends outside of new range
            gst{s}( ~drng(s,gst{s}) ) = []; gnd{s}( ~drng(s,gnd{s}) ) = [];
        end
        
        %Enforce minimum data range
        if nnz(drng(s,:)) < p.dlengthmin; drng(s,:) = false; continue; end
        %Fill gaps with linear interpolation
        for ss = 1:min(numel(gnd{s}),numel(gst{s}))  %FOR each gap
            for sc = 1:nchan
                in{sc}(s, gst{s}(ss)-1:gnd{s}(ss)+1) = linspace( ...
                    in{sc}(s, gst{s}(ss)-1), in{sc}(s, gnd{s}(ss)+1), ...
                    gnd{s}(ss) - gst{s}(ss) + 3 ); %+3 for overlap and each end
            end
        end
    elseif nnz(drng(s,:)) < p.dlengthmin; drng(s,:) = false; continue;  
        %If no gaps, but not a long enough segment, skip
    end  %IF any gaps exist
    
end  %FOR each cell trace

%Retain only useable segments of data with time indices
%   If desired packed as cell array, separate cell tracks and truncate
if p.packascell
    [stord, dst] = find(drng & ~[false(nCells,1), drng(:,1:end-1)]);
    [ndord, dnd] = find(drng & ~[drng(:,2:end), false(nCells,1)]);
    [hasdata, sti] = sort(stord);  [~, ndi] = sort(ndord);
    timei = [dst(sti), dnd(ndi)];   %Time start and end
    
    in = cellfun( @(x)cellfun(@(xx,yy)xx(yy(1):yy(2)), ...
        num2cell(x(hasdata,:),2), num2cell(timei,2), 'UniformOutput', false), ...
        in, 'UniformOutput', false);
    out.time = timei;
else %  IF desired as NaN-padded array, NaNify and clear bad cells
    hasdata = any(drng,2);   out.cellindex = find(hasdata);
    for sc = 1:nchan; in{sc}(~drng) = NaN; in{sc} = in{sc}(hasdata,:); end
end


%Manage output type, structure or cell
if exist('infnames','var')
    %Re-pack data from structure with any non-numeric fields
    repack(~skp) = in;    repack(skp) = inpass;
    %Store as output structure
    if p.packascell;  repack = [repack{:}]; end
    out.data = cell2struct( repack, infnames, 2);
else
    out.data = in;
end

end



%Subfunction: NaN padded array from cell array
function x = sub_nanpad(x)
mx = max(cellfun(@numel, x));
x = cellfun(@(y)[y,nan(1,mx-numel(y))], x, 'UniformOutput', false);
x = cat(1,x{:});
end