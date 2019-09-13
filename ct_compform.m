%CT_COMPFORM
%   Convert CellTrace data structures to (and from) a compressed form,
%   allowing easier automated indexing/loop operation across whole
%   datasets.  Compressed form is d{xy} = [nCells x nTime x nChan].
%
%   [d, iv, p] = ct_compform(d)
%       Converts dataset d, and returns iv, an index vector structure with
%   fields 'lin' and 'cel' for linfo and cellindex, respectively.  Output p
%   contains details of the original form, including channel names.
%
%   [d] = ct_compform(d, iv, p)
%       Converts d from compressed form to full form using parameters p and
%   index vectors iv from a previous compression run.
%
%   [d, iv, p, ci] = ct_compform(d, [], [], ci)
%       Converts a channel index, ci, from names to indices, matching the
%   compressed from (i.e. the third dimension of the data arrays).
%   


function [d, iv, p, ci] = ct_compform(d, iv, p, ci)
%Acceptable input formats:
%   d{xy}.data.name = [nC x nT]  (w/ or w/o index fields)
%   d(xy).data.name = [nC x nT]  (w/ or w/o index fields)
%   d{xy} = [nC x nT x nChan] - the compressed format

if ~exist('p','var') || isempty(p)  %IF mapping to compressed form
    %Compress cell array to structure array
    if iscell(d);
        p.gi = ~cellfun('isempty', d);
        d = cat(1,d{:}); p.cell = true; 
    else p.gi = true(size(d));  p.cell = false;
    end
    
    iv = struct('lin',[],'cel',[]);     %Initialize index vectors
    p.cname = fieldnames(d(1).data);    %Collect channel names
    %   Convert channel names to indices as needed
    if exist('ci','var') && iscell(ci) && ~isnumeric(ci{1});
        ci = cellfun(@(cc)find(ismember(p.cname, cc)), ci, 'Un', 0);
    end
    
    %Collect cell and lineage indices
    if isfield(d(1), 'linfo');      iv.lin = cell(size(p.gi));
        iv.lin(p.gi) = {d.linfo};        end
    if isfield(d(1), 'cellindex');  iv.cel = cell(size(p.gi));
        iv.cel(p.gi) = {d.cellindex};    end
    
    %Convert data to compact format
    d = arrayfun(@(x)struct2cell(x.data), d, 'Un', 0);
    d = cellfun(@(x)cat(3, x{:}), d, 'Un', 0);
    %   Include empty cells for bad inidices
    if ~all(p.gi); dt = d; d = cell(size(p.gi)); d(p.gi) = dt; end
    
else  %IF mapping back from compressed form (have p)
    %Ensure iv, civ are present
    if ~exist('iv','var');  iv = struct('lin',[],'cel',[]);  end 
    %Revert to full data structure format
    %   Convert to named channels
    dtmp = cellfun( @(x)cell2struct(num2cell(x,[1,2]), p.cname, 3), ...
        d(p.gi), 'Un', 0 );  d = struct('data',[]);
    %   Add data
    [d(1:numel(dtmp)).data] = deal(dtmp{:});
    %   Repackage cell and lineage indices
    if ~isempty(iv.cel); [d.cellindex] = deal(iv.cel{p.gi});   end
    if ~isempty(iv.lin); [d.linfo] = deal(iv.lin{p.gi});       end
    %   Re-encapsulate in cell array as needed
    if p.cell; dt = d; d = cell(size(p.gi)); d(p.gi) = num2cell(dt); end    
end

end