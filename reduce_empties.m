%REDUCE_EMPTIES
%   Convert empty contents of any size to 0x0 for consistency.
%   Operates recursively through cells and structures.
%
%   y = reduce_empties(x)
%
%   returns, y, the same datatype as x, with any empty contents reduced to
%   0x0.
%

function x = reduce_empties(x)
%Handle container types (cell, structure) recursively
if iscell(x)
    %   Find any empty contents, replace with 0x0 empty
    z = cellfun(@isempty, x);    [x{z}] = deal([]);
    %   Recursively process non-empty cells
    x(~z) = cellfun(@reduce_empties, x(~z), 'un', 0);
    
elseif isstruct(x)
    fn = fieldnames(x);
    %   Find any empty contents, replace with 0x0 empty
    for ss = 1:numel(x)
        z = structfun(@isempty, x(ss));
        for s = find(z)'; x(ss).(fn{s}) = []; end
        %   Recursively process non-empty fields
        for s = find(~z)'
            x(ss).(fn{s}) = reduce_empties(x(ss).(fn{s}));
        end
    end
    
elseif isempty(x)  %Handle in case x itself is empty
    if ischar(x);  x = char([]); else x = [];  end
end

end



