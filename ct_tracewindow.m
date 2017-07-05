%CT_TRACEWINDOW
%   Restrict Cell Trace data to specified time windows
%
%   d = ct_tracewindow(d, tref, tspan, tbnd, outformat)
%
%   tref: reference time point (integer value, 'start', or 'end')
%   tspan: distance from tref to keep (in number of time points, negative
%           for time before reference)
%   tbnd: bounds on valid time [start, end], will not consider traces past
%           outside of these times
%   outformat: format of output data ('cell' or 'num')


function d = ct_tracewindow(d, tref, tspan, tbnd, outformat)

%Parse input data type
if iscell(d);                                       tp = 'cell';
elseif isstruct(d);     d = mat2cell(d);            tp = 'struct';
elseif isnumeric(d);    d = {struct('name', d)};    tp = 'num';
end

%Default time bound
if ~exist('tbnd','var') || isempty(tbnd); tbnd = [-Inf, Inf]; end

%Override output format if requested
if exist('outformat','var'); ovr = true;
    switch lower(outformat); 
        case 'cell'; usenp = false;
        case {'num', 'array'}; usenp = true;
        otherwise; error('Output format must be "cell" or "num".');
    end
else ovr = false;
end

nx = numel(d);  wasnp = false(1,nx);
for sx = 1:nx
    %Get start and end times
    if isfield(d{sx}, 'time'); st = d{sx}.time(:,1); et = d{sx}.time(:,2);
        wasnp(sx) = false;
        %Repack to NaN-padded array
        nc = length(d{sx}.data);    
        fn = fieldnames(d{sx}.data(1));
        for sf = 1:numel(fn);   v.(fn{sf}) = NaN(nc, max(et));
            for s = 1:nc
                v.(fn{sf})(s, st(s):et(s)) = d{sx}.data(s).(fn{sf});
            end
        end
    else    wasnp(sx) = true;
        %  Find start and end of each trace from NaN-pads
        v = d{sx};  fn = fieldnames(v);    gi = ~isnan(v.(fn{1}));
        [a,b] = find ( [gi(:,1), gi(:,2:end) & ~gi(:,1:end-1)] ); 	
            [a,ai] = sort(a, 'ascend');  st = b(ai);
        [a,b] = find ( [gi(:,1:end-1) & ~gi(:,2:end), gi(:,end)] );	
            [a,ai] = sort(a, 'ascend');  et = b(ai);
    end
    nfn = numel(fn);
    
    %Get reference times, ensuring validity
    if ischar(tref)
        switch lower(tref)
            case 'end';     rt = et;
            case 'start';   rt = st;
        end
    else
        rt = ones(size(st))*tref;
        %Check validity of reference times
        isv = rt >= max(st,tbnd(1)) & rt <= min(et,tbnd(2)) & ~isnan(rt);
        %Remove any invalid traces
        rt = rt(isv); st = st(isv); et = et(isv);
        v = structfun(@(x)x(isv,:), v, 'UniformOutput', false);
    end
    
    %Determine window end
    we = rt + tspan;
    %Check validity of window ends
    isv = we >= max(st,tbnd(1)) & we <= min(et,tbnd(2)) & ~isnan(we);
    %Excise and return valid data
    rt = rt(isv); we = we(isv); st = st(isv); et = et(isv);
    v = structfun(@(x)x(isv,:), v, 'UniformOutput', false);
    %   Window excision function
    if tspan > 0;   exfun = @(x,i1,i2)x(i1:i2);
    else            exfun = @(x,i1,i2)x(i2:i1);   end
    
    %   Excise valid data
    nvc = nnz(isv);  vout = nan(nvc, tspan+1);
    for sf = 1:nfn
        for s = 1:nvc
            vout.(fn{sf})(s,:) = exfun(v.(fn{sf})(s,:), rt(s), we(s));
        end
    end
    
    %Repack to input data format
    if (ovr && ~usenp) || (~ovr && ~wasnp(sx))
        dtmp.time = [st,et]; dtmp.data = cell2struct(cell(nvc,nfn), fn, 2);
        for sf = 1:nfn
            stmp = num2cell(vout.(fn{sf}),2);
            [dtmp.data.(fn{sf})] = deal(stmp{:});
        end
            d{sx} = dtmp;
    else    d{sx} = vout;
    end
    
end

%Repack to input layout
switch tp
    case 'cell';    %Do nothing
    case 'struct';  d = cell2mat(d);
    case 'num';     d = d{1}.name;
end

end