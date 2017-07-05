% dataFilt
% Inputs:
%   d - data, output from dataproc. Should be cell array of cell x time matrices
%    d{xy}.data.channel
%   threshval - values for filtering data. 1 or 2 element vector.
%   channel - name of channel to preform filtering on (EX: 'ekar').
%   type -
%           'thresh' - filter out traces containing values above or below
%           values in threshval. Enter threshval as a one element vector
%           containting an upper bound only [UB] , or as a 2 element vector
%           containing upper and lower bounds: [LB,UB]
%           'diff' - filter out traces containing derivatives above the
%           threshold defined in threshval. Derivative calculated using
%           diff.
%           'time' - filter traces for length within a desired time window.
%           Enter start and end times in threshval [start,end]. This filter
%           will keep traces which have more than 50 time points in this
%           window. To change the length filter enter a value for
%           'lengthfilt' in variable arguments.

function [df] = dataFilt(d,threshval,channel,type,varargin)
p.lengthfilt = 50; p.diffparam = 2;
%Input option pair parsing:
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end
%Splits pairs to a structure
for s = 1:2:nin;   p.(lower(varargin{s})) = varargin{s+1};   end
%%
fn = fieldnames(d{1}.data);
% find indices of cells containing data
datacells = 1:numel(d);
datacells = datacells(~cellfun(@isempty,d));

if size(d{1}.data,1) == 1
    switch type
        case 'thresh'
            % make sure threshval is provided
            if isempty(threshval) %#ok<ALIGN>
                error('Must Provide Threshval')
            else end
            % if 2 values provided for threshold filter data max and min
            if numel(threshval) == 2
                for i = datacells
                    idx = sum(d{i}.data.(channel) < threshval(1),2) == 0 & ...
                        sum(d{i}.data.(channel) > threshval(2),2) == 0 ;
                    for ia = 1:numel(fn)
                        df{i,1}.data.(fn{ia}) = d{i}.data.(fn{ia})(idx,:);
                    end
                    if isfield(d{i},'cellindex')
                        df{i,1}.cellindex = d{i}.cellindex(idx,:);
                    end
                end
                % if only one value provided for thresh only filter data max.
            else
                for i = datacells
                    idx = sum(abs(d{i}.data.(channel)) > threshval,2) == 0 ;
                    for ia = 1:numel(fn)
                        df{i,1}.data.(fn{ia}) = d{i}.data.(fn{ia})(idx,:);
                    end
                    if isfield(d{i},'cellindex')
                        df{i,1}.cellindex = d{i}.cellindex(idx,:);
                    end
                end
            end
            
        case 'diff'
            % removed cells whos derivative is above threshval
            for i = datacells
                idx = sum(diff(d{i}.data.(channel),p.diffparam) > threshval,2) == 0;
                if isfield(d{i},'cellindex')
                    df{i,1}.cellindex = d{i}.cellindex(idx,:);
                end
                for ia = 1:numel(fn)
                    df{i,1}.data.(fn{ia}) = d{i}.data.(fn{ia})(idx,:);
                end
            end
            
        case 'time'
            % filters cells which have less than 'lengthfilt' timepoints in
            % time window specified by thresval(1):threshval(2) interval
            
            % if time interval is less than lengthfilt change
            % lengthfilt to length of interval.
            if length(threshval(1):threshval(2)) < p.lengthfilt
                p.lengthfilt = length(threshval(1):threshval(2));
            end
            for i = 1:numel(d)
                idx = sum(~isnan(d{i}.data.(channel)(:,threshval(1):threshval(2))),2) ...
                    >= p.lengthfilt;
                if isfield(d{i},'cellindex')
                    df{i,1}.cellindex = d{i}.cellindex(idx,:);
                end
                for ia = 1:numel(fn)
                    df{i,1}.data.(fn{ia}) = d{i}.data.(fn{ia})(idx,:);
                end
            end
            
        otherwise warning('not an accepted input for type')
    end
    % DEPRECATED DATA TYPE: old data proc / 'packascell' data type
else
    if strcmpi(type,'thresh')
        if isempty(threshval)
            threshval = 5;
        else end
        
        for i = 1:numel(d)
            j=0;
            for ia = 1:size(d{i}.data,1)
                if sum(d{i}.data(ia).(channel) > threshval) == 0
                    j=j+1;
                    df{i,1}.data(j,:) = d{i}.data(ia,:);
                    df{i,1}.time(j,:) = d{i}.time(ia,:);
                    %             df{i,1}.cellindex(j,:) = d{i}.cellindex(ia,:);
                else end
            end
        end
    elseif strcmpi(type,'diff')
        for i = 1:numel(d)
            dm = cell2mat(arrayfun(@(x)max(diff(d{i}.data(x).(channel),2)),...
                1:size(d{i}.data,1),'UniformOutput',false)');
            idx = dm > threshval;
            df{i,1}.data = d{i}.data(~idx,:);
            df{i,1}.time = d{i}.time(~idx,:);
        end
        
    elseif strcmpi(type,'time')
        for i = 1:numel(d)
            m = 0;
            for ia = 1:size(d{i}.data,1)
                trk = d{i}.time(ia,1):d{i}.time(ia,2);
                if sum(trk > threshval(1) & trk < threshval(2)) > p.lengthfilt
                    m=m+1;
                    df{i,1}.data(m,:) = d{i}.data(ia,:);
                    df{i,1}.time(m,:) = d{i}.time(ia,:);
                else
                end
            end
        end
        
    else warning('not an accepted input for type')
    end
end