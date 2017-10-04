function [] = plotby_singlecells(plotby,data,channel,pmd,varargin)
%% Check input
plotby = lower(plotby);

if any(strcmp(plotby,{'treatment','treatments','tx','txs'}))
    plotby = 'treatment';
elseif any(strcmp(plotby,{'cell','celltype','celltypes'}))
    plotby = 'celltype';
else error('Accepted plotby types: treatment or celltype. plotby must be a string')
end
%% Input parsing
p.lines = [];
%Input option pair parsing:
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end
%Splits pairs to a structure
for s = 1:2:nin;   p.(lower(varargin{s})) = varargin{s+1};   end
%% get cell type names
c1 = pmd.Cell(:,:,1);
c2 = c1(cellfun(@isstr,c1)); celltypes = unique(c2); % find unique names of celltypes
celltypes = celltypes(~cellfun(@isempty,celltypes)); % discard empty name fields
% get rid of @ density if present for fieldname
cellfn = cellfun(@(x)x{1},regexp(celltypes,'@','split'),'un',0);
% remove any spaces in name
cellfn = arrayfun(@(x)strrep(cellfn{x},' ',''),1:size(cellfn),'un',0);
% make list of xys for each celltype
for s = 1:numel(cellfn)
    %   Assemble list of xys with the current cname
    idx.Cell.(cellfn{s}) = [pmd.xy{strcmp(celltypes{s},c1)}];
end

%% get treatments and find unique combos
catTxnames = {'pTx','Tx1','Tx2','Tx3'}; % possible treatment fields
[~,I]= intersect(catTxnames, fieldnames(pmd)); % tx fields in experiment
I = sort(I); catTxnames = catTxnames(I);

Txcat = cell(size(pmd.Cell,1),size(pmd.Cell,2)); % initialize cell array
linetp = cell(size(pmd.Cell,1),size(pmd.Cell,2)); % initialize cell array

for s = 1:numel(catTxnames) % cat each treatment in a loop
    Txcat = cat(3,Txcat, cellfun(@(x)cat(2,x{:}),...
        cat(3,num2cell(pmd.(catTxnames{s})(:,:,:),3)),'Un',0));
    linetp = cat(3,linetp, cellfun(@(x)cat(2,x{:}),...
        cat(3,num2cell(pmd.(catTxnames{s})(:,:,5),3)),'Un',0));
end
linetp = linetp(:,:,2:end); % remove empty layer from initialization
linetp = num2cell(linetp,3);
% color sequence for treatment lines
linecolor = parula(numel(catTxnames)+1);


stridx = cellfun(@isstr,Txcat);
spacer =  repmat({'+'},size(Txcat)); % add + between txs for readability
spacer(~stridx) = {[]};

% add + between treatment names
Txcat2 = Txcat(:,:,2);
for s = 2:size(Txcat,3)-1
    Txcat2 = cat(3,Txcat2,spacer(:,:,s+1),Txcat(:,:,s+1));
end
% cat Tx name and space together into one long string
Txcat2 = cellfun(@(x)cat(2,x{:}),cat(3,num2cell(Txcat2,3)),'Un',0);
stridx = cellfun(@isstr,Txcat2);
% find the unique treatment combos
txs = unique(Txcat2(cellfun(@isstr,Txcat2)));
% remove extra spaces in name caused by NaNs to make pretty names
treatments = arrayfun(@(x)regexprep(txs{x},char(0),''),1:size(txs),'un',0)';
% remove double +s
treatments = arrayfun(@(x)regexprep(treatments{x},'++','+'),1:size(treatments),'un',0)';
% remove @##tp
treatments = arrayfun(@(x)regexprep(treatments{x},'\d*tp',''),1:size(treatments),'un',0)';

% vector of xys for indexing during plotter
xymat = cell2mat(pmd.xy(stridx));
linetp = linetp(stridx);
%% Plotting
switch plotby
    case 'celltype'
        [smin,smax] = sminmax(data,channel);
        [nrow,ncol,pos] = subplot_sizer(numel(txs)); % optimal subplot dimensions
        for s = 1:numel(cellfn)
            figure,
            for sa = 1:numel(txs) % new subplot for each treatment
                ah = subplot(nrow,ncol,sa); hold on
                xy = xymat(strcmp(txs{sa},Txcat2(stridx))); % find xy
                xy = intersect(xy, idx.Cell.(cellfn{s})); % find xy for treatment/celltype
                xy = xy(1);
                % actual plotting
                stackplot(ah,data{xy}.data.(channel), smin, smax)
                title(treatments{sa})
                
                if p.lines % if lines is on
                    yy = ylim; % get ylims of plot for line
                    for sb = 1:size(linetp{xymat == xy},3)
                        xx = str2double(linetp{xymat == xy}{sb});
                        if ~isnan(xx)
                            ob{sb} = line([xx,xx],[yy(1),yy(2)],'Color',...
                                [linecolor(sb,:),0.7]);
                        else
                        end
                    end
                end
            end
            if p.lines
            % take only txt for which lines exist
            legtxt = regexp(txs{sa},'+','split');
            if exist('ob','var') && sum(~cellfun(@isempty,ob)) ~= numel(legtxt)
                % split tx names for legend if lines is on
                tps = linetp{xymat == xy}; tps = tps(~cellfun(@isempty,tps));
                txidx = zeros(1,numel(legtxt));
                for se = 1:numel(tps) % find strings associated with the correct time points
                    txidx = txidx + ~cellfun(@isempty,strfind(legtxt,tps{se}));
                end
                txidx = logical(txidx); % make logical index
                legtxt = legtxt(txidx); % trim legend names to string assoc. w/ tp #
                % remove @##tp
            
            legtxt  = arrayfun(@(x)regexprep(legtxt{x},'\d*tp',''),1:numel(legtxt),'un',0)'; % remove tp for label
            ob = ob(~cellfun(@isempty,ob)); % remove empties from line objects
            legend([ob{:}],legtxt,'Orientation','vertical','Location','BestOutside') % make legend
            clear ob legtxt
            end
            end
            suptitle(cellfn{s})
            set(gcf,'Position',pos)
        end
        
    case 'treatment'
        [smin,smax] = sminmax(data,channel);
        [nrow,ncol,pos] = subplot_sizer(numel(cellfn)); % optimal subplot dimensions
        for s = 1:numel(txs)
            figure,
            for sa = 1:numel(cellfn) % new subplot for each treatment
                ah = subplot(nrow,ncol,sa); hold on
                xy = xymat(strcmp(txs{s},Txcat2(stridx))); % find xy
                xy = intersect(xy, idx.Cell.(cellfn{sa})); % find xy for treatment/celltype
                xy = xy(1); % take only 1 xy for plotting
                % actual plotting
                stackplot(ah,data{xy}.data.(channel), smin, smax)
                title(cellfn{sa})
                
                if p.lines 
                    yy = ylim; % get ylims of plot for line
                    for sb = 1:size(linetp{xymat == xy},3)
                        xx = str2double(linetp{xymat == xy}{sb});
                        if ~isnan(xx)
                            ob{sb} = line([xx,xx],[yy(1),yy(2)],'Color',...
                                [linecolor(sb,:),0.7]);
                        else
                        end
                    end
                end
            end
            if p.lines
            % take only txt for which lines exist
            legtxt = regexp(txs{s},'+','split');
            if exist('ob','var') && sum(~cellfun(@isempty,ob)) ~= numel(legtxt)
                % split tx names for legend if lines is on
                tps = linetp{xymat == xy}; tps = tps(~cellfun(@isempty,tps));
                txidx = zeros(1,numel(legtxt));
                for se = 1:numel(tps) % find strings associated with the correct time points
                    txidx = txidx + ~cellfun(@isempty,strfind(legtxt,tps{se}));
                end
                txidx = logical(txidx); % make logical index
                legtxt = legtxt(txidx); % trim legend names to string assoc. w/ tp #
                % remove @##tp
            
            legtxt  = arrayfun(@(x)regexprep(legtxt{x},'\d*tp',''),1:numel(legtxt),'un',0)'; % remove tp for label
            ob = ob(~cellfun(@isempty,ob)); % remove empties from line objects
            legend([ob{:}],legtxt,'Orientation','vertical','Location','BestOutside') % make legend
            clear ob legtxt
            end
            end
            suptitle(treatments{s})
            set(gcf,'Position',pos)
        end
end
end
function [] = stackplot(ah,dat, smin, smax)
nTime = size(dat,2);
offst = 0;
nTracks = 5;
datlength = sum(~isnan(dat),2);
[~,I] = sort(datlength,'descend');
dat = dat(I,:);
for i = 1:nTracks
    plot(ah,1:nTime,dat(i,:)+ offst,'LineWidth',1,'Color','k'); hold on;
    line([1,.05*nTime],[offst+smax, offst+smax], 'Color','k','LineWidth',1);
    offst = offst + (smax - smin);
end
axis([1,nTime,smin*0.95,smin+(smax-smin)*nTracks*1.05]);
set(gca,'YTick',[smin,smax],'YTicklabel',...
    {num2str(smin,'%3.1f'),num2str(smax,'%3.1f')});
box on;
end
%% find min and max of all data to set scales
function [smin,smax] = sminmax(data,channel)
mx = []; mn = [];
for s = 1:numel(data)
    mx = [mx;quantile(data{s}.data.(channel)(:),.95)];
    mn = [mn;quantile(data{s}.data.(channel)(:),.05)];
end
smin = min(mn); smax = max(mx);
end
%% function to select figure dimensions
function [nrow,ncol,pos] = subplot_sizer(numplots)
% determines optimal subplot sizing
if numplots > 6
    ncol = ceil(sqrt(numplots));
    nrow = ceil(numplots/ncol);
else
    ncol = numplots; nrow = 1;
end
hp = ncol*450; if hp > 1920, hp = 1920; end
vp = nrow * 450; if vp > 964, vp = 964; end
pos = [10,550,hp,vp];
end
