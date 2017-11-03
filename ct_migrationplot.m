%CT_MIGRATIONPLOT
%   Plots X and Y coordinate data over time to visualize migration.
%
%   ct_migrationplot(X,Y,...)
%       plots the cell traces in space from the x and y coordinates (X,Y).
%       Additional parameters may be provided as Name/Value pairs.
%
%   ct_migrationplot(AH,X,Y,...)
%       plots the traces on axes AH
%
%   Parameters:
%   noflippy    - Logical, TRUE to prevent inversion of Y axis as is
%                   typical for MATLAB image plots.
%   framesize   - 2 element array, specifying number of pixels in the
%                   original frame [X, Y]
%   linewidth   - Linewidth at which to plot traces. Default is 1.5\
%   colorbar    - Logical, TRUE (default) to show a colorbar
%   tsamp       - Sampling interval (time between samples)
%   tunit       - Units for sampling interval (e.g. 's', 'min', 'hr')


function [] = ct_migrationplot(varargin)
%% Input handling
p.noflipy = false;
p.linewidth = 1.5;
p.colorbar = true;
p.tsamp = 1;
p.tunit = 'au';

%Get axis handle, if provided
if ishandle(varargin{1}); ah = varargin{1}; varargin = varargin(2:end); 
else figure; ah = axes; end

%Re-assign x and y from input
[x,y] = deal(varargin{1:2});%Assume provided X, then Y
%Ensure X and Y sizes match exactly
assert(all(size(x)==size(y)), 'CT:Size', 'X and Y must be the same size.');

%Default frame size, based on trace locations
p.framesize = 10.*ceil([max(x(:)), max(y(:))]./10);

%Input option pair parsing:
varargin(1:2) = [];  nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end
%Split pairs to the parameter structure
for s = 1:2:nin;    p.(lower(varargin{s})) = varargin{s+1};     end


%% Preliminaries
%Flip Y, unless directed not to
if ~p.noflipy;      y = p.framesize(2) - y;     end

%Get size of data
[nc,nt] = size(x);
%Set up relative time vector
tm = (0:nt-1)./(nt-1);

%Set up color mapping
cmp = interp1([0,0.25,0.75,1], ...
    [0,0,0.75; 0.75,0,0.75; 1,0,0; 1,1,0], (0:63)/63);


%% Plotting
%   Plots gradient colored lines by relative time

zt = zeros(2,nt);    %Pre-fill uniform Z value
%   Initialize figure, with black background
axes(ah);  hold on;  set(ah, 'Color', 'k');
%   Loop through cell traces, need to make a surface for each
for s = 1:nc
    surface([x(s,:);x(s,:)], [y(s,:);y(s,:)], zt, [tm;tm],...
        'EdgeColor','interp', 'LineWidth', p.linewidth, ...
        'MeshStyle','row', 'FaceColor','none');
end
%   Set axis bounds etc.
colormap(ah, cmp);  axis([0, p.framesize(1), 0, p.framesize(2)]);
set(ah, 'YDir', 'Reverse', 'XTick', [], 'YTick', []);

if p.colorbar
    %   Set colormap labeling (first and last time value only)
    cmt = [0, 1];   cml = {'0', num2str(nt*p.tsamp)};
    %   Insert colorbar and label
    cbh = colorbar('peer', ah, 'Ticks', cmt, 'TickLabels', cml);
    cbh.Label.String = ['Time (',p.tunit,')'];
end

end




%   Keeping these plotting bits for any future use
% %Colored dots along the line
% figure; ah = axes;
% plot(ah, x',y', 'w'); hold on;
% scatter(ah, x(:), y(:), 2, tm(:), 'filled');
% colormap(ah, 'hot');
% set(ah, 'Color', 'k');
% axis([0, p.framesize(1), 0, p.framesize(2)]);
% set(ah, 'YDir', 'Reverse'); caxis(ah, [0,1]);
% 
% 
% %Dot at end of line
% figure; ah = axes;  hold on; cm = colormap('winter'); cs = size(cm,1);
% for s = 1:nc
%     ci = mod(s-1,cs)+1;
%     plot(ah, x(s,:), y(s,:), 'Color', cm(ci,:));
%     scatter(ah, x(s,end), y(s,end), 10, cm(ci,:), 'filled', 'Marker', 'o');
% end
% % colormap(ah, 'hsv');
% axis([0, p.framesize(1), 0, p.framesize(2)]);
% set(ah, 'YDir', 'Reverse');