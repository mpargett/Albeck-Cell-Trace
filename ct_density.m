%CT_DENSITY
%   Cell density analysis for Cell Trace data.
%
%Usage:
%   D = ct_density(X, Y)
%       returns D, a structure containing density estimates, as detailed in
%       the 'Output' section below.    
%   D = ct_density(X, Y, 'PARAM1', VALUE, ...) uses the specified parameter
%       values (see 'Parameters' section below). 
%
%Inputs:
%   X - Abscissa (X) coordinates, as a column vector
%   Y - Ordinate (Y) coordinates, as a column vector
%
%Output:  
%   D.avg   scalar average density (i.e. total cells/total area)
%    .pdf   probability density function values at each cell site
%    .dloc  local density based on gaussian kernel
%    .dvor  local density based on Voronoi tesselation
%    .emin  minimum edge length per cell, based on Delaunay triangulation
%    .emed  median edge length per cell
%    .evar  variance in edge length per cell
%
%Parameters:
%   SCALE - 1x2 double providing scaling factors for coordinate values
%

function d = ct_density(x, y, varargin)
%Define default parameters
p.scale = [1,1];

%Parse inputs
p = ct_input(varargin, p);

%Scale data
xy = bsxfun(@times, [x(:),y(:)], p.scale);

nc = size(xy,1);  %Number of cells
%Get region bounds from data
rg = [min(xy,[],1); max(xy,[],1)];      

%Average density, using incoming data units
d.avg = nc/prod(rg(2,:)-rg(1,:));


%% Local kernel density
%Evalulate kernal density, with gaussian and automatic bandwidth
[lkd, ~, u] = ksdensity(xy,xy);
%Store the PDF of kernel-based density 
d.pdf = lkd;    

%Estimate the uniform spacing density matched at each node
%   Define approximate model mapping kernel density value to local density
dmap = cfit(fittype(@(a,b,c,x)((x-c)./a).^(1./b)), 2.5, -2, -0.077);
%       Model was fit with code at bottom of this file (commented)

%Evaluate sum of point contributions to node
%   Reverse scaling used in kernel density estimate
%   i.e. bandwidth area, weights, 2D normal kernel scaling
%   The forward process is: Get x and y distances, Scale each by bandwidth
%   (dx/ux, dy/uy), Evaluate kernel on each dimension, Take product across
%   dimensions, Multiply by weights (which are divided by sum(weights)),
%   Scale by bandwidth area (divide by prod(u)) 
consum = lkd*prod(u)*nc*sqrt(2*pi);
%Evaluate the density as 1/area around each point in the uniformly spaced
%   approximation, scaled by the bandwidth used in each dimension
d.dloc = 1./(dmap(consum).^2.*prod(u));


%% Triangulation based metrics
%Get triangulation
dt = delaunayTriangulation(xy);     

% ------------------------------------------------------------------------
%Estimate local density by specific area around each cell (1/density)
% ------------------------------------------------------------------------
%Get Voronoi diagram
[v,r] = voronoiDiagram(dt);
%Validate regions by completeness in frame
%   Determine bad vertices as being out of region bounds
bv = find(~all(bsxfun(@lt, v, rg(2,:)) & bsxfun(@gt, v, rg(1,:)), 2));
%   Determine bad regions as having a bad vertex, and filter them out
br = cellfun(@(rr)any(ismember(bv,rr)), r);     r = r(~br);

%   Define function to take line integrals around a region
lint = @(idx)abs(sum( (v(idx,2) + v(circshift(idx,[0,-1]),2))...
            .*diff(v([idx,idx(1)],1)) ))/2;
%   Evaluate line integrals around each valid region (for specific area)
%       and take reciprocal for local density estimate
d.dvor = 1./cellfun(lint, r);

% ------------------------------------------------------------------------
%Estimate neighbor cell-cell distances by Delaunay edge lengths
% ------------------------------------------------------------------------
%Get edge lengths and edges per cell index
ed = edges(dt);     %Collect edges
%   Get edge lengths
edl = sqrt(sum((dt.Points(ed(:,2),:) - dt.Points(ed(:,1),:)).^2,2));
%   Get edge indices per cell, use like: edl(epc(:,1))
epci = squeeze(any(bsxfun(@eq, ed, shiftdim(1:nc,-1)),2));
%   Get edge lengths per cell
epc = arrayfun(@(x)edl(epci(:,x)), 1:nc, 'UniformOutput', 0);

%   Report edge length statistics per cell
d.emin = cellfun(@min, epc);
d.emed = cellfun(@median, epc);
d.evar = cellfun(@var, epc);


%FIXME - any additions?
%Evaluate any features of region shape?
%Estimate local density variance? (say, by values for nearby cells)

end


%NOTES
% %To fit approximately power law relationship between uniform grid spacing
% %   and normal kernal density value:
% %Define a reference spacing
% dxx = 0.1;     ss = dxx:dxx:6;  clear yout;
% for s = 1:numel(ss)
% dx = ss(s);	 n = max(6,ceil(6./dx));	%Spacing and number of nodes
% A = ones(n,1)*((0:n-1)-(n-1)/2);        %Node array (whole space)
% rspc = ( exp(-0.5 * (A.^2 + A'.^2).*dx.^2) ./ sqrt(2*pi) ); %Kernel vals
% yout(s) = sum(rspc(:));  %Summed kernal density
% end
% %   Compare against density (1/area-per-point)
% xa = ss;
% %   Fit power law model
% w = ctm_infer(xa,yout,'runinfo',false,'models', 'power2');
% figure; scatter(xa, yout); hold on; plot(xa, w.m(2).model(xa));
% w.m(2).model    %Display fit model with coefficients
% %   Invert model for reverse power law mapping (show coefficients)
% c.b = 1./(-1.997);  c.a = (1./(2.513)).^c.b;  c.c = (0.07656.*c.a).^c.b
% dmap = cfit(fittype('power2'), c.a, c.b, c.c);   %Density model
% dmap = cfit(fittype(@(a,b,c,x)((x-c)./a).^(1./b)), ...
%   w.m(2).model.a, w.m(2).model.b, w.m(2).model.c);
% plot(dmap(yout), yout); %Test against previously plotted data


