function [h,ax,BigAx,hhist,pax] = plotmatrix(varargin)
%PLOTMATRIX Scatter plot matrix.
%   PLOTMATRIX(X,Y) scatter plots the columns of X against the columns
%   of Y.  If X is P-by-M and Y is P-by-N, PLOTMATRIX will produce a
%   N-by-M matrix of axes. PLOTMATRIX(Y) is the same as PLOTMATRIX(Y,Y)
%   except that the diagonal will be replaced by HISTOGRAM(Y(:,i)).
%
%   PLOTMATRIX(...,'LineSpec') uses the given line specification in the
%   string 'LineSpec'; '.' is the default (see PLOT for possibilities).
%
%   PLOTMATRIX(AX,...) uses AX as the BigAx instead of GCA.
%
%   [H,AX,BigAx,P,PAx] = PLOTMATRIX(...) returns a matrix of handles
%   to the objects created in H, a matrix of handles to the individual
%   subaxes in AX, a handle to big (invisible) axes that frame the
%   subaxes in BigAx, a vector of handles for the histogram plots in
%   P, and a vector of handles for invisible axes that control the
%   histogram axes scales in PAx.  BigAx is left as the CurrentAxes so
%   that a subsequent TITLE, XLABEL, or YLABEL will be centered with
%   respect to the matrix of axes.
%
%   Example:
%       x = randn(50,3); y = x*[-1 2 1;2 0 1;1 -2 3;]';
%       plotmatrix(y)

%   Copyright 1984-2020 The MathWorks, Inc.

% Parse possible Axes input
[cax,args,nargs] = axescheck(varargin{:});
if nargs < 1
    error(message('MATLAB:narginchk:notEnoughInputs'));
elseif nargs > 3
    error(message('MATLAB:narginchk:tooManyInputs'));
end
nin = nargs;

sym = '.'; % Default scatter plot symbol.
dohist = 0;

if matlab.graphics.internal.isCharOrString(args{nin})
    args{nin} = char(args{nin});
    sym = args{nin};
    [~,~,~,msg] = colstyle(sym);
    if ~isempty(msg), error(msg); end
    nin = nin - 1;
end

if nin==2 % plotmatrix(y)
    rows = size(args{1},2); cols = rows;
    x = args{1}; y = args{1};
    dohist = 1;
    
elseif nin==3 % plotmatrix(x,y)
    rows = size(args{2},2); cols = size(args{1},2);
    x = args{1}; y = args{2};
else
    %error(message('MATLAB:plotmatrix:InvalidLineSpec'));
end

metric_names = args{numel(args)};

colors = linspecer(round(size(x,2)/2));
    colors = [colors;colors];
x = matlab.graphics.chart.internal.datachk(x,'numeric');
y = matlab.graphics.chart.internal.datachk(y,'numeric');

% Don't plot anything if either x or y is empty
hhist = gobjects(0);
pax = gobjects(0);
if isempty(rows) || isempty(cols)
    if nargout>0, h = gobjects(0); ax = gobjects(0); BigAx = gobjects(0); end
    return
end

if ~ismatrix(x) || ~ismatrix(y)
    error(message('MATLAB:plotmatrix:InvalidXYMatrices'))
end
if size(x,1)~=size(y,1) || size(x,3)~=size(y,3)
    error(message('MATLAB:plotmatrix:XYSizeMismatch'));
end

if isempty(cax)
    parent = handle(get(groot, 'CurrentFigure'));
    if ~isempty(parent) && ~isempty(parent.CurrentAxes)
        parent = parent.CurrentAxes.Parent;
    end
else
    parent = cax.Parent;
end

inLayout=isa(parent,'matlab.graphics.layout.Layout');
if inLayout && isempty(cax)
    cax=gca;
end

% Error if AutoResizeChildren is 'on'
if ~isempty(parent) && isprop(parent,'AutoResizeChildren') && strcmp(parent.AutoResizeChildren,'on')
    error(message('MATLAB:plotmatrix:AutoResizeChildren'))
end

% Create/find BigAx and make it invisible
BigAx = newplot(cax);
fig = ancestor(BigAx,'figure');
hold_state = ishold(BigAx);
set(BigAx,'Visible','off','color','none')
disableDefaultInteractivity(BigAx);

% Customize the AxesToolbar
% Call the internal function to silently work in web graphics
[~, btns] = axtoolbar(BigAx, {'pan', 'zoomin', 'zoomout', 'datacursor', 'brush'});
set(btns, 'Visible', 'on');

% Getting the position will force an auto-calc and update the axes if
% necessary. This should be done before calling GetLayoutInformation.
BigAxesPos = get(BigAx,'InnerPosition');
BigAxUnits = get(BigAx,'Units');

% Create and plot into axes
ax = gobjects(rows,cols);
width = BigAxesPos(3)/cols;
height = BigAxesPos(4)/rows;
space = .02; % 2 percent space between axes
pos = BigAxesPos + [space*[width height] 0 0];
m = size(y,1);
k = size(y,3);
xlim = zeros([rows cols 2]);
ylim = zeros([rows cols 2]);
BigAxHV = get(BigAx,'HandleVisibility');
BigAxParent = get(BigAx,'Parent');
if inLayout
    % In TiledChartLayouts BigAx will be parented to the layout, but the
    % visible axes (which use BigAxParent below) will be parented to the
    % container holding the layout.
    BigAxParent=ancestor(BigAx, 'matlab.ui.internal.mixin.CanvasHostMixin','node');
end

% Use the pixel width/height of the axes to determine the marker size.
if any(sym=='.')
    pixelPos = BigAx.GetLayoutInformation.Position;
    markersize = max(1,min(15,round(15*min(pixelPos(3:4))/max(1,size(x,1))/max(rows,cols))));
else
    markersize = get(0,'DefaultLineMarkerSize');
end

paxes = findobj(fig,'Type','axes','tag','PlotMatrixScatterAx');

for i=rows:-1:1
    for j=cols:-1:1
        axPos = [pos(1)+(j-1)*width pos(2)+(rows-i)*height ...
            width*(1-space) height*(1-space)];
        findax = findaxpos(paxes, axPos);
        if isempty(findax)
            ax(i,j) = axes('Units',BigAxUnits,'InnerPosition',axPos,'HandleVisibility',BigAxHV,'parent',BigAxParent, 'Toolbar', []);
            set(ax(i,j),'visible','on');
        else
            ax(i,j) = findax(1);
        end
        
        hh(i,j,:) = plot(reshape(x(:,j,:),[m k]), ...
            reshape(y(:,i,:),[m k]),sym,'parent',ax(i,j),'Color',colors(i,:))';
                lsline
                
                [R,P] = corr(x(:,i,:),y(:,j,:));
                %text(0,0,num2str(R));
                text(0,.9,num2str(R),'Units','normalized')
                
                if P<=.05 & i~=j
               % ,'Color','r')
                set(ax(i,j),'XColor','r', 'YColor','r')
                set(ax(i,j),'LineWidth',2)
                end
                
                if P<=.1 & P>.05 & i~=j
               % ,'Color','r')
                set(ax(i,j),'XColor','y', 'YColor','y')
                set(ax(i,j),'LineWidth',2)
                end
                
                

        set(hh(i,j,:),'markersize',markersize);
        set(ax(i,j),'xlimmode','auto','ylimmode','auto','xgrid','off','ygrid','off')
        xlim(i,j,:) = get(ax(i,j),'xlim');
        ylim(i,j,:) = get(ax(i,j),'ylim');
        
        if j == 1
                    set(get(ax(i,j),'Ylabel'),'String',metric_names{i})
                    yticklabels([])
                     hYLabel = get(ax(i,j),'YLabel');
 set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment', 'right','Color','k')
        end
                
        if i == rows
                    set(get(ax(i,j),'Xlabel'),'String',metric_names{j})
                    xticklabels([])
                     hXLabel = get(ax(i,j),'XLabel');
 set(hXLabel,'rotation',0,'Color','k');%,'VerticalAlignment','middle','HorizontalAlignment', 'right')
                end
        
        % Disable AxesToolbar
        currAx = ax(i,j);
        currAx.Toolbar = [];
    end
end

xlimmin = min(xlim(:,:,1),[],1); xlimmax = max(xlim(:,:,2),[],1);
ylimmin = min(ylim(:,:,1),[],2); ylimmax = max(ylim(:,:,2),[],2);
lsline
% Try to be smart about axes limits and labels.  Set all the limits of a
% row or column to be the same and inset the tick marks by 10 percent.
inset = .15;
for i=1:rows
    set(ax(i,1),'ylim',[ylimmin(i,1) ylimmax(i,1)])
    dy = diff(get(ax(i,1),'ylim'))*inset;
    set(ax(i,:),'ylim',[ylimmin(i,1)-dy ylimmax(i,1)+dy])
end
dx = zeros(1,cols);
for j=1:cols
    set(ax(1,j),'xlim',[xlimmin(1,j) xlimmax(1,j)])
    dx(j) = diff(get(ax(1,j),'xlim'))*inset;
    set(ax(:,j),'xlim',[xlimmin(1,j)-dx(j) xlimmax(1,j)+dx(j)])
end

set(ax(1:rows-1,:),'xticklabel','')
set(ax(:,2:cols),'yticklabel','')
set(BigAx,'XTick',get(ax(rows,1),'xtick'),'YTick',get(ax(rows,1),'ytick'), ...
    'YLim',get(ax(rows,1),'ylim'),... %help Axes make room for y-label
    'tag','PlotMatrixBigAx')
set(ax,'tag','PlotMatrixScatterAx');

if dohist % Put a histogram on the diagonal for plotmatrix(y) case
    histAxes = gobjects(1, rows);
    paxes = findobj(fig,'Type','axes','tag','PlotMatrixHistAx');
    pax = gobjects(1, rows);
    for i=rows:-1:1
        axPos = get(ax(i,i),'InnerPosition');
        findax = findaxpos(paxes, axPos);
        if isempty(findax)
            axUnits = get(ax(i,i),'Units');
            histax = axes('Units',axUnits,'InnerPosition',axPos,'HandleVisibility',BigAxHV,'parent',BigAxParent, 'Toolbar', []);
            set(histax,'visible','on');
            histAxes(i) = histax;
        else
            histax = findax(1);
        end
        hhist(i) = histogram(histax,y(:,i,:),[min(y(:,i,:)):.1*range(y(:,i,:)):max((y(:,i,:)))],'FaceColor',colors(i,:),'EdgeColor','none','FaceAlpha',.75);
        set(histax,'xtick',[],'ytick',[],'xgrid','off','ygrid','off');
        set(histax,'xlim',[xlimmin(1,i)-dx(i) xlimmax(1,i)+dx(i)])
        set(histax,'tag','PlotMatrixHistAx');
        % Disable the AxesToolbar
        histax.Toolbar = [];
        pax(i) = histax;  % ax handles for histograms
    end
else
    histAxes = []; % Make empty for listener
end

% Set Title and X/YLabel visibility to on and strings to empty
set([BigAx.Title_I, BigAx.XAxis.Label, BigAx.YAxis.Label], ...
    'String','','Visible','on')

BigAx.UserData = {ax, histAxes, BigAxesPos};
addlistener(BigAx, 'MarkedClean', @(~,~) matlab.graphics.internal.resizePlotMatrix(BigAx, rows, cols));

% Make BigAx the CurrentAxes and set up cla behavior
set(fig,'CurrentAx',BigAx)
addlistener(BigAx, 'Cla', @(~,~) deleteSubAxes(ax, histAxes));


if ~hold_state && ~inLayout 
    try %#ok<TRYNC>
        set(fig,'NextPlot','replacechildren')
    end
end

if nargout~=0
    h = hh;
end

end

function deleteSubAxes(gridOfAxes, histogramAxes)
delete(gridOfAxes);
delete(histogramAxes);
end


function findax = findaxpos(ax, axpos)
tol = eps;
findax = [];
for i = 1:length(ax)
    axipos = get(ax(i),'InnerPosition');
    diffpos = axipos - axpos;
    if (max(max(abs(diffpos))) < tol)
        findax = ax(i);
        break;
    end
end

end

% function lineStyles = linspecer(N)
% This function creates an Nx3 array of N [R B G] colors
% These can be used to plot lots of lines with distinguishable and nice
% looking colors.
% 
% lineStyles = linspecer(N);  makes N colors for you to use: lineStyles(ii,:)
% 
% colormap(linspecer); set your colormap to have easily distinguishable 
%                      colors and a pleasing aesthetic
% 
% lineStyles = linspecer(N,'qualitative'); forces the colors to all be distinguishable (up to 12)
% lineStyles = linspecer(N,'sequential'); forces the colors to vary along a spectrum 
% 
% % Examples demonstrating the colors.
% 
% LINE COLORS
% N=6;
% X = linspace(0,pi*3,1000); 
% Y = bsxfun(@(x,n)sin(x+2*n*pi/N), X.', 1:N); 
% C = linspecer(N);
% axes('NextPlot','replacechildren', 'ColorOrder',C);
% plot(X,Y,'linewidth',5)
% ylim([-1.1 1.1]);
% 
% SIMPLER LINE COLOR EXAMPLE
% N = 6; X = linspace(0,pi*3,1000);
% C = linspecer(N)
% hold off;
% for ii=1:N
%     Y = sin(X+2*ii*pi/N);
%     plot(X,Y,'color',C(ii,:),'linewidth',3);
%     hold on;
% end
% 
% COLORMAP EXAMPLE
% A = rand(15);
% figure; imagesc(A); % default colormap
% figure; imagesc(A); colormap(linspecer); % linspecer colormap
% 
%   See also NDHIST, NHIST, PLOT, COLORMAP, 43700-cubehelix-colormaps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Jonathan Lansey, March 2009-2013 � Lansey at gmail.com               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%% credits and where the function came from
% The colors are largely taken from:
% http://colorbrewer2.org and Cynthia Brewer, Mark Harrower and The Pennsylvania State University
% 
% 
% She studied this from a phsychometric perspective and crafted the colors
% beautifully.
% 
% I made choices from the many there to decide the nicest once for plotting
% lines in Matlab. I also made a small change to one of the colors I
% thought was a bit too bright. In addition some interpolation is going on
% for the sequential line styles.
% 
% 
%%
function lineStyles=linspecer(N,varargin)
if nargin==0 % return a colormap
    lineStyles = linspecer(128);
    return;
end
if ischar(N)
    lineStyles = linspecer(128,N);
    return;
end
if N<=0 % its empty, nothing else to do here
    lineStyles=[];
    return;
end
% interperet varagin
qualFlag = 0;
colorblindFlag = 0;
if ~isempty(varargin)>0 % you set a parameter?
    switch lower(varargin{1})
        case {'qualitative','qua'}
            if N>12 % go home, you just can't get this.
                warning('qualitiative is not possible for greater than 12 items, please reconsider');
            else
                if N>9
                    warning(['Default may be nicer for ' num2str(N) ' for clearer colors use: whitebg(''black''); ']);
                end
            end
            qualFlag = 1;
        case {'sequential','seq'}
            lineStyles = colorm(N);
            return;
        case {'white','whitefade'}
            lineStyles = whiteFade(N);return;
        case 'red'
            lineStyles = whiteFade(N,'red');return;
        case 'blue'
            lineStyles = whiteFade(N,'blue');return;
        case 'green'
            lineStyles = whiteFade(N,'green');return;
        case {'gray','grey'}
            lineStyles = whiteFade(N,'gray');return;
        case {'colorblind'}
            colorblindFlag = 1;
        otherwise
            warning(['parameter ''' varargin{1} ''' not recognized']);
    end
end      
% *.95
% predefine some colormaps
  set3 = colorBrew2mat({[141, 211, 199];[ 255, 237, 111];[ 190, 186, 218];[ 251, 128, 114];[ 128, 177, 211];[ 253, 180, 98];[ 179, 222, 105];[ 188, 128, 189];[ 217, 217, 217];[ 204, 235, 197];[ 252, 205, 229];[ 255, 255, 179]}');
set1JL = brighten(colorBrew2mat({[228, 26, 28];[ 55, 126, 184]; [ 77, 175, 74];[ 255, 127, 0];[ 255, 237, 111]*.85;[ 166, 86, 40];[ 247, 129, 191];[ 153, 153, 153];[ 152, 78, 163]}'));
set1 = brighten(colorBrew2mat({[ 55, 126, 184]*.85;[228, 26, 28];[ 77, 175, 74];[ 255, 127, 0];[ 152, 78, 163]}),.8);
% colorblindSet = {[215,25,28];[253,174,97];[171,217,233];[44,123,182]};
colorblindSet = {[215,25,28];[253,174,97];[171,217,233]*.8;[44,123,182]*.8};
set3 = dim(set3,.93);
if colorblindFlag
    switch N
        %     sorry about this line folks. kind of legacy here because I used to
        %     use individual 1x3 cells instead of nx3 arrays
        case 4
            lineStyles = colorBrew2mat(colorblindSet);
        otherwise
            colorblindFlag = false;
            warning('sorry unsupported colorblind set for this number, using regular types');
    end
end
if ~colorblindFlag
    switch N
        case 1
            lineStyles = { [  55, 126, 184]/255};
        case {2, 3, 4, 5 }
            lineStyles = set1(1:N);
        case {6 , 7, 8, 9}
            lineStyles = set1JL(1:N)';
        case {10, 11, 12}
            if qualFlag % force qualitative graphs
                lineStyles = set3(1:N)';
            else % 10 is a good number to start with the sequential ones.
                lineStyles = cmap2linspecer(colorm(N));
            end
        otherwise % any old case where I need a quick job done.
            lineStyles = cmap2linspecer(colorm(N));
    end
end
lineStyles = cell2mat(lineStyles);
end
% extra functions
function varIn = colorBrew2mat(varIn)
for ii=1:length(varIn) % just divide by 255
    varIn{ii}=varIn{ii}/255;
end        
end
function varIn = brighten(varIn,varargin) % increase the brightness
if isempty(varargin),
    frac = .9; 
else
    frac = varargin{1}; 
end
for ii=1:length(varIn)
    varIn{ii}=varIn{ii}*frac+(1-frac);
end        
end
function varIn = dim(varIn,f)
    for ii=1:length(varIn)
        varIn{ii} = f*varIn{ii};
    end
end
function vOut = cmap2linspecer(vIn) % changes the format from a double array to a cell array with the right format
vOut = cell(size(vIn,1),1);
for ii=1:size(vIn,1)
    vOut{ii} = vIn(ii,:);
end
end
%%
% colorm returns a colormap which is really good for creating informative
% heatmap style figures.
% No particular color stands out and it doesn't do too badly for colorblind people either.
% It works by interpolating the data from the
% 'spectral' setting on http://colorbrewer2.org/ set to 11 colors
% It is modified a little to make the brightest yellow a little less bright.
function cmap = colorm(varargin)
n = 100;
if ~isempty(varargin)
    n = varargin{1};
end
if n==1
    cmap =  [0.2005    0.5593    0.7380];
    return;
end
if n==2
     cmap =  [0.2005    0.5593    0.7380;
              0.9684    0.4799    0.2723];
          return;
end
frac=.95; % Slight modification from colorbrewer here to make the yellows in the center just a bit darker
cmapp = [158, 1, 66; 213, 62, 79; 244, 109, 67; 253, 174, 97; 254, 224, 139; 255*frac, 255*frac, 191*frac; 230, 245, 152; 171, 221, 164; 102, 194, 165; 50, 136, 189; 94, 79, 162];
x = linspace(1,n,size(cmapp,1));
xi = 1:n;
cmap = zeros(n,3);
for ii=1:3
    cmap(:,ii) = pchip(x,cmapp(:,ii),xi);
end
cmap = flipud(cmap/255);
end
function cmap = whiteFade(varargin)
n = 100;
if nargin>0
    n = varargin{1};
end
thisColor = 'blue';
if nargin>1
    thisColor = varargin{2};
end
switch thisColor
    case {'gray','grey'}
        cmapp = [255,255,255;240,240,240;217,217,217;189,189,189;150,150,150;115,115,115;82,82,82;37,37,37;0,0,0];
    case 'green'
        cmapp = [247,252,245;229,245,224;199,233,192;161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27];
    case 'blue'
        cmapp = [247,251,255;222,235,247;198,219,239;158,202,225;107,174,214;66,146,198;33,113,181;8,81,156;8,48,107];
    case 'red'
        cmapp = [255,245,240;254,224,210;252,187,161;252,146,114;251,106,74;239,59,44;203,24,29;165,15,21;103,0,13];
    otherwise
        warning(['sorry your color argument ' thisColor ' was not recognized']);
end
cmap = interpomap(n,cmapp);
end
% Eat a approximate colormap, then interpolate the rest of it up.
function cmap = interpomap(n,cmapp)
    x = linspace(1,n,size(cmapp,1));
    xi = 1:n;
    cmap = zeros(n,3);
    for ii=1:3
        cmap(:,ii) = pchip(x,cmapp(:,ii),xi);
    end
    cmap = (cmap/255); % flipud??
end

