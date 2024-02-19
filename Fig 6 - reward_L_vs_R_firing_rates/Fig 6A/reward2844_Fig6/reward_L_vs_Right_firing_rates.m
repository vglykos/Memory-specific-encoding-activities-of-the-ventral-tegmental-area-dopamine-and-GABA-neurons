clear all;
close all; 
clc
%%
load("reward2844" , 'merge')
unit_list=[2844]; % the unit id or a list of units
%% load variables
unit_log = merge.units.unit_log; 
wnd = [0 1]; %time window :wnd(1) time before 1st lick / wnd(2) time after 1st lick
Nbins = 1000; %num of histogram bins
xaxis = linspace(0 , diff(wnd) , diff(wnd) * Nbins); % create the time axis

%% more parameters
p=struct;

%gaussian filter smoothing
p.conv.x=linspace(-1 , 1, numel(xaxis)/5); 
p.conv.mi=0; 
p.conv.sigma=0.5; 
p.conv.pad=0;

%camera sampling frequency
p.Fcam=1;

%number of permutations
p.Nshuffles=500;

%statistical thresholds
p.Ppoint=0.05;  
p.Pglobal=0.05;  


%%
%(1) for every trial extract the timestamps of the first lick which activates
%the water pump.
%(2) interpolate the 1 sec window from the ts of the 1st lick until a
%second later
%(3) the number of samples of interpolations is N
first_lick_L = merge.reward.first_lick_L;
first_lick_R = merge.reward.first_lick_R;
L_ts = first_lick_L(:,3); edges_L = [L_ts-wnd(1) L_ts + wnd(2)];
R_ts = first_lick_R(:,3); edges_R = [R_ts-wnd(1) R_ts + wnd(2)];

edges_L=[];
for ii=1:1:size(L_ts,1)
    edges_L = [edges_L ; linspace(L_ts(ii) - wnd(1) , L_ts(ii) + wnd(2) , numel(xaxis)) ];
end

edges_R=[];
for ii=1:1:size(R_ts,1)
    edges_R = [edges_R ; linspace(R_ts(ii) - wnd(1) , R_ts(ii) + wnd(2) , numel(xaxis)) ];
end

%% start the unit analysis
% I made a for-loop because sometimes I analyze a batch of units
for i_list = 1:1:numel(unit_list)
    
    i_unit=find(cell2mat( unit_log(:,1) ) == unit_list(i_list));
    
    unit_ID = unit_log{i_unit,1};
    % extract the timestamps of the unit
    unit_ts = merge.units.spikes_data.spikes{i_unit,2}(:,4);
    %create the histograms
    L_hist  = double(ismember(round(edges_L * Nbins) , round(unit_ts * Nbins)));
    R_hist  = double(ismember(round(edges_R * Nbins) , round(unit_ts * Nbins)));
    hist_i = [L_hist ; R_hist];
    
    %average the original histograms
    L_resp_0 = nanmean(L_hist / xaxis(2)); 
    [out]=convolution_fun(L_resp_0 , p.conv); 
    L_resp_0 = out.data;
    
    R_resp_0 = nanmean(R_hist / xaxis(2));
    [out]=convolution_fun(R_resp_0 , p.conv); 
    R_resp_0 = out.data;
    
    %calculate the  original difference
    D0 = L_resp_0 - R_resp_0;
    
    %% Perform the N permutations
    in.data = hist_i;
    in.factors = [ones(size(L_hist,1),1) ; ones(size(R_hist,1),1)*2];
    in.denom = ones(size(hist_i));
    p.Fcam=1/xaxis(2);
    
    D = shuffling_fun(in,p);
    %NOTE. This may look funny. I was lazy to write a new shuffling code.
    %Thus, I used the one I wrote for the position-arrange spikes. 
    %So, I divided number of spikes by the occupation time. But all
    %histogram bins had occupation time set to 1. So, all of them were
    %divided by the same denomitator with the original histograms
    
    %global & pointwise bands
    in.D0 = D0;
    in.D = D;
    in.Ppoint = p.Ppoint;
    in.Pglobal = p.Pglobal;
    in.Nshuffles = p.Nshuffles;
    in.binranges = xaxis;
    
    %Calculate the statistical thresholds and make the comparions.
    o=stats_fun(in);
    
    sig = o.sig;  
    breaks_global5 = o.breaks_global5;
    breaks_point5 = o.breaks_point5;
    Dpoint = o.Dpoint;
    point5 = o.point5;
    global5 = o.global5;
    
    %dominant vector
    dominant = zeros(size(D0));
    dominant(D0 > 0)=1;
    dominant(D0 < 0)=2;
    
    %% scatterplot_data
    L_scatter = L_hist .* xaxis;
    [r,c,v] = find(L_scatter); [r,I] = sort(r); v = v(I);
    L_scatter = [r,v]; clear r c v
    
    R_scatter = R_hist .* xaxis;
    [r,c,v] = find(R_scatter); [r,I] = sort(r); v = v(I);
    R_scatter = [r,v]; clear r c v
    if ~isempty(L_scatter)
        R_scatter(:,1) = R_scatter(:,1) + max(L_scatter(:,1));
    else
        R_scatter(:,1) = R_scatter(:,1);
    end
    

    %% figure
        
    f1=figure('NumberTitle','off');
    set(f1,'position',[10 150 700 800],'Name',sprintf('unit%d',unit_list(i_list)));
   
    ax1=axes('position',[0.1 0.73 0.80 0.220]);   hold on; 
    plot(xaxis , L_resp_0,'b'); hold on
    plot(xaxis , R_resp_0,'r');
%     title(save_fig_name,'color',[0.23 0.44 0.6]) 
    ylabel('spikes/sec')
    dim=get(ax1,'ylim');
        
    [start,finish,clr] = significance_bar_fun(breaks_global5);
    
    for ii = 1:size(start,1)
        if clr(ii)==1; c='b'; end
        if clr(ii)==2; c='r'; end
        plot([xaxis(start(ii)) xaxis(finish(ii))] , [dim(2)+0.1*dim(2) dim(2)+0.1*dim(2)],'linewidth',5,'color',c); hold on
    end
    
    line([wnd(1) wnd(1)] , [get(ax1,'ylim')],'color','w','linestyle',':');
   
    ax2 = axes('position',[0.1 0.47 0.80 0.220]); hold on
    plot(xaxis,D0,'m'); hold on
    if ~isempty(global5)         
        plot(xaxis , Dpoint(global5,:), 'g-');         
        plot(xaxis , Dpoint(end - global5,:), 'g-');
    end
    plot(xaxis(find(breaks_global5)) , D0(find(breaks_global5)) , 'r*' , 'markersize' , 3.5)
    line([wnd(1) wnd(1)] , [get(ax2,'ylim')],'color','w','linestyle',':')
    ylabel('spikes/sec')
    
    ax3=axes('position',[0.1 0.055 0.80 0.380]); hold on
    scatter(L_scatter(:,2), L_scatter(:,1), 1, ...
        'marker', '.' , 'markerfacecolor','c' , 'markeredgecolor','c'); hold on           
    scatter(R_scatter(:,2), R_scatter(:,1), 1, ...
        'marker', '.' , 'markerfacecolor','r' , 'markeredgecolor','r');
   
    if ~isempty(R_scatter)
        set(ax3, 'ylim', [0 max(R_scatter(:,1))]);
    end
    line([wnd(1) wnd(1)] , [get(ax3,'ylim')],'color','w','linestyle',':')
    ylabel('licks')
    xlabel('time (sec)')
    
    c = [0.23 0.44 0.6];
    set(f1, 'color', 'k')
    set(ax1,'XColor',c,'YColor',c,'color','k','box','on','TickDir','out','FontSize',8);
    set(ax2,'XColor',c,'YColor',c,'color','k','box','on','TickDir','out','FontSize',8);
    set(ax3,'XColor',c,'YColor',c,'color','k','box','on','TickDir','out','FontSize',8,'xlim',[xaxis(1) xaxis(end)]);         
    
end %for i=1:1:Nunits


%%
function out=shuffling_fun(in,p)
% inputs
% 1:Num of shuffles (Nshuffles)
% 2:a matrix (MxN) containing the original data (data)
% 3:a matrix (Mx1) with the factors. (1,1,1,2,2,1) (factors)
h = waitbar(0,'Shuffling in progress....Please wait...', 'position',[-1430  10  270.0000   56.2500]);

Nshuffles=p.Nshuffles;
array = in.data;
rows = in.factors;
occup = in.denom;

out=ones(Nshuffles , size(array,2));

row1=find(rows==1);
row2=find(rows==2);

array = [ array(row1,:) ; array(row2,:) ];

occup = [ occup(row1,:) ; occup(row2,:) ];


L1=length(row1);
L2=length(row2);


for i=1:1:Nshuffles
        waitbar(i/Nshuffles);   %disp(i);
        %%
        I=randperm(L1 + L2);
        array = array(I,:);
        occup = occup(I,:);
        %%
        array1 = array(1:L1,:);    
        if isnan(array1);   array1=ones(N)*NaN;   end %fix bugs
        
        occup1 = occup(1:L1,:); 
        if isnan(occup1);   occup1=ones(N)*NaN;   end %fix bugs
        occup1=occup1/p.Fcam;  occup1(occup1==0)=0.0001; %I dont want occup elements = 0. But again array1 will have 0 elements.
        
        array1=array1./occup1;
        array1=mean(array1,1);
        
        [out1]=convolution_fun(array1 , p.conv);
        array1=out1.data;

        %%
        array2 = array(L1+1:L1+L2,:);    
        if isnan(array2);   array2=ones(N)*NaN;   end %fix bugs
        
        occup2 = occup(L1+1:L1+L2,:); 
        if isnan(occup2);   occup2=ones(N)*NaN;   end %fix bugs
        occup2=occup2/p.Fcam;  occup2(occup2==0)=0.0001; %I dont want occup elements = 0. But again array1 will have 0 elements.
        
        array2=array2./occup2;
        array2=mean(array2,1);
        
        
        [out2]=convolution_fun(array2 , p.conv);
        array2=out2.data;
        
        %%
        D_shuf = array1 - array2;
      
        if size(D_shuf,1) > size(D_shuf,2); D_shuf = D_shuf'; end
        
        out(i,:)=D_shuf;
        

%         figure;plot(D_shuf,'k'); hold on; plot(array1,'b'); plot(array2,'r')
    end %for i=1:1:Nshuffles

close(h); pause(0.1)
end %End of function


%%
function out=stats_fun(in)
%(1)produce the statistical threshold and 
%(2)check if they are violated and 
%(3)if so, in which maze locations

D0=in.D0;
D=in.D;
Ppoint=in.Ppoint;
Pglobal=in.Pglobal; 
Nshuffles=in.Nshuffles;
binranges=in.binranges;

if size(D,1)~=Nshuffles && size(D,2)~=Nshuffles
    error('dimension of Nshuffles and D dont agree') 
end

if size(D,1)~=Nshuffles
    D=D';
end

Dglobal=D; %unsorted
 
Dpoint=sort(D,1);   
Dpoint=flipud(Dpoint); %sorted
    
point5=round(Ppoint*(Nshuffles/2));
global_band=[]; pointwise_band=[];
    
    
for row=point5:-1:1
    Igt=bsxfun(@gt,Dglobal,Dpoint(row,:));
    Ilt=bsxfun(@lt,Dglobal,Dpoint(end-row+1,:));
    ind=sum(Igt,2)+sum(Ilt,2);
    ind=find(ind);
      
    pointwise_band=[pointwise_band ; row/(Nshuffles/2)];
    global_band=[global_band ; length(ind)/Nshuffles];
end

%Find the pointwise band at which less than Pglobal of the shuffled data break it
row_break=find(global_band<Pglobal);     

%if there is no pointwise band which less than 5% of global databreak it, the report error
% if isempty(row_break); 
%     error('Cannot find <5% global data breaking the pointwise band'); 
% end
if isempty(row_break)
    breaks_global5=zeros(size(binranges));
    global5=[];
else
    row_break=row_break(1); %keep the first one
    global5=point5-row_break+1;
    %find if the original data break the 5% global band
    up=D0>Dpoint(global5,:);  up=double(up);
    down=D0<Dpoint(end-global5+1,:);   down=double(down); down(down==1)=-1;
    breaks_global5=[up + down];
end

if sum(breaks_global5)~=0  
    sig=1;
    up=D0>Dpoint(point5,:);  up=double(up);
    down=D0<Dpoint(end-point5+1,:); down=double(down); down(down==1)=-1;
    breaks_point5=[up+down];    
     
else
    sig=0;
    breaks_point5=zeros(size(binranges));   
end

%if there is no pointwise band which less than 5% of global databreak it,
%the report error
if isempty(row_break)
    sig=2;  
end


out.sig=sig;
out.breaks_global5=breaks_global5;
out.breaks_point5=breaks_point5;
out.Dpoint=Dpoint;
out.point5=point5;
out.global5=global5;

end %End of function

%%
function [out]=convolution_fun(in,param)


    %% parameters
    mi = param.mi;
    sigma = param.sigma;
    x = param.x;
    
    %% convolution function
    conv_wnd=2.5*( 1/(sigma*sqrt(2*pi)) ) * exp(-(x-mi).^2 / (2*sigma^2) );
    conv_wnd=conv_wnd/sum(conv_wnd);

    y = conv_e_f(in, conv_wnd);
    if size(y,2) ~= size(in,2)
        y=y';
    end
    out=struct;
    out.data=y;
%     out.data=y';
    out.conv_param.mi=mi;
    out.conv_param.sigma=sigma;
    out.conv_param.x=x;

end %function


%%
function y = conv_e_f(x, h)

% if mod(h/2) == 0; error('second input should be odd'); end 
if size(x,2)>=100; x=x';end
% if size(x,1)<=size(x,2); x=x' ; end
if size(h,1)<=size(h,2); h=h' ; end
    
y = nanconv(x,h/sum(h),'same');

Nv = numel(h);
Nha = floor(numel(h)/2);
Nhb = round(numel(h)/2);

for ii=1:Nha
    if ii==10
        d=0;
    end
    % Left edge
    v1 = h((Nhb-ii+1):end); % Including the center point
    u1 = x(1:numel(v1)); if size(u1,1) < size(u1,2); u1=u1'; end
    v1(isnan(u1)) = NaN;
    v1 = v1 / nansum(v1); if size(v1,1) < size(v1,2); v1=v1'; end
    y(ii) = nansum(u1.*v1);
    
    % Right edge
    v2 = h(1:(Nhb-1+ii)); % Including the center point 
    u2 = x((end-numel(v2)+1):end); if size(u2,1) < size(u2,2); u2=u2'; end
    v2(isnan(u2)) = NaN;
    v2 = v2 / nansum(v2); if size(v2,1) < size(v2,2); v2=v2'; end
    y(end-ii+1) = nansum(u2.*v2);
end
end

%%
function c = nanconv(a, k, varargin)
% NANCONV Convolution in 1D or 2D ignoring NaNs.
%   C = NANCONV(A, K) convolves A and K, correcting for any NaN values
%   in the input vector A. The result is the same size as A (as though you
%   called 'conv' or 'conv2' with the 'same' shape).
%
%   C = NANCONV(A, K, 'param1', 'param2', ...) specifies one or more of the following:
%     'edge'     - Apply edge correction to the output.
%     'noedge'   - Do not apply edge correction to the output (default).
%     'nanout'   - The result C should have NaNs in the same places as A.
%     'nonanout' - The result C should have ignored NaNs removed (default).
%                  Even with this option, C will have NaN values where the
%                  number of consecutive NaNs is too large to ignore.
%     '2d'       - Treat the input vectors as 2D matrices (default).
%     '1d'       - Treat the input vectors as 1D vectors.
%                  This option only matters if 'a' or 'k' is a row vector,
%                  and the other is a column vector. Otherwise, this
%                  option has no effect.
%
%   NANCONV works by running 'conv2' either two or three times. The first
%   time is run on the original input signals A and K, except all the
%   NaN values in A are replaced with zeros. The 'same' input argument is
%   used so the output is the same size as A. The second convolution is
%   done between a matrix the same size as A, except with zeros wherever
%   there is a NaN value in A, and ones everywhere else. The output from
%   the first convolution is normalized by the output from the second 
%   convolution. This corrects for missing (NaN) values in A, but it has
%   the side effect of correcting for edge effects due to the assumption of
%   zero padding during convolution. When the optional 'noedge' parameter
%   is included, the convolution is run a third time, this time on a matrix
%   of all ones the same size as A. The output from this third convolution
%   is used to restore the edge effects. The 'noedge' parameter is enabled
%   by default so that the output from 'nanconv' is identical to the output
%   from 'conv2' when the input argument A has no NaN values.
%
% See also conv, conv2
%
% AUTHOR: Benjamin Kraus (bkraus@bu.edu, ben@benkraus.com)
% Copyright (c) 2013, Benjamin Kraus
% $Id: nanconv.m 4861 2013-05-27 03:16:22Z bkraus $

% Process input arguments
for arg = 1:nargin-2
    switch lower(varargin{arg})
        case 'edge'; edge = true; % Apply edge correction
        case 'noedge'; edge = false; % Do not apply edge correction
        case {'same','full','valid'}; shape = varargin{arg}; % Specify shape
        case 'nanout'; nanout = true; % Include original NaNs in the output.
        case 'nonanout'; nanout = false; % Do not include NaNs in the output.
        case {'2d','is2d'}; is1D = false; % Treat the input as 2D
        case {'1d','is1d'}; is1D = true; % Treat the input as 1D
    end
end

% Apply default options when necessary.
if(exist('edge','var')~=1); edge = false; end
if(exist('nanout','var')~=1); nanout = false; end
if(exist('is1D','var')~=1); is1D = false; end
if(exist('shape','var')~=1); shape = 'same';
elseif(~strcmp(shape,'same'))
    error([mfilename ':NotImplemented'],'Shape ''%s'' not implemented',shape);
end

% Get the size of 'a' for use later.
sza = size(a);

% If 1D, then convert them both to columns.
% This modification only matters if 'a' or 'k' is a row vector, and the
% other is a column vector. Otherwise, this argument has no effect.
if(is1D);
    if(~isvector(a) || ~isvector(k))
        error('MATLAB:conv:AorBNotVector','A and B must be vectors.');
    end
    a = a(:); k = k(:);
end

% Flat function for comparison.
o = ones(size(a));

% Flat function with NaNs for comparison.
on = ones(size(a));

% Find all the NaNs in the input.
n = isnan(a);

% Replace NaNs with zero, both in 'a' and 'on'.
a(n) = 0;
on(n) = 0;

% Check that the filter does not have NaNs.
if(any(isnan(k)));
    error([mfilename ':NaNinFilter'],'Filter (k) contains NaN values.');
end

% Calculate what a 'flat' function looks like after convolution.
if(any(n(:)) || edge)
    flat = conv2(on,k,shape);
else flat = o;
end

% The line above will automatically include a correction for edge effects,
% so remove that correction if the user does not want it.
if(any(n(:)) && ~edge); flat = flat./conv2(o,k,shape); end

% Do the actual convolution
c = conv2(a,k,shape)./flat;

% If requested, replace output values with NaNs corresponding to input.
if(nanout); c(n) = NaN; end

% If 1D, convert back to the original shape.
if(is1D && sza(1) == 1); c = c.'; end

end

%%
function [start,finish,clr] = significance_bar_fun(in)
%%
sz = size(in);
if sz(1) < sz(2)
    in=in';
end

%%
z=[[in ; 0] , [0 ; in]];
              
start_ = all([ any([z(:,1)==-1 , z(:,1)==1],2) , z(:,2)==0 ] , 2);
finish_ = all([ z(:,1)==0 , any([z(:,2)==-1 , z(:,2)==1],2) ] , 2);

start_(end)=[]; start = find(start_);
finish_(1)=[]; finish = find(finish_);

%%
clr = zeros(size(start));
for i=1:1:size(start,1)
    if in(start(i))==1
        clr(i) = 1;
    end
    if in(start(i))==-1
        clr(i) = 2;
    end
end

%bugs
for i=1:1:size(start,1)
    if start(i) == finish(i) && start(i)==1
        start(i) = 1;
        finish(i) = 2;
    end
    if start(i) == finish(i) && start(i)==100    
        start(i) = 99;    
        finish(i) = 100;
    end
    if start(i) == finish(i) && start(i)~=1 && start(i)~=100     
        start(i) = start(i) - 1;
    end
      
end %for i=1:1:size(start,1)


        
end %EOF