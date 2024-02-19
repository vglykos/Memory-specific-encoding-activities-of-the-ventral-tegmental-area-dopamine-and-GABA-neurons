% This script performs the following calculations
% (1) Produces the average firing rate for Left and Right trials within a
% single session
% (2) Produces the firing rate difference between left and right trials
% (3) Produces N permutations of the trial trajectories and then executes step (1) and step (2)
% (4) From the N permuted L vs R rate differences it calculates the
% statistical threshold (2 standard deviations).
% (5) Compares the statistical threshold to the original firing rate difference.
% (6) When the original rate difference is larger or smaller to the
% statistical threshold at a certain maze location, then these location
% points are delignated with red asteriscs (bottom plot)
% 
% 
% Key notes
% (a) To produce the firing rate (spikes/sec) we have to divide the number of spikes occuring in a certain 
%     position, by the time the animal spent in this position (occupation time).
% (b) Error trials were excluded from this analysis


%%
close all; clear; clc
%%
load('CNC845','merge');  % load the mat file
unit_list=[845]; % the unit id or a list of units

%% user settings
p=struct;
%firing rate smoothing parameters (gaussian)
p.conv.x=linspace(-1,1,20); 
p.conv.mi=0; 
p.conv.sigma=0.5; 
p.conv.pad=0;

p.Fcam=merge.user.Fcam; % camera frame rate

% parameters for statistics
p.Ppoint=0.05;  
p.Pglobal=0.05;  
p.Nshuffles=500; %Num of permutations

% saving parameters
p.save_results=0;
p.save_fig=0; 
p.close_fig=0; 
p.fig_position=1; %position of plots / Left=1 , Right=2
p.start=0; %choose sensor to begin or set to zero for all trial

% call the main function
[tests] = psth_comparison_fun(unit_list , merge , p ,'LR')


%%
function [tests] = psth_comparison_fun(list , merge , p , save_name)
%(1) extract the timestamps of spikes of the individual unit
%(2) assign trajectories to the spikes
%(3) calculate the occupation time for every location in every trial from
%the number of camera frames
%(4) call "psth_difference_fun" which produces the original firing rate
%difference along with the permutation
%(5) plot these results and save them (or not)
%(6) the actual code is a lot simpler than it seems, but I had to make it
%longer in order to fix bugs, mainly caused by animal behavior. Also, I had
%to write a lot of code for arranging and saving the results.

%import vars
protocol = merge.protocol;
unit_log = merge.units.unit_log;
psth = merge.units.psth_data.psth; psth(:,2)=[];
spikes = merge.units.spikes_data.spikes;
occup = merge.position.occup_real;
Z = merge.position.Zreal;
Z_sensor = merge.position.Zreal_sensor;
Z_axis = merge.position.Z_axis;
group1 = merge.trials.group1;
group2 = merge.trials.group2;


% create vars
comparison=struct;
comparison.labels={'unit_id', 'psth-1', 'psth-2', 'orig diff', ...
    'breaks-point', 'breaks-global', 'sig'};
comparison.results=cell(size(unit_log,1) , 7);

%%


%% generate unit list to be analyzed
if isempty(list)
    list=cell2mat( unit_log(:,1) );
end

%%
%start unit analysis
for i_list=1:1:numel(list)
    
    i_unit=find(cell2mat( unit_log(:,1) ) == list(i_list));
    
    
    %% is the unit silent? (MANY QUESTIONS WHU I HAVE THIS CODE HERE)
    if isempty(cell2mat(psth(i_unit,2)))
        disp(sprintf('unit%d is silent' , unit_log{i_unit,1}))
        if p.start
            index = Z_axis >= X_sensor(p.start);
        else
            index = Z_axis >= 0;
        end
        NN=numel(find(index)); clear index
        comparison.results(i_unit,:)={list(i_list) , zeros(1,NN) , ...
        zeros(1,NN) , zeros(1,NN) , zeros(1,NN) , zeros(1,NN) , 0};
    else
        
        
        %% extract spikes
        spikes_ = spikes{i_unit,2};
        
        %group spikes by arm choice
        m=ismember(spikes_(:,1),group1);  
        spikes1=[ spikes_(m,1), spikes_(m,7) ]; 
        sz1=size(spikes1,1);
        
        m=ismember(spikes_(:,1),group2);  
        spikes2=[ spikes_(m,1), spikes_(m,7) ]; 
        sz2=size(spikes2,1);
        
        %arrange spikes so left & right trials are plotted in seperatre groups 
        u=unique(spikes1(:,1));
        for i=1:1:numel(u)
            index = spikes1(:,1)==u(i);
            spikes1(index,1) = i;
        end
        M=max(spikes1(:,1));
            
        u=unique(spikes2(:,1));   
        for i=1:1:numel(u)
            index = spikes2(:,1)==u(i);
            spikes2(index,1) = i;   
        end
        spikes2(:,1) = M + spikes2(:,1);

        
    %% psth % occupation time 
    histograms=psth{i_unit,2};
    occup_time=cellfun(@transpose, occup, 'un', 0); 
    occup_time=cell2mat(occup_time);
    
    %% Isolate the trial window I want to analyze
    if p.start
        index = Z_axis >= Z_sensor(p.start);
        zaxis = Z_axis(index);
        z_sensor = Z_sensor(p.start : end);
        ZTickLabel = {'1' '2' '3' '4' '5' '6/8' '7/9'}; 
        ZTickLabel = ZTickLabel(p.start : end);
        histograms = histograms(:,index);
        occup_time = occup_time(:,index);
        %spikes
        index=spikes1(:,2) >= z_sensor(1); spikes1 = spikes1(index,:);
        index=spikes2(:,2) >= z_sensor(1); spikes2 = spikes2(index,:);  
        
    else
        index = Z_axis >= 0;
        zaxis = Z_axis(index);
        z_sensor = Z_sensor;
        ZTickLabel = {'1' '2' '3' '4' '5' '6/8' '7/9'};
        histograms = histograms(:,index);
        occup_time = occup_time(:,index);
        %spikes
        index=spikes1(:,2) >= 0; spikes1 = spikes1(index,:);
        index=spikes2(:,2) >= 0; spikes2 = spikes2(index,:);
    end
    
    
    %% comparisons
    p.zaxis=zaxis;
    out = psth_difference_fun(histograms ,occup_time ,group1, group2, p);
    %update results?   
    comparison.results(i_unit,:)={list(i_list) , out.array1 , ...
        out.array2 , out.D0 , out.breaks_point5 , out.breaks_global5 , out.sig};
      
            %% plot
            if p.fig_position==1; fig_position=[10 150 700 800]; end;   if p.fig_position==2; fig_position=[770 150 700 800]; end
            f1=figure('NumberTitle','off','position', fig_position,'color','k','visible','on');
            fig_name=sprintf('%d-sh%d-clu%d-%s-%s-%s' , unit_log{i_unit,1} , unit_log{i_unit,5} , unit_log{i_unit,6} , unit_log{i_unit,2} , unit_log{i_unit,3} , protocol);  
            save_fig_name = sprintf('%d-%s' , unit_log{i_unit,1} , protocol)
            set(f1,'Name', fig_name);

            set(gca,'position', [0.13 0.96 0.775 0.001],'color','k')
            title(fig_name,'Edgecolor','y','Background','y' , 'Interpreter', 'none')
            
            ax321=axes('position',[0.1 0.73 0.80 0.230]);   hold on; 
            plot(zaxis,out.array1,'color','c');  
            plot(zaxis,out.array2,'color','r');

            for i=1:1:numel(z_sensor)
                line([z_sensor(i) z_sensor(i)] , [0  max([out.array1 , out.array2])] , 'color','y','linestyle','-.')
            end
            
            if p.start
                set(ax321,'XColor','y' , 'YColor','y','color','k','box','on' , 'TickDir','out' , 'xlim',[z_sensor(1) z_sensor(end)])
            else
                set(ax321,'XColor','y' , 'YColor','y','color','k','box','on' , 'TickDir','out' , 'xlim',[0 z_sensor(end)])
            end
            
            set(gca , 'Xtick',z_sensor , 'XtickLabel',ZTickLabel,'FontSize',8, 'Box','on');
            descr={'gr1=blue' ; 'gr2=red'}; text(1.02 ,max(get(ax321,'ylim')), descr,'color','y')
%             title('group1=blue / group2=red','color','y'); ylabel('spikes/sec')
    
    
            ax323=axes('position',[0.1 0.47 0.80 0.230]); hold on
            plot(zaxis , out.D0 , 'w');
            plot(zaxis , out.Dpoint(out.point5,:), 'color',[0 0.5 0], 'linestyle','-.');
            plot(zaxis , out.Dpoint(end - out.point5,:), 'color',[0 0.5 0], 'linestyle','-.');
    
            if ~isempty(out.global5)
                plot(zaxis , out.Dpoint(out.global5,:), 'g-');
                plot(zaxis , out.Dpoint(end-out.global5,:), 'g-');
            end
            plot(zaxis(find(out.breaks_global5)) , out.D0(find(out.breaks_global5)) , 'r*' , 'markersize' , 3.5)
            
            for i=1:1:numel(z_sensor)
                line([z_sensor(i) z_sensor(i)], [-max([out.array1, out.array2])  max([out.array1, out.array2])] , 'color','y','linestyle','-.')
            end
            
            if p.start
                set(ax323,'XColor','y' , 'YColor','y','color','k','TickDir','out' , 'xlim',[z_sensor(1) z_sensor(end)])
            else
                set(ax323,'XColor','y' , 'YColor','y','color','k','TickDir','out' , 'xlim',[0 z_sensor(end)])
            end
                
            set(ax323 , 'Xtick',[z_sensor] , 'XtickLabel',ZTickLabel,'FontSize',8,'Box','on')
            ylabel('spikes/sec')
            
            
            ax325=axes('position',[0.1 0.055 0.80 0.380]); hold on
            scatter(spikes1(:,2), spikes1(:,1), 1,'marker', '.' , 'markerfacecolor','c' , 'markeredgecolor','c')  
            scatter(spikes2(:,2), spikes2(:,1), 1, 'marker', '.' , 'markerfacecolor','r' , 'markeredgecolor','r')
            M=max([ spikes1(:,1) ; spikes2(:,1)]); if isempty(M); M=1; end
            
            if p.start
                set(ax325,'XColor','y' , 'YColor','y','color','k','TickDir','out' , 'xlim',[z_sensor(1) z_sensor(end)]  ,'ylim', [0 M+1] ,'FontSize',8,'Box','on')
            else
                set(ax325,'XColor','y' , 'YColor','y','color','k','TickDir','out' , 'xlim',[0 z_sensor(end)]  ,'ylim', [0 M+1] ,'FontSize',8,'Box','on')
            end
            xlabel('one dimension'); ylabel('trials'); hold on
    
    
    end %end of IF is psth empty
    
    %% save figure if it exists
    if exist('f1')
        % save figure
        if p.save_fig; 
            savefig(f1,save_fig_name); 
            set(f1,'InvertHardCopy','Off' , 'PaperPositionMode','Auto');
            print(f1,save_fig_name,'-dpng' , '-r0')
        end
        
        % close fig
        if p.close_fig
            pause(0.2); close
        end
        
        clear f1
    end %exist('f1')
    
    
    %% update tests
    comparison.group1=group1;  
    comparison.group2=group2;
    eval( sprintf('tests.psth.%s=comparison;', save_name) );

  
end %for i_list=1:1:numel(list)


%% save results
if p.save_results    
    save(protocol, '-append', 'tests')
end

disp(' ');disp(' ');
disp('*** compare_groups complete ****')

end %end of function



%%
function out = psth_difference_fun(array,occup,row1,row2,p)
Npoints=numel(p.zaxis); %zaxis samples

%% Calculate the original firing rate difference between left and right trials
%IMPORTANT: Take into consideration the occupation time

%Group-1 psth
array1 = array(row1,:);  
if isnan(array1);   array1=ones(Npoints)*NaN;   end %fix bugs

occup1 = occup(row1,:);  
if isnan(occup1);   occup1=ones(Npoints)*NaN;   end %fix bugs

occup1=occup1/p.Fcam;  occup1(occup1==0)=0.0001; %I dont want occup elements = 0. But again array1 will have 0 elements.
array1=array1./occup1;

array1=mean(array1,1);

[out]=convolution_fun(array1 , p.conv); %smoothing
array1=out.data;
 

%Group-2 psth
array2 = array(row2,:);  
if isnan(array2);   array2=ones(Npoints)*NaN;   end %fix bugs

occup2 = occup(row2,:); 
if isnan(occup2);   occup2=ones(Npoints)*NaN;   end %fix bugs

occup2=occup2/p.Fcam;  occup2(occup2==0)=0.0001; %I dont want occup elements = 0. But again array1 will have 0 elements.
array2=array2./occup2;

array2=mean(array2,1); %smoothing

[out]=convolution_fun(array2 , p.conv);
array2=out.data;

% Original difference
D0 = array1 - array2;


%% shuffling (Shuffle the Left or Right tags assigned on spikes)
in.factors = [ones(numel(row1),1) ; 2*ones(numel(row2),1)];
in.data = [array(row1,:) ; array(row2,:)];
in.denom = [occup(row1,:) ; occup(row2,:)];

    D=shuffling_fun(in,p);

%global & pointwise bands
in.D0=D0;
in.D=D;
in.Ppoint=p.Ppoint;
in.Pglobal=p.Pglobal;
in.Nshuffles=p.Nshuffles;
in.binranges=p.zaxis;

%NEXT
%(1)produce the statistical threshold and 
%(2)check if they are violated and 
%(3)if so, in which maze locations
o=stats_fun(in); 
        
sig=o.sig;  
breaks_global5=o.breaks_global5;
breaks_point5=o.breaks_point5;
Dpoint=o.Dpoint;
point5=o.point5;
global5=o.global5;



% output variables
out=struct;
out.array1=array1; out.array2=array2;
out.occup1=occup1; out.occup2=occup2;
out.D0=D0;
out.point5=point5;
out.Dpoint=Dpoint;
out.global5=global5;
out.breaks_global5=breaks_global5;
out.breaks_point5=breaks_point5;
out.sig=sig;

end %function out=compare_f(in);


function out = shuffling_fun(in,p)
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


function [out]=convolution_fun(in,param)


    %% parameters
    mi = param.mi;
    sigma = param.sigma;
    x = param.x;
    
    
%     x=linspace(-1,1,20); 
%     mi=0; 
%     sigma=0.5; 
%     param.pad=1;
    
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