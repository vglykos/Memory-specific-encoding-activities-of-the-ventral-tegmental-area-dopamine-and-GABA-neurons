%% This script calculates the response of units to light pulses

%The stimulation-recording protocol consisted of sets of pulses with different amplitude. 
%Usually, we delivered 150 light pulses / amplitude. 
%The number of amplitudes was variable (3 - 5 different sets of amplitude, 1,2,3,4,5mW)
%For example with 4 amplitude sets we delivered 4 sets X 150 pulses = 600 pulses 

% The script produces a figure. 
% In the left side we plot the rasteplots of single APs 
%(red: APs fired during the light pulse / white: APs fired during the dark period)

%In the right side we plot the histograms of the firing activities during light stimulation.
%Also, we plot a horizontal green line. This green line illustrates the
%global threshold of the spontaneous activity. If the hidtogram exceeds the
% green light, then the light-induced response is significant and the unit
% is consdered as light-responsive (DAT unit or Vgat unit).

%%
clear all; close all; clc

%% load unit data
load("unit1539")

%% Define variables or extract variable from the mat file
% the number of shuffles from which we will extract the global statistical threshold
Nshuffles = 500; 

% the global threshold limit (p=0.001). 
P_global = 0.001;   

% histogram bin size
bin = user.bin; 

%in the onset & offset of the light pulse there is 1ms artefact 
%Light pulses had 13ms duration. Therefore, we will eliminate from our
%analysis the first and last ms.
artef = user.impulse_artef;  artef=artef/bin; 

% sampling frequency of the recording equipment
Fs = user.Fs;

%Npulse contais the total number of light pulses (for all amplitudes)
Npulse = stimulus.Npulse; 

%set Amplitudes (e.g. 0.5mW, 1.0mW, 1.5mW, 2.0mW, 2.5mW)
Nset = stimulus.Nset;

% the number of pulses per amplitude
Npulse_per_set = Npulse/Nset;

%light-pulse duration in sec
Tpulse = stimulus.Tpulse; Tpulse = round(Tpulse*1000)/1000; 

%the timestamps of the onset and offset of every stimulus 
stim_on = stimulus.stim_on;

%we defined 13ms windows preceding the light pulses. The rasteplots of
%spontaneous activity are extracted from those windows.
stim_off = stimulus.stim_off;


%the x_axis of the rasteplots and the histograms
xaxis = [0:bin:Tpulse];
%But as I described above, we need to remove the 1ms artefact.
clean_axis = [artef+1 : numel(xaxis)-artef]; 

%% DO I NEED THAT?
% light.results = unit_log;
%     light.labels = {'unitID', 'animal', 'genotype', 'date', 'record', ...
%     'shank', 'clu', 'spikes_per_stim (percent)', 'cv', 'xaxis', ...
%     'light-hist', 'global_thres', 'sig_points', 'sig_hist', ...
%     'sig_record', 'total_sig_points',  'latency', 'point_thres'};
%%
% %% figure
% f1=figure('NumberTitle','off' , 'position',[30 150 800 800] ,'color','k');
% fig_name=sprintf('unit-%d  sh%d-clu%d  %s-%s-%s' , unit_log{i_,1} , unit_log{i_,5} , unit_log{i_,6} , unit_log{i_,2} , unit_log{i_,4} , unit_log{i_,5});  
% set(f1,'Name', fig_name);
% set(gca,'position', [0.13 0.96 0.775 0.001],'color','k'); axis off
% title(fig_name,'Edgecolor','y','Background','y')

%%
% stim_on = stimulus.stim_on;
% % stim_on = stimulus.t_on;
% stim_off = stimulus.stim_off;
% % stim_off = stimulus.t_off;

% Npulse = stimulus.Npulse;
% Nset = stimulus.Nset;


%% arrange the indeces of the y-axis in the rasterplot. 
%It is not part of the analysis. It is only for illustration purposes. 
%If you do not understand it, just ignore it.
wnd=[ [1 : Npulse/Nset : 2*Npulse]' , [Npulse/Nset : Npulse/Nset : 2*Npulse]' ];        
wnd_light=wnd([1:2:2*Nset],:);
wnd_spont=wnd([2:2:2*Nset],:);
 
Yaxis_light=[]; Yaxis_spont=[];    
for i=1:1:Nset
    Yaxis_light = [ Yaxis_light ; [wnd_light(i,1):1:wnd_light(i,2)]' ];
    Yaxis_spont = [ Yaxis_spont ; [wnd_spont(i,1):1:wnd_spont(i,2)]' ];
end
     


%% RASTEPLOT: Extract the timestamps of APs fired during the light pulses

binranges = arrayfun(@(x,y) x:1/2E4:y , stim_on(:,1),stim_on(:,2) , 'UniformOutput',false);
[counts,ind ] = cellfun(@(x) histc(ts, x), binranges, 'UniformOutput',false);
ind = cellfun(@(x) x(x~=0), ind, 'UniformOutput',false);
AP_events = cellfun(@(x,y,z) transpose(y(x)-z) , ind,binranges,num2cell(stim_on(:,1)),  'UniformOutput',false); 
index = cellfun(@(x,y) ones(size(x,1),1)*y, AP_events, num2cell([1:1:size(AP_events,1)]'), 'UniformOutput',false); 
AP_events = [cell2mat(AP_events) , cell2mat(index)];

if ~isempty(AP_events)    
    light = [AP_events , Yaxis_light(AP_events(:,2))];
else    
    light = nan(1,3);
end


%% RASTERPLOT: Extract the timestamps of APs fired during spontaneous activity (i.e. 13ms before every light pulse)  
binranges = arrayfun(@(x,y) x:1/2E4:y , stim_off(:,1),stim_off(:,2) , 'UniformOutput',false);
[counts,ind ] = cellfun(@(x) histc(ts, x), binranges, 'UniformOutput',false);
ind = cellfun(@(x) x(x~=0), ind, 'UniformOutput',false);
AP_events = cellfun(@(x,y,z) transpose(y(x)-z) , ind,binranges,num2cell(stim_off(:,1)),  'UniformOutput',false); 
index = cellfun(@(x,y) ones(size(x,1),1)*y, AP_events, num2cell([1:1:size(AP_events,1)]'), 'UniformOutput',false); 
AP_events = [cell2mat(AP_events) , cell2mat(index)];
    
if ~isempty(AP_events)
    spont=[AP_events , Yaxis_spont(AP_events(:,2))];    
else    
    spont=nan(1,3);
end
       

% keep the data for the figure
raster.light = light;
raster.spont = spont;
raster.Y = wnd_spont(:,2);






%% (1) produce light response histograms, (2) the global statistical thresh from spont activity
% I produce histograms for every amplitude (1mW, 2mW, 3mW, 4mW, 5mW)

%Therefore, I group the indeces of pulses by the mplitude set. I call them
%windows
wnd=[ [1 : Npulse/Nset : Npulse]' , [Npulse/Nset : Npulse/Nset : Npulse]' ];



for i_set=1:1:Nset 
    %break time event vectors into smaller ones
    stim_on_set = stimulus.stim_on( [wnd(i_set,1) : wnd(i_set,2) ] , :);
    % Produce 'stim_off_set'
    stim_off_set = stimulus.stim_off( [wnd(i_set,1) : wnd(i_set,2) ] , :);
%     %dark period
%     dark_period = stimulus.dark_period(1:end-1,:);
%     dark_period_set = dark_period([wnd(i_set,1) : wnd(i_set,2) ] , :);
            
    %To make code run faster I reckon I should break the timestamp
    %'ts' of the whole session in sets. To do so i can use the
    %'stim_on_set' & 'stim_off_set'. 
    ind_ = all( [ts >= stim_off_set(1,1) , ts <= stim_on_set(end,2)] ,2 );
    ts_set = ts(ind_);
            
    %Light-response histogram     
    %create binranges
    binranges = arrayfun(@(x,y) x:bin:y , stim_on_set(:,1),stim_on_set(:,2) , 'Un',0);
    %histogram
    counts = cellfun(@(x) histc(ts_set, x), binranges, 'Un',0);
    %fix bug
    [sz] = cellfun(@size, counts,'Un',0);   sz=cell2mat(sz); 
    gt=find(sz(:,1)>sz(:,2));    
    counts(gt)=cellfun(@(x) transpose(x) , counts(gt), 'Un',0); 
    %bug fixxed
    counts=cell2mat(counts); 
    %final product
    light_hist{i_set,1} = [mean(counts,1)];
%   spikes_per_stim{i_set,1} = 100 * (numel(find(sum(counts,2)==1)) / (Npulse/Nset)) ;
    spikes_per_stim{i_set,1}=sum(counts,2);
    light_cv{i_set,1} = std(sum(counts,2)) / mean(sum(counts,2));  
    clear binranges sz counts gt sz
               
                
    %Spontaneous activity
    % to calculate the statistical threshold from the
    % spontaneous activity, we will jitter the windows of
    % spontaneous activity (stim_off) N times (defined by
    % Nshuffles).
                

    Npulses_per_set = size(stim_off_set,1);
    hist_x_axis_samples = numel(stim_off_set(1,1):bin:stim_off_set(1,2));
            
    %For every stim_off window produce Nshuffles (jittered)
    repmat_ = repmat(stim_off_set(:,2) , Nshuffles , 1);
    jittering_factor = rand(size(repmat_,1) , 1);
    jitter_ = Tpulse + jittering_factor * Tpulse;
                
    %Finally I produce the jittered stim_off windows
    stim_off_shuffle = [repmat_ - jitter_  , repmat_ - jitter_ + Tpulse]; 
    clear jittering_factor
    %create binranges
    binranges = arrayfun(@(x,y) x:bin:y , ...
    stim_off_shuffle(:,1),stim_off_shuffle(:,2) , 'Un',0);
    %histograms
    counts = cellfun(@(x) histc(ts_set, x), binranges, 'Un',0);
    %fix bug
    [sz] = cellfun(@size, counts,'Un',0);   sz=cell2mat(sz); 
    gt = find(sz(:,1)>sz(:,2));        
    counts(gt) = cellfun(@(x) transpose(x) , counts(gt), 'Un',0);
    %bug fixxed
    counts = cell2mat(counts); 
    %mat2cell
    counts = mat2cell(counts,[ones(Nshuffles,1)*Npulses_per_set] , hist_x_axis_samples);
    %produce the mean
    counts = cellfun(@(x) mean(x,1) , counts, 'UniformOutput',false);
    %final product
    spont_hist = cell2mat(counts);   
    clear binranges sz counts gt sz stim_off_shuffle Npulses_per_set ...
    hist_x_axis_samples
                
    %calculate the statistical global threshold, and fid out if the
    %light-induced histogram exceeds it.
    global_ = spont_hist(:,clean_axis);
    global_ = reshape(global_,[],1);
                
    global_ = sort(global_,'descend'); 
    row_max = round( P_global/2 * numel(global_) ); 
    global_thres{i_set,1} = global_(row_max);
    clear global_
    sig_points{i_set,1} = double(light_hist{i_set,1} > global_thres{i_set,1});
                
                
    if sum(sig_points{i_set,1})
        sig_hist{i_set,1}=1;
    else
        sig_hist{i_set,1}=0;
    end
                

 end %for i_set=1:1:Nset
 
%% figure
f1=figure('NumberTitle','off' , 'position',[30 150 800 800] ,'color','k');

%rasterplot    
ax1=axes('position',[0.05 0.05 0.45 0.9]); hold on
set(ax1,'XColor','y' , 'YColor','y','color','k','TickDir','out' , 'xlim',[0 Tpulse] ,'FontSize',8,'Box','on');
xlabel('time (sec)'); ylabel('stimuli'); hold on
scatter(raster.spont(:,1), raster.spont(:,3), 0.5 ,'markerfacecolor','w' , 'markeredgecolor','w')  
scatter(raster.light(:,1), raster.light(:,3), 0.5 , 'markerfacecolor','r' , 'markeredgecolor','r')
for i=1:1:Nset;  line( [0 Tpulse] , [raster.Y(i) raster.Y(i)] , 'color','y'); end
ylim([1 2*Npulse+1])

 for i_set=1:1:Nset
            axes('position',stimulus.axes(i_set,:)); hold on
            set(gca,'XColor','y' , 'YColor','y','color','k','TickDir','out' , ...
                'xlim',[0 Tpulse] ,'FontSize',8,'Box','on'); ylabel('pmf');
            
            if i_set==1
                xlabel('one dimension'); 
            else
                set(gca,'XTickLabel' , []); 
            end
        
            bar(xaxis,light_hist{i_set,1},'barwidth',1, 'facecolor','r', 'edgecolor','y')
            clean_axis=[artef+1 : numel(xaxis)-artef];
            line( [xaxis(clean_axis(1)) xaxis(clean_axis(end))] , [global_thres{i_set,1} global_thres{i_set,1}] , 'color','g')
 end %for Iset=1:1:Nset;