close all
clear 
clc
%%
file ="glm2344"
load(file)
%%
%original average Left firing rate
L0 = meta.L0;
%original average Right firing rate
R0 = meta.R0;


%GLM regressor design
model = ['rate ~ ' ...
        'trials + accuracy + old_accuracy + performance +'  ...
        'speed + speed^2 + speed^3 + speed^4 + speed^5 + speed^6 +'...
        'pos*choice + pos^2*choice + pos^3*choice + pos^4*choice + pos^5*choice + pos^6*choice'];

%fit model with original trajectories
%"T" variable is like a pandas dataframe
mdl_0 = fitglm(T,model,'Distribution','normal');
R2_0 = mdl_0.Rsquared.Adjusted;

%model predictions
L_= T.choice==1;
yhatL = predict(mdl_0,T(L_,:)); yhatL = reshape(yhatL,100,[]); yhatL0 = nanmean(yhatL,2);

R_= T.choice==2;
yhatR = predict(mdl_0,T(R_,:)); yhatR = reshape(yhatR,100,[]); yhatR0 = nanmean(yhatR,2);

% plot
f=figure;plot(L0,'g');hold on;plot(R0,'m');plot(yhatL0,'g.');hold on;plot(yhatR0,'m.'); title(sprintf('%s R2: %4.3f',file, R2_0))
legend({'Lorig','Rorig','Lhat','Rhat'},'Location','northeast') 
