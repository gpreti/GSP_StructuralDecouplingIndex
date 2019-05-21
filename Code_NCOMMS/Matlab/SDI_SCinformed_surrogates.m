%% find coupled/decoupled components of surrogate signals and test significance of structural decoupling index (SDI)

%% 1) FILTER coupled/decoupled signal portions in surrogates

fdata=XrandS;

clear X_hat_surr X_c_surr X_d_surr X_m_surr N_c_surr N_d_surr
for s=1:size(fdata,2)
    for i=1:size(fdata{1},3)
        X_hat_surr{s}(:,:,i)=M'*fdata{s}(:,:,i);
        X_c_surr{s}(:,:,i)=Vlow*X_hat_surr{s}(:,:,i);
        X_d_surr{s}(:,:,i)=Vhigh*X_hat_surr{s}(:,:,i);
        %% norms  of the weights
        for r=1:size(fdata{1},1)
            N_c_surr(r,i,s)=norm(X_c_surr{s}(r,:,i));
            N_d_surr(r,i,s)=norm(X_d_surr{s}(r,:,i));
            
        end
    end
end

%% Find significant SDI
%consider and test mean across subjects
SDI_surr=N_d_surr./N_c_surr; % SDI for every subject and surrogates
SDI_surr_avgsurr=squeeze(mean(SDI_surr,2)); % mean across surr
SDI_surr_avgsurrsubjs=mean(SDI_surr_avgsurr,2); %% MEAN SDI SURR SHOWN IN FIGURE fig. 2B

mean_SDI=mean_d./mean_c; %empirical AVERAGE SDI
SDI=N_d./N_c; %emipirical individual SDI

%find threshold for max
%for every subject, max across surrogates
for s=1:size(SDI_surr,3)
max_SDI_surr(:,s)=max(SDI_surr(:,:,s)')';
end

%find threshold for min
%for every subject, in across surrogates
for s=1:size(SDI_surr,3)
min_SDI_surr(:,s)=min(SDI_surr(:,:,s)')';
end

%% select significant SDI for each subject, across surrogates 
%individual thr, first screening
for s=1:size(fdata,2) %for each subject, I threshold the ratio based on individual ratio's surrogate distribution 
    significant_values_max(s)=size(find(SDI(:,s)>max_SDI_surr(:,s)),1);
    significant_values_min(s)=size(find(SDI(:,s)<min_SDI_surr(:,s)),1);
    SDI_thr_max(:,s)=SDI(:,s)>max_SDI_surr(:,s);
    SDI_thr_min(:,s)=SDI(:,s)<min_SDI_surr(:,s);
    detect_max=sum(SDI_thr_max'); %amounts of detection per region
    detect_min=sum(SDI_thr_min');
end

%%for every region, test across subjects 0.05, correcting for the number of
%%tests (regions), 0.05/360
x=0:1:100;
y=binocdf(x,100,0.05,'upper');
THRsubjects=x(min(find(y<0.05/360))); 
THRsubjects=floor(size(fdata,2)/100*THRsubjects)+1;

SDI_sig_higher=detect_max>THRsubjects;
SDI_sig_lower=detect_min>THRsubjects;

SDI_sig_higher_positions=find(SDI_sig_higher==1);
SDI_sig_lower_positions=find(SDI_sig_lower==1);

SDI_sig_tot_positions=[find(SDI_sig_higher==1),find(SDI_sig_lower==1)];
SDI_sig_tot_positions=sort(unique(SDI_sig_tot_positions));

%%threshold empirical mean ratios
mean_SDI_thr=ones(360,1);
mean_SDI_thr(SDI_sig_tot_positions)=mean_SDI(SDI_sig_tot_positions);


%% Fig. 2B
saturate=1;
CC2=log2(SDI_surr_avgsurrsubjs); 
PlotGraph;title('Fig. 2B')
%% Fig. 2C
saturate=1;
CC2=log2(mean_SDI_thr); 
PlotGraph;title('Fig. 2C')

%% Fig. S3
saturate=0;
CC2=log2(mean_SDI_thr./SDI_surr_avgsurrsubjs); 
PlotGraph;title('Supplementary Fig. S3')
