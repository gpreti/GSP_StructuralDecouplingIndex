%% find SDI for SC-ignorant surrogates
%%filter SC-ignorant surrogates into coupled/decoupled portions (based on real SC harmonics)

fdata=XrandSran;

for s=1:size(fdata,2)
    for i=1:size(fdata{1},3)
        X_hat_surr_ran{s}(:,:,i)=M'*fdata{s}(:,:,i);
        X_c_surr_ran{s}(:,:,i)=Vlow*X_hat_surr_ran{s}(:,:,i);
        X_d_surr_ran{s}(:,:,i)=Vhigh*X_hat_surr_ran{s}(:,:,i);
        %% norms (or variance) of the weights
        for r=1:size(fdata{1},1) %for each region
            N_c_surr_ran(r,i,s)=norm(X_c_surr_ran{s}(r,:,i));
            N_d_surr_ran(r,i,s)=norm(X_d_surr_ran{s}(r,:,i));
            
        end
    end
end

SDI_surr_ran=N_d_surr_ran./N_c_surr_ran; % SDI for every subject and surrogates
SDI_surr_ran_avgsurr=squeeze(mean(SDI_surr_ran,2)); % mean across surr
SDI_surr_ran_avgsurrsubjs=mean(SDI_surr_ran_avgsurr,2); %%AVERAGE SDI OF SC-ignorant SURROGATE SIGNALS, FILTERED BASED ON EMPIRICAL GRAPH (Fig. 2A)

%% plot Fig. 2A
saturate=1;
CC2=log2(SDI_surr_ran_avgsurrsubjs); 
PlotGraph;title('Fig. 2A')