%GRAPH SIGNAL ANALYSIS 

data=zX_RS;

%% compute fMRI HF/LF portions
clear X_hat X_c X_d X_m N_c N_d X_all 
for s=1:size(data,3)
    X_hat(:,:,s)=M'*data(:,:,s);
    X_c(:,:,s)=Vlow*X_hat(:,:,s);
    X_d(:,:,s)=Vhigh*X_hat(:,:,s);
    X_all(:,:,s)=M*X_hat(:,:,s); % reconstruct back full signal to check norm
    
    %% norms  of the weights
    for r=1:n_ROI
        N_c(r,s)=norm(X_c(r,:,s));
        N_d(r,s)=norm(X_d(r,:,s));
    end
end

%% mean across subjects
mean_c=mean(N_c,2); %average coupling
mean_d=mean(N_d,2); %average decoupling
