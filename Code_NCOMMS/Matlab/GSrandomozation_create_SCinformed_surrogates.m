%% Create SC-informed graph signal surrogates by randomization of real harmonics Fourier coefficients
clear XrandS 
nSurr=19;

%% SPATIAL RANDOMIZATION
X_hat_rand=zeros(size(data,1),size(data,2),size(data,3),nSurr);
for s=1:size(data,3)
    X=data(:,:,s);
    for n=1:nSurr
        %randomize sign of Fourier coefficients
        PHIdiag=round(rand(size(M,1),1));
        PHIdiag(PHIdiag==0)=-1;
        PHI=diag(PHIdiag);
        XrandS{s}(:,:,n)=M*PHI*M'*X; % X_hat=M'X, normally reconstructed signal would be Xrecon=M*X_hat=MM'X, instead of M, M*PHI is V with randomized signs
        
    end
end




