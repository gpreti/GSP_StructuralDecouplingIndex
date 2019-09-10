%% Create SC-ignorant graph signal surrogates by randomization of configuration model derived-SC harmonics Fourier coefficients

%% 1) Graph structure "randomization": configuration model used to generate degree-preserving graph 

Wcm=sum(Wnew,2)*sum(Wnew,2)'./sum(sum(Wnew,2));

%% 2) Create surrogate signals based on configuration model graph : XrandSran

Lcm=diag(sum(Wnew,2))-Wcm;
%%Laplacian Decomposition
[Ucm,LambdaLcm] = eig(Lcm);
[LambdaLcm, IndLcm]=sort(diag(LambdaLcm));
Ucm=Ucm(:,IndLcm);
Mcm=fliplr(Ucm);

%% Create Surrogates
nSurr=19;
for s=1:nsubjs_RS
    X=zX_RS(:,:,s);
    for n=1:nSurr
        %randomize sign of Fourier coefficients
        PHIdiag=round(rand(size(Mcm,1),1));
        PHIdiag(PHIdiag==0)=-1;
        PHI=diag(PHIdiag);
        XrandSran{s}(:,:,n)=Mcm*PHI*Mcm'*X; 
        
    end
end




