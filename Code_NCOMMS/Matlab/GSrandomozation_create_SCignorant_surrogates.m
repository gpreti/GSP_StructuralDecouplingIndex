%% Create SC-ignorant graph signal surrogates by randomization of configuration model derived-SC harmonics Fourier coefficients

%% 1) Graph structure "randomization": configuration model used to generate degree-preserving graph 
for n=1:1 %one can also repeat some times, then take the mean
    Wran1{n}=rand_graph(Wnew);
    Wran=Wran1;
end

%% 2) Create surrogate signals based on configuration model graph : XrandSran
%create average config.model SC and decompose it
mean_Wran=zeros(360,360);
for i=1:n
    mean_Wran=mean_Wran+Wran{i};
end
mean_Wran=mean_Wran./n;

%%Laplacian
meanLran=eye(n_ROI)-mean_Wran;
%%Laplacian Decomposition
[Uran,LambdaLran] = eig(meanLran);
[LambdaLran, IndLran]=sort(diag(LambdaLran));
Uran=Uran(:,IndLran);
meanMran=fliplr(Uran);

%%Create Surrogates
nSurr=19;
for s=1:size(data,3)
    X=zX_RS(:,:,s);
    for n=1:nSurr
        %randomize sign of Fourier coefficients
        PHIdiag=round(rand(size(meanMran,1),1));
        PHIdiag(PHIdiag==0)=-1;
        PHI=diag(PHIdiag);
        XrandSran{s}(:,:,n)=meanMran*PHI*meanMran'*X; 
        
    end
end




