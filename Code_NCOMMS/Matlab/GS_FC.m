%% FC of empirical signals
for s=1:size(zX_RS,3)
    CSpacereal(:,:,s)=corr(zX_RS(:,:,s)');
end
mean_CSpacereal=mean(CSpacereal,3);
figure;imagesc(mean_CSpacereal,[0 0.7]);title('Fig. 1D - Empirical FC')

%% FC of SC-informed surrogates
count=1;
for s=1:size(XrandS,2)
    for i=1:nSurr
        CSpace(:,:,count)=corr(XrandS{s}(:,:,i)'); 
        count=count+1;
        
    end
end
mean_CSpacesurr=mean(CSpace,3);
figure;imagesc(mean_CSpacesurr,[0 0.7]);title('Fig. 1D - SC-informed FC')
%% FC of SC-ignorant surrogates
for s=1:size(XrandSran,2)
    for i=1:nSurr
        CSpace3{i}(:,:,s)=corr(XrandSran{s}(:,:,i)');
    end
end
for i=1:nSurr
    mean_CSpacesurr3(:,:,i)=mean(CSpace3{i},3); %mean across subjects
end
figure;imagesc(mean(mean_CSpacesurr3,3),[0 0.7]);title('Fig. 1D - SC-ignorant FC')

%% compute nodestrengths 
ns_real=sum(abs(mean_CSpacereal))'; %empirical FC
ns_surr=sum(abs(mean_CSpacesurr))';%SC-informed FC
ns_surr_surr=sum(abs(mean(mean_CSpacesurr3,3)))';%SC-ignorant FC
ns_W=sum(abs(W))'; %SC

%correlation FC-SC
[h1,p1]=corr(ns_W,ns_surr, 'type','Spearman');%0.94
[h2,p2]=corr(ns_W,ns_real, 'type','Spearman');%0.46
