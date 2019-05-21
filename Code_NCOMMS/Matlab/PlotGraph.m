%% Brain plots
%% Atlas + Selection
N=360; %number of regions in Glasser atlas
CM=zeros(N,N);
CodeBookpath=which('Glasser360_2mm_codebook.mat');
CodeBook=load(CodeBookpath);
CodeBook=CodeBook.codeBook;

%% adjust Cvalues for saturation (to eliminate outliers peaks)
if saturate
thr=1;
CC2new=CC2;
CC2new(find(CC2>thr))=0;
CC2new(find(CC2>thr))=max(CC2new);
CC2new(find(CC2<-thr))=0;
CC2new(find(CC2<-thr))=min(CC2new);
CC2=CC2new;
end

%% plot with normal color scheme 
CC=abs(CC2);
T_conn=0;
Factor_SphereSize=max(CC);
Factor_Col=1;
Exp_Sphere=2;
View=2;

Colormap_edges='jet';
Colormap_nodes='jet';

Gamma=0.5;
LinearWeight=1;
CA=[-1 1]; 


%%
PlotBrainGraph(CM,CC,CC2,CodeBook,T_conn,Factor_SphereSize,...
    Factor_Col,Exp_Sphere,View,Colormap_nodes,Colormap_edges,Gamma,...
    LinearWeight,CA)
