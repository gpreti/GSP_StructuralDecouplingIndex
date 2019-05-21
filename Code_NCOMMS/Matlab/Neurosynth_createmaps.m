%Create maps as inputs for Neurosynth 

%% map back atlas based values to nifti file
function NSmaps=Neurosynth_createmaps(vec_to_consider,refnii)


%% division in percentiles bins like Margulies 2016
for i=1:20
   thr1=percentile(vec_to_consider,(i-1)*5);
   thr2=percentile(vec_to_consider,i*5);

   map1=vec_to_consider>thr1;
   map2=vec_to_consider<=thr2;
   map=map1+map2;
   NSvecs(:,i)=map==2;
   
   f=find(NSvecs(:,i)==1);
   NSmaps{i}=zeros(size(refnii));
   for j=1:size(f,1)
   NSmaps{i}(round(refnii)==f(j))=1;
   end
end

end