
%% Generate maps for Neurosynth analysis (which will be in Python, code from Margulies 2016)
%The maps are saved in the folder data, inside the code folder

%select reference nifti map (MNI)
ref_nifti_path=which('HCP-MMP1_onMNI152_2mm_Glasser360.nii');
refhdr=spm_vol(ref_nifti_path);
refnii=spm_read_vols(refhdr);

%select map to split into percentiles
map_to_consider=mean_SDI;
map_to_consider(isnan(map_to_consider)==1)=0;

NSmaps=Neurosynth_createmaps(map_to_consider,refnii);

mkdir(strcat(mypath,'results/my_masks'));
%save maps
for i=1:9
    refhdr.fname=strcat(mypath,'results/my_masks/56subjs_SDI_0',num2str(i),'.nii');
    spm_write_vol(refhdr,NSmaps{i});
end
for i=10:20
    refhdr.fname=strcat(mypath,'results/my_masks/56subjs_SDI_',num2str(i),'.nii');
    spm_write_vol(refhdr,NSmaps{i});
end

%% for visualization of gradient chopped in 5th percentiles
sequence=[-1:0.1053:1];
sequence=[sequence,1];
sequence_percentiles=prctile(mean_SDI,5:5:100);
for i=1:360 %for each node, which percentile it belongs to?
    [dist,index]=min(abs(mean_SDI(i)-sequence_percentiles));
    if mean_SDI(i)>sequence_percentiles(index)
        mean_SDI_chopped(i)=sequence(index+1);
    else
        mean_SDI_chopped(i)=sequence(index);
    end
end
