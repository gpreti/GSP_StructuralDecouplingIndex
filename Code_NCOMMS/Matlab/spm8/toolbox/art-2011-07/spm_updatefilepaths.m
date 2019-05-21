
function var=spm_updatefilepaths(varargin)
% SPM_UPDATEFILEPATHS(SPM) Checks file paths in SPM structure and asks the user
% to relocate any incorrect file references (e.g. when the data has been moved 
% to a different directory or to a different computer)
%
% SPM_UPDATEFILEPATHS('SPM.mat');
% Updates file SPM.mat in current directory (the file SPM.mat will be
% overwritten; a file named old_SPM.mat will contain the original
% information)
%

% alfnie@gmail.com 01/09

persistent changed begin fullname1 fullname2 fullnamematch changescount;

if nargin==1,changed=0;changescount=0;varname='';var=varargin{1}; cumdisp;
    if ischar(varargin{1}), load(varargin{1},'SPM'); var=SPM; var0=var; varname='SPM'; [filepath,filename,fileext]=fileparts(varargin{1}); cwd=pwd; if ~isempty(filepath), cd(filepath); end; end
else var=varargin{1}; varname=varargin{2}; 
end

cumdisp(['Checking ',varname]);
if ~isstruct(var), return; end
if numel(var)==1 && isfield(var,'fname'),
    if  isempty(dir(var.fname)),
        begin=1;
        filename=strtrim(deblank(var.fname));
        switch(filesep),case '\',idx=find(filename=='/');case '/',idx=find(filename=='\');end; filename(idx)=filesep;
        n=0; 
        while n<size(filename,1),
            ok=dir(deblank(filename(n+1,:)));
            if ~isempty(ok), n=n+1;
            else
                if begin && changed, begin=0;
                else
                    fullname1=deblank(filename(n+1,:));
                    [pathname1,name1,ext1]=fileparts(fullname1);
                    [name2,pathname2]=uigetfile(['*',ext1],['File not found in field ',varname],['*',name1,ext1]);
                    %[name2,pathname2]=uigetfile(['*',ext1],['File not found: ',name1,ext1],[name1,ext1]);
                    if all(name2==0), filename=[]; break; end
                    fullname2=fullfile(pathname2,name2);
                    changed=1;begin=0;
                end
                fullnamematch=strvcat(fliplr(fullname1),fliplr(fullname2));
                m=sum(cumsum(fullnamematch(1,:)~=fullnamematch(2,:))==0);
                m1=max(0,length(fullname1)-m); m2=max(0,length(fullname2)-m);
                filename=strvcat(filename(1:n,:),[repmat(fullname2(1:m2),[size(filename,1)-n,1]),filename(n+1:end,m1+1:end)]);
            end
        end
        if ~isempty(filename),
            cumdisp(['Updating ',varname]);
            V=spm_vol(filename);
            var=V;
            changescount=changescount+1;
        else var=[]; 
        end
    else
        cumdisp(['Checking ',varname,' ok']);
    end
elseif numel(var)>1,
    clear vartemp;
    for nvar=1:numel(var),
        temp=spm_updatefilepaths(var(nvar),[varname,'(',num2str(nvar),')']);
        if ~isempty(temp), vartemp(nvar)=temp; elseif ~isempty(var(nvar)), var=[]; break; end
    end
    if ~isempty(var), var=reshape(vartemp,size(var)); end
else
    varfieldnames=fieldnames(var);
    for nvar=1:length(varfieldnames),
        temp=spm_updatefilepaths(var.(varfieldnames{nvar}),[varname,'.',varfieldnames{nvar}]);
        if ~isempty(temp), var.(varfieldnames{nvar})=temp; elseif ~isempty(var.(varfieldnames{nvar})), var=[]; break; end
    end
end

if nargin==1,
    if isfield(var,'xY')&&isfield(var.xY,'P')&&isfield(var.xY,'VY'),var.xY.P=char({var.xY.VY(:).fname});end
    cumdisp('');cumdisp;
    disp(['Done. ',num2str(changescount),' updates']);
    if ischar(varargin{1})&&changescount>0, [filepath,filename,fileext]=fileparts(varargin{1}); SPM=var0; save(fullfile(filepath,['old_',filename,fileext]),'SPM'); SPM=var; save(varargin{1},'SPM'); cd(cwd); end
end
end

function cumdisp(txt)
% CUMDISP persistent disp
% cumdisp; initializes persistent display
% cumdisp(text); displays persistent text
%

persistent oldtxt;

if nargin<1,
    oldtxt=''; 
    fprintf(1,'\n'); 
else
    fprintf(1,[repmat('\b',[1,length(oldtxt)]),txt]);
    oldtxt=sprintf(txt);
end
end

