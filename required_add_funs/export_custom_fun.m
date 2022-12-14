function export_custom_fun(funname,options)

arguments
    % mandatory
        funname char
    % optional 
        options.writeDir  = 'current' %['C:\Users\',getenv('USERNAME'),'\Desktop'] %{mustBeMember(options.writeDir,{'current','Desktop'})}
        options.zipFolder {mustBeMember(options.zipFolder,[0,1])} = 0
end
struct2CallerWS(options)

if strcmp(writeDir,'current')
    writeDir = pwd;
end

try
    [fList,pList] = matlab.codetools.requiredFilesAndProducts([funname,'.m']);
catch
    [fList,pList] = matlab.codetools.requiredFilesAndProducts([pwd,'\',funname,'.m']);
end
dir_main     = [writeDir,'\exp_',funname];
dir_add_funs = [dir_main,'\required_add_funs'];
warning('off','all');
mkdir(dir_main)
mkdir(dir_add_funs)
warning('on','all');

%% exclude directory
excl_Casadi = 'C:\Users\Admin\OneDrive - Universität des Saarlandes\WORK_FOLDER\SOLVER\CasADi\casadi-windows-matlabR2016a-v3.5.5';
excl_dirs = {excl_Casadi};
inds_del = [];
for ii=1:length(fList)
    fun_ii = fList{ii};
    for jj=1:length(excl_dirs)
        if contains(fun_ii,excl_dirs{jj})
            inds_del = [inds_del, ii];
        end
    end       
end

inds_keep = setdiff(1:length(fList),inds_del);
fList = fList(inds_keep);


%% copy all required functions
for ii=1:length(fList)
    backslash_pos = strfind(fList{ii},'\');
    cur_funname = fList{ii}(backslash_pos(end)+1:end);
    if strcmp(cur_funname,[funname,'.m'])==0
        copyfile(fList{ii},dir_add_funs)
    end
end
try
    copyfile(which([funname,'.m']),dir_main)
end

%% add info lines to main file and add path to required additional functions
export_fun_lines = ['%This file was auto-exported via export_custom_fun.m with all necessary dependancies',newline,...
                    'mfile_name    = mfilename(''fullpath'');',newline,...
                    '[pathstr,~,~] = fileparts(mfile_name);',newline,...
                    'cd(pathstr);',newline,...
                    'clearvars',newline,...
                    'addpath(''.\required_add_funs'')',newline,...
                    repmat('%',1,100),newline,newline,newline];
str = fileread([dir_main,'\',funname,'.m']);
fid = fopen([dir_main,'\',funname,'.m'],'w');  % open file for writing (overwrite if necessary)
fprintf(fid,'%s%s',export_fun_lines,str);  % Write the char array, interpret newline as new line
fclose(fid);                                   % Close the file (important)
clearvars fid

%% info txt file
info_file = 'required_programs_toolboxes.txt';
fid = fopen([dir_main,'\',info_file],'w');
pList_fns = fieldnames(pList);
for ii=1:length(pList_fns)
    data = pList.(pList_fns{ii});
    if ~isa(data,'char')
        data = num2str(data);
    end
    txt_line = [pList_fns{ii},': ',data,newline];
    fprintf(fid,'%s',txt_line);
end
fclose(fid);

if zipFolder==1
    zip([writeDir,'\','exp_',funname,'.zip'],dir_main)
    cmd_rmdir(dir_main);
end

end


