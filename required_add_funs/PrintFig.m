function PrintFig(filename,options)
    %% function argument validation
    arguments
        % mandatory
            filename (1,:) char
        % optional 
            options.Resolution (1,1) double = 400
            options.FigHandle               = gcf
            options.Directory  (1,:) char   = ''
            options.FileFormat     {mustBeMember(options.FileFormat,{'pdf','svg','png','eps','epsc','eps2','epsc2','fig','meta','ps','psc','ps2','psc2','jpeg','tiff','tiffn','meta','bmpmono','bmp','bmp16m','bmp256','hdf','pbm','pbmraw','pcxmono','pcx24b','pcx256','pcx16','pgm','pgmraw','ppm','ppmraw'})} ...
                                                                                                   = 'pdf'
            options.MaximizeFigure {mustBeMember(options.MaximizeFigure,{'true','false'})}         = 'false'
            options.MaximizeOption {mustBeMember(options.MaximizeOption,{'-fillpage','-bestfit'})} = '-fillpage'
    end
    struct2CallerWS(options)
    
    %% set print options and print
    FilenameFormat  = sprintf('%s.%s',filename,FileFormat);
    if strcmp(FileFormat,'fig')==0
        printNameValuePairs = {};
        if strcmp(MaximizeFigure,'false')
            set(FigHandle,'Units','Inches');
            pos = get(FigHandle,'Position');
            set(FigHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        elseif strcmp(MaximizeFigure,'true')
            printNameValuePairs = [printNameValuePairs,MaximizeOption];
        end
        ResSpecifier    = ['-r',int2str(Resolution)];
        FormatSpecifier = ['-d',FileFormat];
        print(ResSpecifier,FigHandle,FormatSpecifier,FilenameFormat,printNameValuePairs{:});
    else
        saveas(FigHandle,filename,'fig');
    end

    %% move the file to Directory
    if ~isempty(Directory)
        if isfolder(Directory)==1
%             [status,msg] = movefile(FilenameFormat,Directory,'f');  
%             if status==0
%                 warning(['Moving the file was not succesful: ',msg])
%             end
        elseif isfolder(Directory)==0
%             warning('Moving the file was not succesful: Folder does not exist.')
            mkdir(Directory)
        end
        [status,msg] = movefile(FilenameFormat,Directory,'f');  
        if status==0
            warning(['Moving the file was not succesful: ',msg])
        end
    else
        % keep file in current directory
    end
end