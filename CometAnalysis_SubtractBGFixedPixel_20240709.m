% % %
% % % Created by GA on 2018-08-10
% % % Last edited by GA on 2018-11-15
%line 319 and 320 - edit axis

% % % Example usage: CometAnalysis_SubtractBGFixedPixel('R1_20uM_P/individualLS',15,25)

function CometAnalysis_SubtractBGFixedPixel_20240709(foldername,solpix,latpix)
try
    
    close all
    clearvars -except foldername solpix latpix
    
    solpix=-1*solpix;
    
    deletefiles=sprintf('%s/BlatSubtracted-*.dat',foldername);
    delete(deletefiles);
    
    fnrfirstLS=sprintf('%s/*_LS0001.dat',foldername);
    FilesFirstLS=dir(fnrfirstLS);
    NFirstLS=length(FilesFirstLS);
    
    for FirstLSno=1:NFirstLS
        close all
        FirstLSName=sprintf('%s/%s',foldername,FilesFirstLS(FirstLSno).name);
        
        kymoprefix = extractBefore(FirstLSName,'_LS0001.dat');
        fprintf('working on %s\n',kymoprefix)
        
        fnr=sprintf('%s_LS*.dat',kymoprefix);
        Files=dir(fnr);
        Nls=length(Files);
        indexmin=0;
        indexmax=0;
        %%% find min and max index
        for lsno=1:Nls
            clearvars -except foldername FilesFirstLS NFirstLS FirstLSno Files Nls lsno indexmin indexmax fnaveraged kymoprefix solpix latpix
            
            FileNames=sprintf('%s/%s',foldername,Files(lsno).name);
            
            fidr=fopen(FileNames,'r');
            tracknotflipped=fscanf(fidr,'%f %f',[2 Inf]);
            fclose(fidr);
            tracknotflipped=tracknotflipped';
            
            track=zeros(size(tracknotflipped,1),2);
            for i=1:size(tracknotflipped,1)
                track(i,1)=(-1.0)*tracknotflipped(i,1);
                track(i,2)=tracknotflipped(i,2);
            end
            
            if (indexmin>min(track(:,1)))
                indexmin=min(track(:,1)); % get the mimimum x-value
            end
            if (indexmax<max(track(:,1)))
                indexmax=max(track(:,1)); % get the maximum x-value
            end
        end
        
        
        indexminint=indexmin-mod(indexmin,1); % is it's not integer, make it an integer
        indexmaxint=indexmax-mod(indexmax,1)+1; % is it's not integer, make it an integer
        
        %fprintf('found:\nmin index=%d\nmaxindex=%d\n',indexminint,indexmaxint)
        
        if (indexminint<solpix)
            indexminint=solpix;
        end
        if (indexmaxint>latpix)
            indexmaxint=latpix;
        end
        %fprintf('using:\nmin index=%d\nmaxindex=%d\n',indexminint,indexmaxint)
        
        f1=figure;
        LS=zeros((indexmaxint-indexminint+1),3);
        %%% this loop just sums
        for k=indexminint:indexmaxint
            %clc
            %fprintf('Kymo %d of %d: %f percent completed\n',FirstLSno,NFirstLS,100*(k-indexminint)/(indexmaxint-1-indexminint))
            for lsno=1:Nls % this is linescan loop
                clear tracknotflipped track
                FileNames=sprintf('%s/%s',foldername,Files(lsno).name);
                fidr=fopen(FileNames,'r');
                %             [filepath,name,ext] = fileparts(FileNames);
                %             newStr = extractAfter(name,'K');
                %             kymonoextracted = extractBefore(newStr,'_');
                %             kymono=str2num(kymonoextracted);
                %             fprintf('%s\t%s\t%d\n',name,kymonoextracted,kymono)
                tracknotflipped=fscanf(fidr,'%f %f',[2 Inf]);
                fclose(fidr);
                tracknotflipped=tracknotflipped';
                
                
                track=zeros(size(tracknotflipped,1),2);
                for i=1:size(tracknotflipped,1)
                    track(i,1)=(-1.0)*tracknotflipped(i,1);
                    track(i,2)=tracknotflipped(i,2);
                end
                
                Nx=size(track,1);
                
                if (k==indexminint)
                    set(0, 'CurrentFigure', f1)
                    plot(track(:,1),track(:,2))
                    xlim([solpix latpix])
                    
                    hold on
                    drawnow
                end
                
                for i=1:Nx % this is position loop
                    %%% here if contains:
                    %%% 1) time values are not NaN, (outside of the clicked points were assigned to be NaN)
                    %%% 2) position values are not NaN, (outside of the clicked points were assigned to be NaN)
                    %%% 3) if time value is greater than or equal to k,
                    %%% 3) if time value is smaller than k,
                    if ((~isnan(track(i,1)))&&(~isnan(track(i,2)))&&(track(i,1)>=k)&&(track(i,1)<k+1))
                        LS((k-indexminint+1),1)=LS((k-indexminint+1),1)+track(i,1);
                        LS((k-indexminint+1),2)=LS((k-indexminint+1),2)+track(i,2);
                        LS((k-indexminint+1),3)=LS((k-indexminint+1),3)+1;
                    end
                    
                    
                end
            end
        end
        
        %%% this loop now divides by the count to get the average
        for i=1:size(LS,1) % this is position loop
            if (LS(i,3)~=0)
                LS(i,1)=LS(i,1)/LS(i,3);
                LS(i,2)=LS(i,2)/LS(i,3);
            else
                LS(i,1)=NaN;
                LS(i,2)=NaN;
            end
        end
        
        
        %%% calculate SEM
        
        stdevsumLS=zeros((indexmaxint-indexminint+1),3);
        for k=indexminint:indexmaxint
            %clc
            %fprintf('Kymo %d of %d: %f percent completed\n',FirstLSno,NFirstLS,100*(k-indexminint)/(indexmaxint-1-indexminint))
            for lsno=1:Nls % this is linescan loop
                clear tracknotflipped track
                FileNames=sprintf('%s/%s',foldername,Files(lsno).name);
                fidr=fopen(FileNames,'r');
                %             [filepath,name,ext] = fileparts(FileNames);
                %             newStr = extractAfter(name,'K');
                %             kymonoextracted = extractBefore(newStr,'_');
                %             kymono=str2num(kymonoextracted);
                %             fprintf('%s\t%s\t%d\n',name,kymonoextracted,kymono)
                tracknotflipped=fscanf(fidr,'%f %f',[2 Inf]);
                fclose(fidr);
                tracknotflipped=tracknotflipped';
                
                
                track=zeros(size(tracknotflipped,1),2);
                for i=1:size(tracknotflipped,1)
                    track(i,1)=(-1.0)*tracknotflipped(i,1);
                    track(i,2)=tracknotflipped(i,2);
                end
                
                Nx=size(track,1);
                
                for i=1:Nx % this is position loop
                    %%% here if contains:
                    %%% 1) time values are not NaN, (outside of the clicked points were assigned to be NaN)
                    %%% 2) position values are not NaN, (outside of the clicked points were assigned to be NaN)
                    %%% 3) if time value is greater than or equal to k,
                    %%% 3) if time value is smaller than k,
                    if ((~isnan(track(i,1)))&&(~isnan(track(i,2)))&&(track(i,1)>=k)&&(track(i,1)<k+1))
                        stdevsumLS((k-indexminint+1),1)=stdevsumLS((k-indexminint+1),1)+ track(i,1);
                        stdevsumLS((k-indexminint+1),2)=stdevsumLS((k-indexminint+1),2)+ ( (track(i,2)-LS((k-indexminint+1),2))*(track(i,2)-LS((k-indexminint+1),2)) );
                        stdevsumLS((k-indexminint+1),3)=stdevsumLS((k-indexminint+1),3)+ 1;
                    end
                    
                    
                end
            end
        end
        for i=1:size(stdevsumLS,1) % this is position loop
            
            % internal check
            if (stdevsumLS(i,3)~=LS(i,3))
                fprintf('count doesnt match for i=%d\n',i)
            end
            
            % internal check
            if ( (stdevsumLS(i,3)~=0) && (stdevsumLS(i,1)/stdevsumLS(i,3))~=LS(i,1) ) % since LS is already divided
                fprintf('position doesnt match for i=%d\t%f\t%f\n',i,stdevsumLS(i,1)/stdevsumLS(i,3),LS(i,1))
            end
            
            if (stdevsumLS(i,3)~=0)
                if (stdevsumLS(i,3)>1)
                    LS(i,4)=sqrt( stdevsumLS(i,2)/(stdevsumLS(i,3)-1) ); % This is standard deviation
                    LS(i,4)=LS(i,4)/sqrt(stdevsumLS(i,3)); % This is SEM
                elseif (stdevsumLS(i,3)==1)
                    LS(i,4)=0;
                end
            else
                LS(i,4)=NaN;
            end
        end
        
        
        %%% finish calculate SEM
        
        %%% remove any NaN values that might be assigned during the
        %%% division above
        LS1 = LS(any(LS,2),:);
        LS=LS1;
        
        
        aveI=zeros(size(LS,1),4);
        for i=1:size(LS,1)
            aveI(i,1)=LS(i,1);
            aveI(i,2)=LS(i,2);
            aveI(i,3)=LS(i,3);
            aveI(i,4)=LS(i,4);
        end
        
        
        
        
        
        
        %%% plot the average curve
        set(0, 'CurrentFigure', f1)
        %plot(aveI(:,1),aveI(:,2),'LineWidth',5,'Color','g')
        errorbar(aveI(:,1),aveI(:,2),aveI(:,4),'LineWidth',5,'Color','g')
        title('lattice bg not substracted')
        hold off
        drawnow
        saveas(f1,sprintf('%s-Average.png',kymoprefix))
        
        
        
        % for double-checking, if the max intensity pixel has zero x-value
        [maxval,index]=max(aveI(:,2)); %returns index of maximum values in vector index
        
        Blat=mean(aveI(end-5:end,2));
        
        
        %%% now subtract average Blat and save
        f2=figure;
        for lsno=1:Nls
            clearvars -except foldername FilesFirstLS NFirstLS FirstLSno Files Nls lsno fnaveraged Blat Bsol f2 kymoprefix aveI solpix latpix index indexminint indexmaxint
            clear tracknotflipped track
            FileNames=sprintf('%s/%s',foldername,Files(lsno).name);
            
            fidr=fopen(FileNames,'r');
            tracknotflipped=fscanf(fidr,'%f %f',[2 Inf]);
            fclose(fidr);
            tracknotflipped=tracknotflipped';
            
            
            FileNames=sprintf('%s/BlatSubtracted-%s',foldername,Files(lsno).name);
            fidw=fopen(FileNames,'w');
            
            track=zeros(size(tracknotflipped,1),2);
            zeroindex=NaN;
            for i=1:size(tracknotflipped,1)
                track(i,1)=(-1.0)*tracknotflipped(i,1);
                track(i,2)=tracknotflipped(i,2)-Blat;
                if (track(i,1)==0)
                    zeroindex=i;
                end
            end
            
            
            writestart=zeroindex-indexmaxint; % since it's flipped
            writeend=zeroindex-indexminint; % since it's flipped
            
            if (writestart<1)
                writestart=1;
            end
            if (writeend>size(track,1))
                writeend=size(track,1);
            end
            
            for i=writestart:writeend
                fprintf(fidw,'%f\t%f\n',track(i,1),track(i,2));
            end
            
            fclose(fidw);
            
            
            set(0, 'CurrentFigure', f2)
            plot(track(:,1),track(:,2))
            xlim([solpix latpix])
            hold on
            drawnow
            
        end
        
%         if (aveI(index,1)~=0)
%             fprintf('max intensity pixel is not at zero.\n');% Terminating...\n');
%             fidnotsaved=fopen(sprintf('%s/NotBGsubtractedList.dat',foldername),'a');
%             fprintf(fidnotsaved,'%s\n',kymoprefix);
%             %return;
%         else
            savename=sprintf('%s-BGsubstracted-Average.dat',kymoprefix);
            fidw=fopen(savename,'w');
            for i=1:size(aveI,1)
                aveI(i,2)=aveI(i,2)-Blat;
                
                fprintf(fidw,'%f\t%f\t%f\t%f\n',aveI(i,1),aveI(i,2),aveI(i,3),aveI(i,4));
            end
            fclose(fidw);
%         end
        
        set(0, 'CurrentFigure', f2)
        %plot(aveI(:,1),aveI(:,2),'LineWidth',5,'Color','g')
        errorbar(aveI(:,1),aveI(:,2),aveI(:,4),'LineWidth',5,'Color','g')
        hold on
        %xlim([-15 25])
        %ylim([-50 200])
        title('lattice bg subtracted')
        hold off
        saveas(f2,sprintf('%s-BGsubstracted.png',kymoprefix))
        exportgraphics(f2,sprintf('%s-BGsubstracted.pdf',kymoprefix))
        fprintf('Saved %s-BGsubstracted-Average.dat\n',kymoprefix)
        
    end
    
    
    
    
catch mefn
    fprintf('Error running the function:\n%s\n',mefn.message)
end
end