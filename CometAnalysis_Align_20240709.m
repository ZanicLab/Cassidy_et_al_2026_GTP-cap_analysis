% % %
% % % Created by GA on 2018-07-04
% % % Last edited by GA on 2021-04-15

% % % This code only aligns the linescans at the max intensity:
% % % It reads the excel file of clicked points and the kymograph file. The
% % % comet location is estimated using the clicked points from excel file.
% % % Then the pixel with max intensity is found and its x-coordinate is
% % % placed at zero. After aligning each linescan, the mean intensity and
% % % its standard deviation are calculated. The linescans that have max
% % % intensity lower than mean intensity - standard deviation are discarded.
% % % At this point, a linear fit is performed using not aligned tip
% % % location to estimate the growth velocity. Outliers are determined if
% % % the residual of a given point is outside of meanresidual ±
% % % 1.5*stdevresidual. After removing velocity outliers each linescan is
% % % saved as individual data file. If velocity will be binned, there will
% % % be separate folders to save the data. If velocity will not be binned,
% % % there will be only one folder to save the data.

%line 173 - change points on velocity plot
%line 208 - change points on kymo overlay

% % % Example usage CometAnalysis_Align_20210429('/Users/zanicimac/Desktop/Personal_Folders/Todd/Control/','R1_20uM_P','Points.xlsx',5,160,2)

function CometAnalysis_Align_20240709(pathtofile,outputfoldersuffix,filename,windowsize,pospix,timepix)
%try

close all
clearvars -except pathtofile outputfoldersuffix filename windowsize pospix timepix binwidth
ft='A + (x*B)';
foldername=sprintf('%s/individualLS_%s',pathtofile,outputfoldersuffix);
if exist(foldername, 'dir')
    rmdir(foldername, 's')
end
mkdir(foldername);

%%% velocity data output
velocityname=sprintf('%s/fitted-velocities_%s.dat',pathtofile,outputfoldersuffix);
fidvelocity=fopen(velocityname,'w');

%%% Read the excel file that contains kymograph info
T=readtable(sprintf('%s/%s',pathtofile,filename));
M=table2struct(T);

NK=size(M,1); % 2 times number of kymographs, since each kymograph has two lines of info
for kymono=1:2:NK % loop over kymos, skip 2 since each kymograph has two lines of info
    clearvars -except pathtofile filename windowsize M NK kymono pospix timepix ft fidvelocity foldername
    close all
    forfitcounter=0;
    Nincluded=0;
    
    %%% read kymograph
    fprintf('\tworking on %s\t%d\n',M(kymono).FileName,kymono)
    kymoname=sprintf('%s',M(kymono).FileName);
    [K,map]=imread(sprintf('%s/%s',pathtofile,kymoname));
    
%     %%% Save the image that is being used, just to double check
%     f0=figure(1000+(2*kymono)-1);
%     imshow(K,map)
%     saveas(f0,sprintf('%s/%s-imageread.png',pathtofile,M(kymono).FileName))
%     close(1000+(2*kymono)-1);
    %%% get size of the kymograph
    [Nt,Nx]=size(K);
    
    include=zeros(Nt,1);
    maxintensityvaluesattip=NaN(Nt,1);
    forfit=NaN(Nt,2);
    %%% read the coordinates of the clicked points on the kymograph
    %%% this will be used to estimate the location of the comet
    %%% again, in our kymographs x is x, y is time
    xi=floor(M(kymono+0).X)+1;
    xf=floor(M(kymono+1).X)+1;
    ti=floor(M(kymono+0).Y)+1;
    tf=floor(M(kymono+1).Y)+1;
    
    if ( (xi>=5) && ((xf<=(Nx-10))) && (ti>0) && (tf<=Nt) )
        
        %%% calculate the slope to estimate comet location
        %%% this is basically the growth velocity of the kymograph
        slope=(xf-xi)/(tf-ti);
        
        
        %%%
        %%% At this point, reading and initializing has finished
        %%%
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% pixel level alignment based on maximum intensity %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        LSallaligned=zeros(Nx,2*Nt);
        for j=1:Nt % this is time loop
            
            %%% skip if the line is outside of clicked points
            if (( j<ti )||( j>tf ))
                LSallaligned(:,((2*j)-1))=NaN;
                LSallaligned(:,((2*j)-0))=NaN;
                continue; % this is j loop
            end
            %%% skip of close to sides of the kymo
            
            
            %%% estimate comet location: x=x0 + (v*t)
            t=j-ti; % this is deltaT
            ep= floor( xi + (slope*t) ); % slope was deltaX/deltaT, then deltaX=x-xi=slope*deltaT=x=xi+(slope*deltaT)
            
            %%% estimate lower and upper bounds for finding local maximum
            %%% for the comet
            %%% lower bound for comet tip location during fit
            if ((ep-windowsize)>1)
                lower=ep-windowsize;
            else
                lower=1;
            end
            %%% upper bound for comet tip location during fit
            if ((ep+windowsize)<(Nx-1))
                upper=ep+windowsize;
            else
                upper=Nx;
            end
            
            
            %%% pixel level alignment
            [maxval,index]=max(K(j,lower:upper));
            index=index+lower-1;
            
            for i=1:Nx % this is position loop
                LSallaligned(i,((2*j)-1))=i-index; % put the pixel with max intensity at x=0
                LSallaligned(i,((2*j)-0))=K(j,i);
            end
            %%% pixel level alignment finished
            
            
            forfitcounter=forfitcounter+1;
            forfit(j,1)=j*timepix; % time
            forfit(j,2)=index*pospix; % position
            maxintensityvaluesattip(j,1)=LSallaligned(index,((2*j)-0));
            Nincluded=Nincluded+1;
            include(j,1)=1;
            
        end % end time loop for pixel level alignment
        
        forfitplot=forfit; % to show the ones that removed with both intensity and velocity filtering
        
        for j=1:Nt
            if ( isnan(maxintensityvaluesattip(j,1)) )
                include(j,1)=0;
                forfit(j,1)=NaN;
                forfit(j,2)=NaN;
            end
        end
        
        %%% remove all NaN rows
        forfit1 = forfit(any(forfit,2),:);
        forfit=forfit1;
        
        if (all(isnan(LSallaligned(:))))
            fprintf('Every line is skipped for %s. Going to next kymograph.\n',kymoname);
            continue;
        else
            
            
            % determine velocity
            [fo, gof , fitoutput] = fit(forfit(:,1),forfit(:,2), ft,'StartPoint', [1 1]);
            TFresidualoutlier=isoutlier(fitoutput.residuals);
            
            if (gof.rsquare>=0)
                f0=figure(2000+(2*kymono)-1);
                plot(forfit(~TFresidualoutlier,1),forfit(~TFresidualoutlier,2),'go', 'MarkerSize', 5, "MarkerFaceColor", 'g')
                hold on
                plot(forfit(TFresidualoutlier,1),forfit(TFresidualoutlier,2),'bo')
                hold on
                plot(fo)
                hold off
                title(sprintf('%s',kymoname))
                exportgraphics(f0,sprintf('%s/VelocityFit-K%4.4d.pdf',foldername,kymono));
                close(2000+(2*kymono)-1);
                
                Coeffs = coeffvalues(fo);
                Intercept=Coeffs(1);
                Velocity=Coeffs(2);
                ci = confint(fo);
                Interceptci=ci(:,1);
                Velocityci=ci(:,2);
                fprintf(fidvelocity,'%f\t%f\t%f\t%f\t%s\tK%4.4d\t%d\n',Velocity,Velocityci(1),Velocityci(2),gof.rsquare,kymoname,kymono,sum(TFresidualoutlier));
                
                %%% save each linescan
                figure(2000+(2*kymono)-1)
                imshow(K,map)
                hold on
                linescancounter=0;
                for j=ti:tf
                    if (include(j,1)==1)
                        linescancounter=linescancounter+1;
                        lswritename=sprintf('%s/K%4.4d_LS%4.4d.dat',foldername,kymono,linescancounter);
                        %fprintf('%s\n',lswritename)
                        fidls=fopen(lswritename,'w');
                        
                        for iii=1:size(LSallaligned,1)
                            fprintf(fidls,'%f\t%f\n',LSallaligned(iii,((2*j)-1)),LSallaligned(iii,((2*j)-0)));
                        end
                        fclose(fidls);
                        figure(2000+(2*kymono)-1)
                        plot(forfitplot(j,2)/pospix,j,'go', "MarkerSize", 2)
                        hold on
                    else
                        figure(2000+(2*kymono)-1)
                        plot(forfitplot(j,2)/pospix,j,'ro')
                        hold on
                    end
                end
                f0=figure(2000+(2*kymono)-1);
                plot(xi,ti,'b*') % x on in x-axis, t is on y-axis in the kymograph
                hold on
                plot(xf,tf,'b*') % x on in x-axis, t is on y-axis in the kymograph
                hold on
                title(sprintf('%s',kymoname))
                exportgraphics(f0,sprintf('%s/K%4.4d.pdf',foldername,kymono));
                fprintf('Saved %s\t%d\n',M(kymono).FileName,kymono)
            else
                fprintf('Bad velocity fit. Skipped %s\t%d\n',M(kymono).FileName,kymono)
            end
        end
    else
       fprintf('Too close to edge. Skipped %s\t%d\n',M(kymono).FileName,kymono)
    end
    
end %loop over kymos

fclose(fidvelocity);

% catch mefn
%     fprintf('Error running the function:\n%s\n',mefn.message)
% end
end