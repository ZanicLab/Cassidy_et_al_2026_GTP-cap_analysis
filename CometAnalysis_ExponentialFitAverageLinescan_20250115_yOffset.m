%%% Created by Goker Arpag on May 6, 2018
%%% Last edited on: 2024-05-21 by AC
%%% Read a data file that contains the average intensity profile of the
%%% comets and fit to an exponential decay convolved with a Gaussian.
%%% SW edits: sigma and background fixed, A constrained to be
%%% between .65 and 1.05. x0 and lambda are free
%%% parameters. use if statements and for loops to determine what Blat
%%% and Bsol should be. plots entire linescan but r-squared/fits depend on
%%% only the 20 pixels. for 70 pixels.
%%% usage: CometAnalysis_Fit_v1('folder','filename')
%%% i.e.: CometAnalysis_Fit_v1('./.','20uMCy5_200nMEB1GFP_005_488_K-1-average.dat')

% % % Example useage for superaveraged linescan: CometAnalysis_ConvolutionFitAverageLinescan('R1_20uM_P/individualLS','superaveraged.dat',15,20)
% % % Example usage for individual average from one growth event: CometAnalysis_ConvolutionFitAverageLinescan('R1_20uM_P/individualLS','K0011-BGsubstracted-Average.dat',15,20)

function CometAnalysis_ExponentialFitAverageLinescan_20250115_yOffset(path,subfolder,filename,outputfilename,fitstartpix,fitendpix)
try
    close all
    %%%pixel size MUST be updated 
    pixelsize=158.3;
    
    f1=figure;
    inputname=sprintf('%s/%s/%s',path,subfolder,filename);
    [filepath,name,ext] = fileparts(filename);
    aveI=dlmread(inputname); %reads numeric data from the ASCII delimited file
        
    
    [maxval,index]=max(aveI(:,2)); %returns index of maximum values in vector index

        aveI(:,1)=aveI(:,1)-fitstartpix;

        %[maxval,index]=max(aveI(:,2));
        
        fitrange=(aveI(:,1)>=0)&(aveI(:,1)<=fitendpix)
        
        set(0, 'CurrentFigure', f1)
        plot(aveI(:,1)+fitstartpix, aveI(:,2),'ko')
        hold on
        plot(aveI(fitrange,1)+fitstartpix, aveI(fitrange,2),'go')
        hold on
        

        
        %ratio = (maxval - B);
        %%% fit to an exponential decay convolved with a gaussian
        %below, x0 is a free parameter, B and sigma are fixed, A is constrained, free parameters are lambda, and xzero (for now)
        % Define the new fit type with an offset (B)
        ft = fittype('A*exp(-x/lambda) + B', 'coefficients', {'A', 'lambda', 'B'});
        try
            % Perform the fit with an initial guess for A, lambda, and B
            [fo, gof] = fit(aveI(fitrange,1), aveI(fitrange,2), ft, 'StartPoint', [maxval 1 mean(aveI(fitrange,2))], 'Lower', [0 0 -Inf], 'Upper', [Inf Inf Inf]);
        catch ME
            fprintf('Cannot fit function: %s\n',ME.message);
        end
        set(0, 'CurrentFigure', f1)
        %plot(fo)

        xfitresult = linspace(0, 25);
        yfitresult = fo(xfitresult);
        plot(xfitresult+fitstartpix, yfitresult)

        hold off
        %xlim([-2 25])
        %ylim([-10 50])
        %saveas(f1,sprintf('%s/%s/ExponentialFitted_StartPoint%d-%s.png',path,subfolder,fitstartpix,name));
        pbaspect([1 1 1])
        exportgraphics(f1,sprintf('%s/%s/ExponentialFitted_StartPoint%d-%s.pdf',path,subfolder,fitstartpix,name),'Resolution',500);
        
        %%% get fit coefficients
        Coeffs = coeffvalues(fo);
        Afit=Coeffs(1);
        lambdafit=Coeffs(2);
        lambdafitnm=Coeffs(2) * pixelsize;
        Bfit = Coeffs(3); % New coefficient for the offset
        
        ci = confint(fo);
        Aci=ci(:,1);
        lambdaci=ci(:,2) * pixelsize;
        Bci = ci(:,3); % Confidence intervals for the offset
        
        lambdacirange=(lambdaci(2)-lambdaci(1))/lambdafitnm;
        
%%% individual fit
fid=fopen(sprintf('%s/%s',path,outputfilename),'a');
%kymoindex = extractBetween(filename,'Sorted-ALL-','.dat');
kymoindex = extractBetween(filename,'Selected-','-BGsubstracted-Average.dat');
%kymoindex = extractBefore(filename,'-BGsubstracted-Average.dat');
%if ( ~isnan(lambdaci(1)) && ~isnan(lambdaci(2)) && ~isnan(Aci(1)) && ~isnan(Aci(2)) && (gof.rsquare>=0) && (gof.adjrsquare>=0) )

fprintf(fid, ['%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.5f\t%.5f\t%f\t%s\n'], ...
    lambdafitnm, lambdaci(2)-lambdafitnm, lambdafitnm-lambdaci(1), ...
    Afit, Aci(2)-Afit, Afit-Aci(1), Bfit, Bci(2)-Bfit, Bfit-Bci(1), ...
    gof.rsquare, gof.adjrsquare, gof.rmse, string(kymoindex))
        %fprintf(fid,'%s\t%f\t[%f - %f]\n',filename,lambdafit,round(lambdaci(1)),round(lambdaci(2)));
        fclose(fid);
    %end
    

    
catch mefn
    fprintf('%s\n',mefn.message)
end
end