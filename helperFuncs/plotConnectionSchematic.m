function [] = plotConnectionSchematic(plotMat, frex, plotOpt, regions, ...
    keyRegIdx, colOpt, regColors, wantLegend, fn1, fn2, fn3)

s = [85,25,86];
m = [186,21,77]; 
e = [249,205,15]; 

s2w2y = [[linspace(s(1),m(1),128)'; linspace(m(1),e(1),128)'], ...
         [linspace(s(2),m(2),128)'; linspace(m(2),e(2),128)'], ...
         [linspace(s(3),m(3),128)'; linspace(m(3),e(3),128)'], ...
         ] / 255;
test = s2w2y; 
s2w2y(1:85,:) = repmat([96,62,149]/255, [85,1]);
s2w2y(86:171,:) = repmat([0,157,161]/255, [86,1]);
s2w2y(171:256,:) = repmat([249,205,15]/255, [86,1]);


globalSize = 150; 
%acc, dlpfc, hip, mtl, ppfc (order as in keyRegIdx)
reg_x = [0, .35, -2.15, -1.25, 1.35]; 
reg_y = [0, 1, -.65, -.9, -1]; 


rad = .15; 
rad2 = 3.5; 
circAngles = [2.3, 1.3, 3.2, 4.6, 0];
figure('visible', false, 'position', [0,0,800,500]) %try to draw one connection at a time
img = imread('R:\MSS\Johnson_Lab\dtf8829\publicationFigureData\brain.png');
img(img==0) = 255; 
hold off
image([-4.15, 1.85],[2.15, -2.75], img)
hold on 
scatter(reg_x, reg_y, 600, regColors(keyRegIdx, :),  'filled')
xlim([-2.65, 1.85])
ylim([-1.95, 1.35])
set(gca, 'ydir', 'normal')
hold on 

alreadyPlotted = zeros(5); 
for ii = 1:5
    for jj = 1:5
        skip = false; 
        %get the single values for freq, hit, miss, t
        %get a mask of p values to use to screen for significance
        pMask = squeeze(plotMat(ii, jj, :, :, 4)); 
        if sum(pMask<.05, 'all') >0
           
             
        %get mean hit value
        hitMat = squeeze(plotMat(ii, jj, :,:, 1)); 
        meanH = mean(hitMat(pMask<.05), 'all'); 
        if meanH < 0 
            meanH = 0; 
        end
        meanH = sqrt(meanH); 

        %get mean miss value
        missMat = squeeze(plotMat(ii, jj, :,:, 2)); 
        meanM = mean(missMat(pMask<.05), 'all'); 
        if meanM < 0 
            meanM = 0; 
        end
        meanM = sqrt(meanM); 

        %get mean t value
        tMat = squeeze(plotMat(ii, jj, :,:, 3)); 
        meanT = mean(tMat(pMask<.05), 'all'); 
        if meanT < 0 
            meanT = 0; 
            disp(['negative connection:  ' ...
                regions{keyRegIdx(ii)} ' ' regions{keyRegIdx(jj)}])
            skip = true; 
        end
       
        %get mean freq
        sigVal = sum(pMask<.05); 
        meanF = sum(frex.*sigVal) / sum(sigVal);

        switch plotOpt %what determines thickness? 

            case "hitVal"
                scaleFact = meanH; 
            case "missVal"
                scaleFact = meanM; 
            case "tVal"
                scaleFact = meanT^4 / 1000; 
            otherwise
                scaleFact = 0; 
        end

        switch colOpt %what determines color? 

            case "freq"
                indices = logspace(log10(frex(1)), log10(frex(12)), 256);
                coli = find(indices>=meanF, 1);
                if isempty(coli)
                    coli = 255; 
                end
                plotCol = s2w2y(coli, :); 
            case "reg"
                plotCol = regColors(keyRegIdx(ii), :); 
            otherwise
                scaleFact = 0; 
        end

        if ~skip 
            if alreadyPlotted(ii,jj) == 1 && strcmp(colOpt, 'freq')
                'skipping! '
            else

        if ii == jj
            %draw a line to the orig point, and establish its angle
            theta = circAngles(ii);
           
            xC = reg_x(ii) + rad*cos(theta); 
            yC = reg_y(ii) + rad*sin(theta); 
            
            theta = linspace(-pi, pi, 1000); 

            scatter(xC + rad*cos(theta), yC + rad*sin(theta), ...
                50+(scaleFact * globalSize)^2, plotCol, 'filled')

        else
        
            if ii==1 || jj==1 %don't make arcs involving acc
               
                plx = linspace(reg_x(ii), reg_x(jj), 1000); 
                ply = linspace(reg_y(ii), reg_y(jj), 1000); 
                if alreadyPlotted(ii,jj) == 1 && ...
                                strcmp(colOpt, "reg")
                    d = sqrt((reg_x(ii) - reg_x(jj))^2 + ...
                             (reg_y(ii) - reg_y(jj))^2); 
                    elim = []; 
                    %use to make dashed lines if a connection 
                    % is already plotted
                    for ei = 1:1000
                        if mod(round(ei/((scaleFact*1500)/d)), 4)
                            elim = [elim, ei]; 
                        end
                    end
                    plx(elim) = []; 
                    ply(elim) = []; 
                end
                scatter(plx, ply, ...
                50+(scaleFact * globalSize)^2, plotCol, 'filled')
               
            else
                %find the center point of the circle including both points
                
                M = [(reg_x(ii) + reg_x(jj))/2 ,(reg_y(ii) + reg_y(jj))/2];
                if (ii == 5 && jj == 4) ||(jj == 5 && ii == 4)
                    rad2 = sqrt((reg_x(ii) - reg_x(jj))^2 + ...
                        (reg_y(ii) - reg_y(jj))^2) *1;
                else
                    rad2 = sqrt((reg_x(ii) - reg_x(jj))^2 + ...
                        (reg_y(ii) - reg_y(jj))^2) *.6;
                end
                s = (reg_y(ii) - reg_y(jj)) / (reg_x(ii) - reg_x(jj));
                sperp = -1/s; %perp bisect slope
                intercept = M(2) - sperp*M(1); %perp bisect intercept
                
                %arms of a right triangle with hypotenus equal to radius
                %a = half of the line through the two points
                %b = from midpoint between points to center of circle
                a = sqrt((M(1) - reg_x(ii))^2 + (M(2) - reg_y(ii))^2); 
                b = sqrt(rad2^2 - a^2);
% %                 maxy = max([reg_y(ii), reg_y(jj)]);
%                 %make a split line for which direction to make the arc
%                 %bulge, based on phg(4) point and just below the dlpfc(2) 
%                 u = reg_y(2) - .5; 
%                 bs = (reg_y(4)*reg_x(2) -u*reg_x(4)) / ...
%                     (reg_x(2) - reg_x(4)); 
%                 ms = (u - bs) / reg_x(2);  
% 
% %                 if M(2) > M(1)*ms + bs
                    %start at M and walk positive
                    testX = [M(1):.01:2]; 
                    testY = testX*sperp + intercept; 
                    dists = sqrt((M(1) - testX).^2 + (M(2) - testY).^2);

                    disti = find(dists - b > 0, 1);
                    xC1 = testX(disti); 
                    yC1 = testY(disti); 
% 
%                 
% %                 else 
                    %start at negative and walk to M
                    testX = [-3:.01:M(1)]; 
                    testY = testX*sperp + intercept; 
                    dists = sqrt((M(1) - testX).^2 + (M(2) - testY).^2);

                    disti = find(dists - b < 0, 1);
                    xC2 = testX(disti); 
                    yC2 = testY(disti); 
% 
% %                 end
                  
                  %look for circle centers in both directions: 
                  %using quadratic equation: 
%                   a = 1 + sperp^2; 
%                   c = -(b^2 - M(1)^2 - intercept^2 + 2*M(2)*intercept - M(2)^2);
%                   b = -2*M(1) + 2*sperp*intercept + 2*sperp;
% 
%                   xC1 = (-b + sqrt(b^2 - 4*a*c)) / 2*a; 
%                   xC2 = (-b - sqrt(b^2 - 4*a*c)) / 2*a; 
% 
%                   yC1 = xC1*sperp + intercept; 
%                   yC2 = xC2*sperp + intercept; 
                  

                  %find if xC2, yC2 or xC1, yC1 is closer to origin, then
                  %use other. 
                  if sqrt(xC1^2 + yC1^2) <  sqrt(xC2^2 + yC2^2)
                    xC = xC1; 
                    yC = yC1; 
                  else
                      xC = xC2; 
                      yC = yC2; 
                  end


                    theta = linspace(-pi, pi, 10000); %[-pi:pi/100:pi];
                    arcX = xC+rad2*cos(theta); 
                    arcY = yC+rad2*sin(theta); 

                    theta1 = atan2(reg_y(ii) - yC, reg_x(ii) - xC);
                    theta2 = atan2(reg_y(jj) - yC, reg_x(jj) - xC);
                    
                    theta_points = atan2(arcY - yC, arcX - xC);
                    if theta1 > theta2
                        points_between = (theta_points <= theta1) & ... 
                                        (theta_points >= theta2);
                    else
                        points_between = (theta_points >= theta1) & ... 
                                        (theta_points <= theta2);
                    end
                    
                    if sum(points_between) < 5000
                        points_between = find(points_between);
                        if alreadyPlotted(ii,jj) == 1 && ...
                                strcmp(colOpt, "reg")
                            L = length(points_between);
                            elim2 = []; 
                            for ei = 1:L
                                if mod(round(ei/(L/5)), 2)
                                    elim2 = [elim2, ei]; 
                                end
                            end
                            points_between(elim2) = []; 
                        end
                        scatter(arcX(points_between), ...
                            arcY(points_between), ...
                            50+(scaleFact * globalSize)^2, plotCol, 'filled')
                    else
                        points_between = find(~points_between);
                        if alreadyPlotted(ii,jj) == 1 && ...
                                strcmp(colOpt, "reg")
                            L = length(points_between);
                            elim2 = []; 
                            for ei = 1:L
                                if mod(round(ei/(L/10)), 2)
                                    elim2 = [elim2, ei]; 
                                end
                            end
                            points_between(elim2) = []; 
                        end
                        scatter(arcX(points_between), ...
                            arcY(points_between), ...
                            50+(scaleFact * globalSize)^2, plotCol, 'filled')
                    end



            end

        end
        
        alreadyPlotted(jj, ii) = 1;
            end
        end
        end
    end
end

scatter(reg_x, reg_y, 1000, regColors(keyRegIdx, :),  'filled')
% text(reg_x, reg_y, regions(keyRegIdx))
set(gcf,'color','w');
box off;
set(gca, 'color', 'none');
% set(gca,'XColor', 'none','YColor','none')
export_fig(fn1, '-r300')


if wantLegend

   
        figure('visible', false, 'position', [0,0,800,600]) 
        ylim([-1.85, 1.35])
        xlim([-.35, 1.2])
        yVals = linspace(-1.55, 1.25, 7); 
        hold on 
        frexKey = logspace(log10(frex(1)), log10(frex(12)), size(s2w2y,1));
        idx = round(linspace(1, size(s2w2y,1), 7)); 
        for ii = 1:length(yVals)
            plx = [0:.001:1];
            ply = yVals(ii)*ones(length(plx), 1); 
            scatter(plx, ply, ...
                        50+(3.5.^4 ./ 1000 * globalSize)^2,...
                        s2w2y(idx(ii),:), 'filled')
        end
        yticks(yVals); 
        yticklabels(round(frexKey(idx), 2))
        set(gcf,'color','w');
        box off;
        export_fig(fn2, '-r300')


  

    switch plotOpt %what determines thickness? 

        case "hitVal"
            scaleFact = [.01:.03:.21]; 
        case "missVal"
            scaleFact = [.01:.03:.21]; 
        case "tVal"
            scaleFact = linspace(2,4,7).^4 ./ 1000; 
        otherwise
            scaleFact = 0; 
    end

    
    figure('visible', false, 'position', [0,0,800,600]) 
    ylim([-1.65, 1.45])
    xlim([-.35, 1.2])
    yVals = linspace(-1.55, 1.25, length(scaleFact)); 
    hold on 
    for ii = 1:length(scaleFact)
        plx = [0:.001:1];
        ply = yVals(ii)*ones(length(plx), 1); 
        scatter(plx, ply, ...
                    50+(scaleFact(ii) * globalSize)^2, 'k', 'filled')
    end
    yticks(yVals); 
    yticklabels(round(linspace(2,4,7), 1))
    set(gcf,'color','w');
    box off;
    export_fig(fn3, '-r300') 
    
    
   
end

end