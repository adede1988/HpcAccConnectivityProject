



output = zeros(100,2); 

for ii = 1:100

    N = 2; 
    
    x = randsample([-pi:.0001:pi], N, true); 
    y = randsample([-pi:.0001:pi], N, true);
    
    difs = x - y; 
    
    %ISPC calculation
    output(ii,1) = abs(mean(exp(1i * (difs))));
    
    %PPC calculation
    output(ii,2) = ...
        mean(...
            cell2mat(arrayfun(@(j) ...
                cell2mat(arrayfun(@(k) ...
                    cos(difs(j))*cos(difs(k)) + sin(difs(j))*sin(difs(k)), ...
                j+1:N, 'uniformoutput', false)), ...
            1:N-1, 'uniformoutput', false))...
        );



end

disp(['ispc: ' num2str(round(ispc,2)) ' ppc: ' num2str(round(ppc,2))])

