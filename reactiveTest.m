function [reactive] = reactiveTest(HFB)

reactive = zeros(4,1); 
cnd = fieldnames(HFB); 

 %enc onset
    cc = [1 2]; 
    comboHFB = []; 
    for con = 1:length(cc)
        comboHFB = [comboHFB HFB.(cnd{cc(con)})];
    end

    test = mean(comboHFB,2); 
    test = checkForThreshold(test, HFB.encMulTim, [-50,2000] ); 
    reactive(1) = test;

    %enc RT
    cc = [5 6]; 
    comboHFB = []; 
    for con = 1:length(cc)
        comboHFB = [comboHFB HFB.(cnd{cc(con)})];
    end

    test = mean(comboHFB,2); 
    test = checkForThreshold(test, HFB.encMulTim, [-1500,500] ); 
    reactive(2) = test;

    %ret on
    cc = [9 10 11 12]; 
    comboHFB = []; 
    for con = 1:length(cc)
        comboHFB = [comboHFB HFB.(cnd{cc(con)})];
    end

    test = mean(comboHFB,2); 
    test = checkForThreshold(test, HFB.encMulTim, [-50,2000] ); 
    reactive(3) = test;

    %ret rt
    cc = [15 16 17 18]; 
    comboHFB = []; 
    for con = 1:length(cc)
        comboHFB = [comboHFB HFB.(cnd{cc(con)})];
    end

    test = mean(comboHFB,2); 
    test = checkForThreshold(test, HFB.encMulTim, [-1500,500] ); 
    reactive(4) = test;


end