function [allData] = getAllDataStruct(datFolder, masterSheet, task)




datDirs = dir(datFolder);
datDirs = datDirs([datDirs.isdir]==true);
datDirs = datDirs(3:end); 

allData = struct; 
si = 1; 
%loop on data collection sites
for ii = 1:length(datDirs)
    curDir = dir([datDirs(ii).folder '\'  datDirs(ii).name]);
    curDir = curDir([curDir.isdir]==true);
    curDir = curDir(3:end); 
    %loop on subjects
    for jj = 1:length(curDir)
        %check if subject is ready
        masteri = find(cellfun(@(x) ~isempty(x), ...
            (cellfun(@(x) find(strcmp(curDir(jj).name, x)), masterSheet.subID, 'UniformOutput', false))));

        subDir = dir([curDir(jj).folder '\' curDir(jj).name '\Tasks']);
        %check if the subject has done the target task && has ready data
        if sum(cellfun(@(x) strcmp(task, x), {subDir.name}))==1 && strcmp('ready', masterSheet.iEEG(masteri))
            ti = find(strcmp(task, {subDir.name}));

            taskDir = dir([subDir(ti).folder '\'  subDir(ti).name '\BPM']);
            if isempty(taskDir)
                taskDir = dir([subDir(ti).folder '\'  subDir(ti).name '\BPR']);
            end
            

            allData(si).site = datDirs(ii).name; 
            allData(si).subID = curDir(jj).name; 
            allData(si).dataDir = taskDir(1).folder; 
            ei = find(cellfun(@(x) ~isempty(x), strfind({taskDir.name},{'enc'}) ) );
            allData(si).encDatFn = taskDir(ei).name;
            ri = find(cellfun(@(x) ~isempty(x), strfind({taskDir.name},{'ret'}) ) );
            allData(si).retDatFn = taskDir(ri).name;
            
            try
                allData(si).elecNotes = readtable([curDir(jj).folder '\' curDir(jj).name '\' curDir(jj).name '_Elec_Notes.xlsx']);
            catch
                allData(si).elecNotes = nan; 
            end

            allData(si).type = masterSheet.type(masteri); 
            allData(si).age = masterSheet.age(masteri); 
            allData(si).sex = masterSheet.sex(masteri); 
            allData(si).datNote = masterSheet.data(masteri); 
            allData(si).expNote = masterSheet.experimenter(masteri); 
            
            si = si+1;

        end



    end




end

end