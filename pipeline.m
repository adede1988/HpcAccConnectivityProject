

clear

%% initialize 

task = 'MemDev';
datFolder = 'R:\MSS\Johnson_Lab\DATA\';
masterSheet = readtable('R:\MSS\Johnson_Lab\dtf8829\memDevDat.csv');

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
            allData(si).retDatFn = taskDir(ei).name;
            
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


% MNI TEMPLATE BRAIN: 
% https://nist.mni.mcgill.ca/mni-average-brain-305-mri/
%points in a volume: 
% https://www.mathworks.com/matlabcentral/fileexchange/37856-inpolyhedron-are-points-inside-a-triangulated-volume

%% find all of the subject directories that have memDev data


for ii = 1:length(allData)
ii
dat = load([allData(ii).dataDir '\' allData(ii).encDatFn]).data;

if ~isfield(dat.elec, 'elecpos')
    dat.elec.elecpos = dat.elec.chanpos; 
end


[allData(ii).labels, allData(ii).labErrors] = getLabs(dat.elec, allData(ii).elecNotes);


end