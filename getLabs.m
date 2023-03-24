function [labels, errorReport] = getLabs(varargin) 

switch nargin
    case 1
        elec = varargin{1}; 
        notes = [];
        hasNotes = false; 
    case 2
        elec = varargin{1};
        notes = varargin{2};
        hasNotes = true; 
end


%% check if there's anything in the meeting notes
try
  
    if ~iscell(notes.LocMeeting(1))
        hasNotes = false; 
    end
    
catch

    try 
        if  ~iscell(notes.Loc_meeting(1))

            hasNotes = false; 
        else

            notes = renamevars(notes, "Loc_meeting", "LocMeeting");
        end
    catch
        if isnan(notes)
            hasNotes = false; 
        end
    end
    
end

%% check for dumbly named electrodes

try
if strcmp("A'1", notes.Electrode(1))
    for li = 1:length(notes.Electrode)
        cur = split(notes.Electrode(li), "'");
        if length(cur)==2
        notes.Electrode(li) = {[cur{1} cur{2} '_left']};
        end
    end
elseif length(split(notes.Electrode(1), " " ) )>1
    for li = 1:length(notes.Electrode)
        cur = split(notes.Electrode(li), " ");
        if length(cur)==2
        notes.Electrode(li) = {[cur{1} cur{2}]};
        end
    end

elseif strcmp(notes.Properties.VariableNames(2), 'originalLabel')
    notes.Electrode = notes.originalLabel;
    for li = 1:length(elec.label)
        curLab = elec.label{li}; 
        [ch, ~] = split(curLab, '-');
        if length(ch{1}) == 2 && length(ch{2}) == 2
            elec.label{li} = [ch{1}(1) '0' ch{1}(2) '-' ch{2}(1) '0' ch{2}(2) ];
        elseif length(ch{1}) ==2
            elec.label{li} = [ch{1}(1) '0' ch{1}(2) '-' ch{2} ];
        elseif length(ch{2}) ==2
            elec.label{li} = [ch{1} '-' ch{2}(1) '0' ch{2}(2) ];
        end

    end
    for li = 1:length(notes.Electrode)
        if length(notes.Electrode{li}) == 2
            notes.Electrode{li} = [notes.Electrode{li}(1) '0' notes.Electrode{li}(2)];
        end
    end

end
catch
disp(['misnamed note'])
end

%% check for nan electrodes and eliminate

if sum(isnan(elec.elecpos(:,1)))>0
    deletei = isnan(elec.elecpos(:,1)); 
    elec.label(deletei) = [];
    elec.elecpos(deletei,:) = [];
    elec.chanpos(deletei,:) = []; 
end


%% hard code a bunch of codes: 

ignore = {'wm', 'ofb', 'oob', 'white matter', 'ventricle', 'csf', 'WM', 'out of brain', 'white mat', 'lateral ventracle', 'corpus collosum', 'white', 'lesion'}; 
sfs = {'sfs', 'SFS', 'superior frontal sulcus', 'SFC'};
acc = {'acc', 'aCC', 'ant cingulate', 'anterior cingulate', 'anterior cingulate cortex', 'corpus collosum/cing', 'cingulate'}; 
hip = {'hip', 'ca1', 'CA1', 'dg', 'DG', 'CA3', 'ca3', 'Hip', 'subiculum'}; 
sfg = {'sfg', 'SFG', 'sup frontal', 'superior frontal gyrus'};
mfg = {'mfg', 'med frontal', 'medial frontal gyrus', 'mid frontal gyrus', 'MFG', 'middle frontal gyrus', 'post med frontal/supp motor'}; 
ifg = {'ifg', 'inf frontal', 'IFG', 'inferior frontal gyrus', 'Inferior frontal gyrus', 'inferior lateral frontal', 'inferior temporal cortex'}; 
mfs = {'mfs', 'medial frontal sulcus', 'MFS', 'MF sulcus', 'middle frontal sulcus'}; 
cls = {'cls', 'claustrum'}; 
ins = {'ins', 'insula', 'Insula', 'parietal operculum'}; 
mpfc= {'MPFC', 'mpfc', 'mPFC', 'medial PFC', 'medial frontal cortex'}; 
ofc = {'ofc', 'OFC', 'orbital frontal cortex', 'OFG', 'orbital frontal'}; 
bg = {'striatum', 'str', 'Caudate', 'caudate', 'BG'};
ifs = {'IFS', 'inferior frontal sulcus', 'ifs', 'IFG sulcus', 'inferior/frontal sulcus'}; 
phg = {'ERC', 'EC', 'PHG', 'phg', 'lat EC', 'perirhinal', 'parahippocampal gyrus', 'perrhinal cortex', 'entorhinal cortex'}; 
mts = {'MTS', 'mts', 'mt sulcus', 'medial temporal sulcus', 'middle temporal sulcus', 'mid temp sulcus'};
stg = {'STG', 'stg', 'superior temporal gyrus', 'sup temporal gyrus'};
mtg = {'MTG', 'middle temp gyrus', 'middle temporal gyrus', 'mid temp gyrus' 'medial temporal gyrus', 'mtg', 'mt gyrus'};
fus = {'fusiform', 'Fusiform'};
itg = {'itg', 'ITG', 'inferior temporal gyrus', 'inf temporal gyrus'}; 
amy = {'amy', 'amygdala', 'Amygdala'}; 
ats = {'anterior temporal sulcus'};
tp =  {'anterior temporal', 'ant temporal pole', 'temporal pole'};
its = {'inferior temporal sulcus', 'inferior temp sulcus', 'It sulcus', 'ITS'};
ipl = {'IPL', 'ipl','angular gyrus', 'Inferior Pariatal Lobule', 'IPS', 'supramarginal'}; 
prcg = {'Precentral gyrus', 'precentral gyrus', 'precentral', 'M1', 'dorsal motor cortex'}; 
pocg = {'postcentral gyrus', 'Postcentral gyrus', 'S1', 's1'}; 
sts = {'STS', 'sts', 'inferior bank of sts', 'superior temporal sulcus'}; 
poc = {'Parietal occipital cortex', 'Parietal occipital sulcus'};
mpa = {'Medial Parietal occipital', 'medial parietal', 'precuneus'}; 
occ = {'anterior occipital', 'occipital cortex', 'V2', 'lingual', 'occ', 'medial occipital', 'occipital', 'calcarine'};
toc = {'temporal occipital'};
tpj = {'temporal parietal junction', 'TPJ', 'temporal/parietal cortex', 'temporal/parietal junction'}; 
mcc = {'mCC', 'mid/post cing', 'mcc', 'MCC', 'middle cingulate cortex', 'medial cingulate cortex'};
cs = {'superior central sulcus', 'central sulcus'}; 
spl ={'superior parietal lobule', 'SPL sulcus'};
par = {'parietal cortex', 'parietal'}; 
thal = {'thalamus', 'STN', 'Thalamus'};





%% now find labels
%bug for later: spit out some record of second channel conflicts where they
%exist 

errorReport = struct; 
ei = 1; 
labels = cell(size(elec.label,1),3); 
for li = 1:length(elec.label)
    curLab = elec.label{li}; 
    [ch, ~] = split(curLab, '-');
    
    if hasNotes
        
        try
            ch1i = find(cellfun(@(x) ~isempty(x), strfind(notes.Electrode, ch{1})),1);
            ch2i = find(cellfun(@(x) ~isempty(x), strfind(notes.Electrode, ch{2})),1);
        catch
            ch1i = find(cellfun(@(x) ~isempty(x), strfind(notes.Label, ch{1})),1);
            ch2i = find(cellfun(@(x) ~isempty(x), strfind(notes.Label, ch{2})),1);
        end

%         if strcmp('originalLabel', notes.Properties.VariableNames{2})
%             ch1i = find(cellfun(@(x) ~isempty(x), strfind(notes.originalLabel, ch{1})),1);
%             ch2i = find(cellfun(@(x) ~isempty(x), strfind(notes.originalLabel, ch{2})),1);
%         end



%         ch1 = notes.LocMeeting{ch1i};
%         ch2 = notes.LocMeeting{ch2i};

        if isempty(notes.LocMeeting{ch1i}) && isempty(notes.LocMeeting{ch2i})
            labels{li, 1} = {'NO NOTES'};
            labels{li, 2} = {'NO NOTES'};
        else
%         
%         %switch channel 2 into channel 1 position if it's got nothing,
%         %otherwise default to using channel 1 location meeting note
%         if sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(x, ch1), ignore, 'uniformoutput', false))) > 0
%             ch1 = ch2; 
%             if sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(x,ch1), ignore, 'uniformoutput', false))) > 0
%                 errorReport(ei).elec = li; 
%                 errorReport(ei).flag = 'bothBad';
%                 ei = ei+1; 
%                 labels{li} = 'ZZZ';
%             end
%         end

        chi = [ch1i, ch2i]; 
        for ii = 1:2
            ch = notes.LocMeeting{chi(ii)};

        if sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), sfs, 'uniformoutput', false))) > 0
                labels{li,ii} = 'sfs'; 
                %superior frontal sulcus
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), acc, 'uniformoutput', false))) > 0
                labels{li,ii} = 'acc'; 
                %anterior cingulate cortex
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), hip, 'uniformoutput', false))) > 0
                labels{li,ii} = 'hip'; 
                %anterior cingulate cortex  
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), sfg, 'uniformoutput', false))) > 0
                labels{li,ii} = 'sfg'; 
                %anterior cingulate cortex  
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), mfg, 'uniformoutput', false))) > 0
                labels{li,ii} = 'mfg'; 
                %medial frontal gyrus 
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), ifg, 'uniformoutput', false))) > 0
                labels{li,ii} = 'ifg'; 
                %inferior frontal gyrus  
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), mfs, 'uniformoutput', false))) > 0
                labels{li,ii} = 'mfs'; 
                %middle frontal sulcus  
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), cls, 'uniformoutput', false))) > 0
                labels{li,ii} = 'cls'; 
                %claustrum  
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), ins, 'uniformoutput', false))) > 0
                labels{li,ii} = 'ins'; 
                %insula 
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), mpfc, 'uniformoutput', false))) > 0
                labels{li,ii} = 'mpfc'; 
                %medial prefrontal cortex 
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), ofc, 'uniformoutput', false))) > 0
                labels{li,ii} = 'ofc'; 
                %orbital frontal cortex
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), bg, 'uniformoutput', false))) > 0
                labels{li,ii} = 'bg'; 
                %striatum
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), ifs, 'uniformoutput', false))) > 0
                labels{li,ii} = 'ifs'; 
                %inferior frontal sulcus
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), phg, 'uniformoutput', false))) > 0
                labels{li,ii} = 'phg'; 
                %parahippocampal gyrus
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), mts, 'uniformoutput', false))) > 0
                labels{li,ii} = 'mts'; 
                %medial temporal sulcus
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), stg, 'uniformoutput', false))) > 0
                labels{li,ii} = 'stg'; 
                %superior temporal gyrus
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), mtg, 'uniformoutput', false))) > 0
                labels{li,ii} = 'mtg'; 
                %middle temporal gyrus
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), fus, 'uniformoutput', false))) > 0
                labels{li,ii} = 'fus'; 
                %fusiform gyrus
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), itg, 'uniformoutput', false))) > 0
                labels{li,ii} = 'itg'; 
                %inferior temporal gyrus
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), amy, 'uniformoutput', false))) > 0
                labels{li,ii} = 'amy'; 
                %amygdala
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), ats, 'uniformoutput', false))) > 0
                labels{li,ii} = 'ats'; 
                %anterior temporal sulcus
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), tp, 'uniformoutput', false))) > 0
                labels{li,ii} = 'tp'; 
                %temporal pole
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), its, 'uniformoutput', false))) > 0
                labels{li,ii} = 'its'; 
                %inferior temporal sulcus
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), ipl, 'uniformoutput', false))) > 0
                labels{li,ii} = 'ipl'; 
                %inferior parietal lobule
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), prcg, 'uniformoutput', false))) > 0
                labels{li,ii} = 'prcg'; 
                %precentral gyrus
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), pocg, 'uniformoutput', false))) > 0
                labels{li,ii} = 'pocg'; 
                %postcentral gyrus
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), sts, 'uniformoutput', false))) > 0
                labels{li,ii} = 'sts'; 
                %superior temporal sulcus
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), poc, 'uniformoutput', false))) > 0
                labels{li,ii} = 'poc'; 
                %parietal occipital cortex
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), mpa, 'uniformoutput', false))) > 0
                labels{li,ii} = 'mpa'; 
                %medial parietal cortex
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), occ, 'uniformoutput', false))) > 0
                labels{li,ii} = 'occ'; 
                %occipital lobe
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), toc, 'uniformoutput', false))) > 0
                labels{li,ii} = 'toc'; 
                %temporal occipital 
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), tpj, 'uniformoutput', false))) > 0
                labels{li,ii} = 'tpj'; 
                %temporal parietal junction 
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), mcc, 'uniformoutput', false))) > 0
                labels{li,ii} = 'mcc'; 
                %medial cingulate cortex
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), cs, 'uniformoutput', false))) > 0
                labels{li,ii} = 'cs'; 
                %central sulcus
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), par, 'uniformoutput', false))) > 0
                labels{li,ii} = 'par'; 
                %parietal cortex
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), spl, 'uniformoutput', false))) > 0
                labels{li,ii} = 'spl'; 
                %superior parietal lobule
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), thal, 'uniformoutput', false))) > 0
                labels{li,ii} = 'thal'; 
                %thalamus/midbrain
        elseif sum(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(ch, x), ignore, 'uniformoutput', false))) > 0
                labels{li,ii} = 'ZZZ'; 
                %white/out of brain




        else
    
                errorReport(ei).elec = li; 
                errorReport(ei).flag = 'unknown'; 
                errorReport(ei).inName = ch; 
                ei = ei+1;



        end
        end
        end


    else %anatomy based! 

        labels{li, 1} = {'NO NOTES'};
        labels{li, 2} = {'NO NOTES'};


    end



   
   




end


    
    


 %make third position label: 
 for li = 1:length(elec.label)
    ch1 = labels{li,1}; 
    ch2 = labels{li,2};
    if strcmp('ZZZ', ch1) || strcmp('NO NOTES', ch1)
        labels{li,3} = ch2; 
    else
        labels{li,3} = ch1; 
    end


end













end