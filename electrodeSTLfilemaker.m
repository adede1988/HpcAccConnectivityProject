

addpath('C:\Users\dtf8829\Documents\GitHub\myFrequentUse')


load('R:\MSS\Johnson_Lab\DATA\LCH\LC04\Recon\Recon_Jan_2023\FT_Pipeline\Electrodes\LC04_elec_acpc_f.mat')


X = []; 
Y = []; 
Z = []; 

for ii = 1:size(elec_acpc_f.elecpos,1)

    [x,y,z] = sphere; 
    x = x+elec_acpc_f.elecpos(ii,1); 
    y = y-elec_acpc_f.elecpos(ii,2); 
    z = z+elec_acpc_f.elecpos(ii,3); 
    
    X = [X, x]; 
    Y = [Y, y]; 
    Z = [Z, z]; 
end


surf2stl_elec('C:\Users\dtf8829\Documents\GitHub\HpcAccConnectivityProject\nativeElecs.stl', X, Y, Z);



load('R:\MSS\Johnson_Lab\DATA\LCH\LC04\Recon\Recon_Jan_2023\FT_Pipeline\Electrodes\LC04_elec_mni_frv.mat')


X = []; 
Y = []; 
Z = []; 

for ii = 1:size(elec_mni_frv.elecpos,1)

    [x,y,z] = sphere; 
    x = x+elec_mni_frv.elecpos(ii,1); 
    y = y-elec_mni_frv.elecpos(ii,2); 
    z = z+elec_mni_frv.elecpos(ii,3); 
    
    X = [X, x]; 
    Y = [Y, y]; 
    Z = [Z, z]; 
end


surf2stl_elec('C:\Users\dtf8829\Documents\GitHub\HpcAccConnectivityProject\mniElecs.stl', X, Y, Z);