
% Make a table with rxn'IDs and subsystems 
% The info are extrapolated from Recon3DModel_301.mat original model

m = readCbModel('models/Recon3D_301/Recon3DModel_301.mat');

ss ={};

for i = 1:length(m.rxns)  
    
    temp = [m.subSystems{i}];
    ss{i} = join(temp, ' * ');
    
end
ss = transpose(ss);

rxns_subSystems = [m.rxns,ss];
rxns_subSystems = cell2table(rxns_subSystems);

writetable(rxns_subSystems, 'rxns_subSystems.txt','Delimiter', '\t')

