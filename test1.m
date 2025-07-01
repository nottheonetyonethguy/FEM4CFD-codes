model = createpde("thermal","steadystate");

T1 = [2,4,0,0.866,0.866,0,0,0.5,1.5,2]';
g = decsg(T1);

geometryFromEdges(model,g);
mesh = generateMesh(model,"GeometricOrder","linear","Hmax",0.1);
% pdegplot(model,'Edgelabels','on')
pdeplot(model,"NodeLabels","on")

%% BC setting

thermalBC(model,"Edge",2,"ConvectionCoefficient",30,"AmbientTemperature",30);
thermalBC(model,"Edge",4,"Temperature",0);
thermalBC(model,"Edge",[1,3],"HeatFlux",0);
internalHeatSource(model,100)
thermalProperties(model,"ThermalConductivity",50)

temperatures = solve(model).Temperature;
conductivity = 50;
heatSource = 100;
convectionCoefficient = 30;
ambientTemperature = 30;
%%

NL = mesh.Nodes';
EL = mesh.Elements';
EL_size = length(EL);
NL_size = length(NL);
EL_list = zeros(EL_size,1);
NL_list = zeros(NL_size,1);
for i = 1:EL_size
    EL_list(i) = i;
    EL_to_print(i,1:6) = [EL_list(i) EL(i,1:3) conductivity heatSource];
end
for i = 1:NL_size
    NL_list(i) = i;
end


%% writing element data to the dat file
fileID = fopen('data.dat','a+');
fprintf(fileID,'TITLE = %d-elem-problem\n',EL_size);
fprintf(fileID,'\nELEMENTS = %d\n',EL_size);
writematrix(EL_to_print,'data.dat','Delimiter',' ',WriteMode='append');

%% writing node data to the dat file
fprintf(fileID,'\nNODE_COORDINATES = %d', NL_size);

writematrix([NL_list NL],'data.dat','Delimiter',' ',WriteMode='append');
%% known edge temperatures

NL_x = NL(:,1);
c = 1;
for i = 1:NL_size
    if NL_x(i) == 0
        known_temp_node(c,1) = i;
        known_temp_node(c,2) = temperatures(i);
        c = c + 1;
    end
end

fprintf(fileID,'\nNODES_WITH_PRESCRIBED_TEMPERATURE = %d',c-1);
writematrix(known_temp_node,'data.dat','Delimiter',' ',WriteMode='append');

%% convection edges

c_2 = 1;
for i = 1:NL_size
    if NL_x(i) == 0.866
        known_convec_node(c_2,1) = i;
        known_convec_node(c_2,2:3) = NL(i,:);
        c_2 = c_2 + 1;
    end
end
sorted_convec_nodes = sortrows(known_convec_node,3);

for i = 1:length(sorted_convec_nodes)
    if i < length(sorted_convec_nodes)
        edges_to_print(i,1:4) = [sorted_convec_nodes(i,1) sorted_convec_nodes(i+1,1) convectionCoefficient ambientTemperature];
    end
end


fprintf(fileID,'\nEDGES_WITH_PRESCRIBED_CONVECTION = %d',c_2-2);
writematrix(edges_to_print,'data.dat','Delimiter',' ',WriteMode='append');

%%
fclose(fileID)
