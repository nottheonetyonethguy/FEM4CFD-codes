function writeoutputvtk(ndim,nelem,nnode,npelem,ndof,node_coords,elem_node_conn,soln,fileName)
%
% \brief Write the mesh and the solution in VTK legacy format
% \param  ndim         - dimension of the mesh
% \param  nElem        - number of elements in the mesh
% \param  nNode        - number of nodes in the mesh
% \param  npElem       - nodes per element (assumed to be the same for every element)
% \param  ndof         - number of dof per node
% \param  coords       - coordinates of the nodes
% \param  elem_node_conn - elemenet to node connectivity
% \param  elem_procid  - the id of the processor to which the element belongs to
% \param  soln         - solution at all the nodes

% \author Dr Chennakesava Kadapa
% \date 01 March 2025
%
% Open the file and write to it
fileID = fopen(fileName,'w');

% Directives
fprintf(fileID, "# vtk DataFile Version 4.0\n");
fprintf(fileID, "FEM example\n");

% ASCII or Binary
fprintf(fileID, "ASCII\n");

% Type of dataset : Structured/Unstructured/Rectilinear/Polydata...
fprintf(fileID, "DATASET UNSTRUCTURED_GRID\n" );

% Coordinates of the points (nodes)
fprintf(fileID, '%s \t %d \t %s \n', "POINTS ", nnode, " float");

if(ndim == 2)
  for ii=1:nnode
    fprintf(fileID, '%12.8f \t %12.8f \t %12.8f \n', node_coords(ii,1), node_coords(ii,2), 0.0);
  end
else
  for ii=1:nnode
    fprintf(fileID, '%12.8f \t %12.8f \t %12.8f \n', node_coords(ii,1), node_coords(ii,2), coords(ii,3));
  end
end

% Element<->Nodes connectivity
% In VTK terminology, Cell<->Points
% <number of nodes> <n1> <n2> <n3> ...
% Starting index is 0
ind = nelem*(npelem+1);
fprintf(fileID, '%s \t %d \t %d \n', "CELLS ", nelem, ind);

if(ndim == 2)
    if(npelem == 3) % Linear Triangular element
      n1 = 5;
      for ee=1:nelem
        fprintf(fileID, '%d \t %d \t %d \t %d\n', npelem, elem_node_conn(ee,1)-1, elem_node_conn(ee,2)-1, elem_node_conn(ee,3)-1);
      end
    elseif(npelem == 6) % Quadratic Triangular element
      n1 = 22;
      for ee=1:nelem
        fprintf(fileID, '%d \t %d \t %d \t %d \t %d \t %d \t %d\n', npelem, elem_node_conn(ee,1)-1, elem_node_conn(ee,2)-1, elem_node_conn(ee,3)-1,elem_node_conn(ee,4)-1,elem_node_conn(ee,5)-1,elem_node_conn(ee,6)-1);
      end
    else % Quadrilateral element
      n1 = 9;
      for ee=1:nelem
        fprintf(fileID, '%d \t %d \t %d \t %d \t %d\n', npelem, elem_node_conn(ee,1)-1, elem_node_conn(ee,2)-1, elem_node_conn(ee,3)-1,elem_node_conn(ee,4)-1);
      end
    end
else % ndim=3
    if(npelem == 4) % Tetrahedral element
      n1 = 10;
      for ee=1:nelem
        fprintf(fileID, '%d \t %d \t %d \t %d \t %d\n', npelem, elem_node_conn(ee,1)-1, elem_node_conn(ee,2)-1, elem_node_conn(ee,3)-1,elem_node_conn(ee,4)-1);
      end
    elseif(npelem == 8) % Hexahedral element
      n1 = 12;
      for ee=1:nelem
        fprintf(fileID, '%d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d\n', npelem, elem_node_conn(ee,1)-1, elem_node_conn(ee,2)-1, elem_node_conn(ee,3)-1,elem_node_conn(ee,4)-1,elem_node_conn(ee,5)-1,elem_node_conn(ee,6)-1,elem_node_conn(ee,7)-1,elem_node_conn(ee,8)-1);
      end
    end
end

% Cell type, as per VTK
fprintf(fileID, '%s \t %d \n', "CELL_TYPES", nelem);
for ee=1:nelem
    fprintf(fileID, '%d\n', n1);
end

% Point data
fprintf(fileID, '%s \t %d \n', "POINT_DATA", nnode);
if(ndof == 1)
    fprintf(fileID, '%s\n', "SCALARS solution float 1");
    fprintf(fileID, '%s\n', "LOOKUP_TABLE default");
    for ii=1:nnode
      fprintf(fileID, '%12.8f\n', soln(ii));
    end
else
    fprintf(fileID, '%s\n', "VECTORS solution float");
    if(ndof == 2)
      for ii=1:nnode
        ind = (ii-1)*ndof;
        fprintf(fileID, '%12.8f \t %12.8f \t %12.8f\n', soln(ind+1),soln(ind+2),0.0);
      end
    else
      for ii=1:nnode
        ind = (ii-1)*ndof;
        fprintf(fileID, '%12.8f \t %12.8f \t %12.8f\n', soln(ind+1),soln(ind+2),soln(ind+3));
      end
    end
end

% close the file
fclose(fileID);


