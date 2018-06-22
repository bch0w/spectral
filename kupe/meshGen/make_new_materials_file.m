%
% make_new_materials_file.m
%
% This script will reassign each 'material index' for each element in a
% mesh to be a negative number that corresponds to an external xyz
% tomographic file on a regular grid.
%
% The purpose of this script is to allow for multiple xyz files to be used
% -- or perhaps a combination of xyz files with constant-valued properties
% (e.g., a uniform value for the upper mantle).
%
% This should be run from the directory containing your geocubit-produced
% files. (Use a symbolic link for convenience.)
%
% Be sure to copy materials_file_tomo to materials_file at the end, since
% SPECFEM3D is expecting the file name materials_file.
%

% load nodes file and mesh file
% THESE FILE NAMES ARE DEFAULT FROM GEOCUBIT
% (prior to this, you need to do 'mv materials_file materials_file_orig')
[inodes,nx,ny,nz] ...
    = textread('nodes_coords_file','%f%f%f%f','headerlines',1);
[imesh,in1,in2,in3,in4,in5,in6,in7,in8] ...
    = textread('mesh_file','%f%f%f%f%f%f%f%f%f','headerlines',1);
[imat,matind] = textread('materials_file','%f%f');

disp('loaded all files');

% check (technically we can just use imesh, not imat)
sum(imesh-imat)

node_inds_for_each_element = [in1 in2 in3 in4 in5 in6 in7 in8];
n = length(imat);

% The values for zval ought to agree with the values used in
% run_get_tomodd.m for creating the multiple xyz files, each with different
% resolutions.

filename = './materials_file_tomo';
fid = fopen(filename,'w');
for ii=1:n
    % get the z values of all the nodes that define the element
    inode = node_inds_for_each_element(ii,:);
    zvals = nz(inode);
    zval = min(zvals);  % find the DEEPEST node on the element
    kk = imat(ii);

%     % homogeneous halfspace test
%     if zval >= -10*1e3
%         fprintf(fid,'%i %i\n',kk,-1);
%     elseif and( zval >= -50*1e3, zval < -10*1e3)
%         fprintf(fid,'%i %i\n',kk,-2);
%     else
%         fprintf(fid,'%i %i\n',kk,-3);
%     end
  
    % New Zealand (nz) -- Moho fixed at 33 km
    if matind(kk)==1            % mantle
        fprintf(fid,'%i %i\n',kk,-3);
    else
        if zval >= -8*1e3       % shallow
            fprintf(fid,'%i %i\n',kk,-1);
        else                    % crust
            fprintf(fid,'%i %i\n',kk,-2);
        end
    end  
    
%     % nenana_basin
%     if zval >= -8*1e3
%         fprintf(fid,'%i %i\n',kk,-1);
%     else
%         fprintf(fid,'%i %i\n',kk,-2);
%     end 

%     % cook_basin
%     if zval >= -8*1e3
%         fprintf(fid,'%i %i\n',kk,-1);
%     elseif and( zval >= -50*1e3, zval < -8*1e3)
%         fprintf(fid,'%i %i\n',kk,-2);
%     else
%         fprintf(fid,'%i %i\n',kk,-3);
%     end    
    
%     % new zealand
%     if zval >= -8*1e3
%         fprintf(fid,'%i %i\n',kk,-1);
%     elseif and( zval >= -50*1e3, zval < -8*1e3)
%         fprintf(fid,'%i %i\n',kk,-2);
%     else
%         fprintf(fid,'%i %i\n',kk,-3);
%     end
    
%     % alaska 1964gulf for 2018 earthquake (sak)
%     if zval >= -8*1e3
%         fprintf(fid,'%i %i\n',kk,-1);
%     elseif and( zval >= -50*1e3, zval < -8*1e3)
%         fprintf(fid,'%i %i\n',kk,-2);
%     else
%         fprintf(fid,'%i %i\n',kk,-3);
%     end
    
%     % alaska 1964gulf for 2018 earthquake (sak)
%     if zval >= -8*1e3
%         fprintf(fid,'%i %i\n',kk,-1);
%     elseif and( zval >= -100*1e3, zval < -8*1e3)
%         fprintf(fid,'%i %i\n',kk,-2);
%     else
%         fprintf(fid,'%i %i\n',kk,-3);
%     end

%     % alaska 1964 (sak)
%     if zval >= -10*1e3
%         fprintf(fid,'%i %i\n',kk,-1);
%     elseif and( zval >= -50*1e3, zval < -10*1e3)
%         fprintf(fid,'%i %i\n',kk,-2);
%     else
%         fprintf(fid,'%i %i\n',kk,-3);
%     end

%     % san joaquin
%     if zval >= -11*1e3
%         fprintf(fid,'%i %i\n',kk,-1);
%     else
%         fprintf(fid,'%i %i\n',kk,matind(ii));
%     end

%     % socal
%     if zval >= -15*1e3
%         fprintf(fid,'%i %i\n',kk,-1);
%     else
%         fprintf(fid,'%i %i\n',kk,-2);
%     end
    
%     % socal
%     if zval >= -15*1e3
%         fprintf(fid,'%i %i\n',kk,-1);
%     elseif and( zval >= -45*1e3, zval < -15*1e3)
%         fprintf(fid,'%i %i\n',kk,-2);
%     else
%         fprintf(fid,'%i %i\n',kk,-3);
%     end
    
end
fclose(fid);

%==========================================================================
