close all;
clear all;
addpath('./1. BasicModules/kd_tree');
addpath('./1. BasicModules');
% addpath('./2. Non_rigid_registration');
%% Example
if ~exist('VS','var') 
    [VS, FS, NS] = read_obj_file('E:/data/model661.obj');
    [VS2, FS2, NS2] = read_obj_file('E:/data/model771_4.obj');
    [VT, FT, NT] = read_obj_file('E:/data/nt1.obj');    
end
corres = cell(size(FT,1),1);
for i=1:size(FT,1)
    temp = i;
    corres{i} =temp;
end
[ x, nx ] = deformation_transfer(VS, FS, VT, FT, VS2, FS2, corres);
write_obj_file(x, FT, nx, 'E:/data/4.obj');








fprintf('End of demo..\n');


clear VS VT S_factor T_factor FS FT NS NT maker;