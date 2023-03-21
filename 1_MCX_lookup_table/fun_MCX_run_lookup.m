%{
Build a lookup table by change mus
Using MCX with 5-layer MRI model
Change the param of scalp, skull, GM, mus of WM = 3* GM

Benjamin Kao
Last update: 2020/08/05
Version: 4.41
%}

function fun_MCX_run_lookup(target_name)

clc;clearvars -except target_name;close all;
%% param

sim_index_set=load('thisPC_sim_wl_index.txt');

%% init
global sim_version;
sim_version=4.41;
setting_version=load(fullfile(target_name,'sim_version.txt'));
assert(setting_version==sim_version,'ERROR: setting and main version not consistent!');

mus_table=load(fullfile(target_name,'mus_table.txt'));

load(fullfile(target_name,'sim_set.mat'));

model_folder='models';
MRI_voxel_file=fullfile(model_folder,['headModel' target_name '_EEG.mat']);
load(MRI_voxel_file);

if exist('model_version','var')==0
    model_version=1;
    do_sinus=0;
end

p_pos=load(fullfile(target_name,[target_name '_probe_pos.txt']));
p_dir=load(fullfile(target_name,[target_name '_probe_dir.txt']));
if model_version>=2
    p_pos=p_pos(:,[2 1 3]); % swap the x and y
    p_dir=p_dir(:,[2 1 3]);
end

%% check
assert(sim_index_set(1)<=sim_index_set(2),'Error: sim_index_set');
assert(sim_index_set(2)<=size(mus_table,1),'Error: sim_index_set');

%% check the previous sim result
simed_index_arr=zeros(size(mus_table,1),1);
for i=1:size(mus_table,1)
    if exist(fullfile(target_name,['sim_' num2str(i)],'sim_summary_1.json'),'file') || exist(fullfile(target_name,['sim_' num2str(i)],'old_sim_summary_1.json'),'file')
        simed_index_arr(i)=1;
    end
end

%% simulate
orig_sim_start_index=sim_index_set(1);
now_sim_index=sim_index_set(1);
while now_sim_index<=sim_index_set(2)
    if sim_index_set(1)~=orig_sim_start_index % if the sim start have been set to different value, than start sim from it
        now_sim_index=sim_index_set(1);
        orig_sim_start_index=sim_index_set(1);
    end
    if sim_index_set(2)<now_sim_index % the final index had been change to a smaller value
        break;
    end
    
    if simed_index_arr(now_sim_index)==0
        fprintf('sim %d th, %d / %d\n',now_sim_index,now_sim_index-sim_index_set(1)+1,sim_index_set(2)-sim_index_set(1)+1);
        output_folder=fullfile(target_name,['sim_' num2str(now_sim_index)]);
        mkdir(output_folder);
        
        temp_mus_arr=mus_table(now_sim_index,:);
        temp_param=zeros(1,10);
        temp_param(1,2:2:10)=temp_mus_arr;
        temp_param(1,1:2:10)=[0.6 0.55 0.4 0.6 0.3]; % 1/cm, for 700~1064nm, for check CV
%        temp_param(1,1:2:10)=[0.2 0.15 0.05 0.5 0.5]; % 1/cm, for 700~1064nm, for check CV
        if do_sinus==1
            temp_param=[temp_param 0 0.000001];
        end
        save(fullfile(output_folder,'mu.txt'),'temp_param','-ascii','-tabs');
        
        fun_MCX_sim_dist2axis(temp_param,vol,voxel_size,p_pos,p_dir,sim_set,output_folder,0,1,0,1,0,10000);
        
        simed_index_arr(now_sim_index)=1;
    end
    now_sim_index=now_sim_index+1;
    
    sim_index_set=load('thisPC_sim_wl_index.txt'); % load the setting to prevent edit
end

disp('Done!');

end
