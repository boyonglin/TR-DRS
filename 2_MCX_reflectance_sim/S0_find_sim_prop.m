%{
Find simulation layer prop

Clancy Lin
Last update: 2022/12/25
%}

clc;clear;close all;

%% param
input_folder='sim_2E10_n1457_diffNA_16'; % the simulation folder
subject_name_arr={'CS'}; % the name of subjects
sim_index_set=load('thisPC_sim_wl_index.txt');

subject_folder=fullfile(input_folder,subject_name_arr{1});
load(fullfile(subject_folder,'mus_table.txt'))
mus_arr_temp=[];
mus_arr_temp(:,1)=[760 806 850];

for subSim_index=sim_index_set(1):sim_index_set(2)
    if mus_table(subSim_index,1:3)==[168 152 26]
        mus_arr_temp(1,end+1)=subSim_index;
    elseif mus_table(subSim_index,1:3)==[158 144 23]
        mus_arr_temp(2,end+1)=subSim_index;
    elseif mus_table(subSim_index,1:3)==[151 139 21]
        mus_arr_temp(3,end+1)=subSim_index;
    end
end

t = mus_arr_temp~=0;
n = sum(t,2);
m = max(n);
mus_arr = nan(size(mus_arr_temp,1),m);
for k = 1:size(mus_arr_temp,1)
    mus_arr(k,m-n(k)+1:m) = mus_arr_temp(k,t(k,:));
end

save('find_sim_prop.txt','mus_arr','-ascii','-tabs')

to_save=[mus_arr(:,1)'];
save('find_wl.txt','to_save','-ascii','-tabs')

disp('Done!')