%{
Calculate the average pathlength, reflectance, and timegate of the simulation

Clancy Lin
Last update: 2022/12/24
%}

clc;clear;close all;

%% param
input_folder='sim_2E10_n1457_diffNA_16'; % the simulation folder
output_dir='cal_PL';
sim_subDir='sim_'; % the prefix of the subDir
subject_name_arr={'CS'}; % the name of subjects
num_SDS=1; % number of SDS
num_layer=6; % number of layers
to_output_layers=1:5;
mua_change_rate=0.001; % the change rate of mua to calculate average PL
num_gate=25; % number of timegates
sim_index_set=load('thisPC_sim_wl_index.txt');

%% init
if exist(fullfile(input_folder,output_dir),'dir')==0
    mkdir(fullfile(input_folder,output_dir));
end

for sbj_i=1:length(subject_name_arr)
    subject_folder=fullfile(input_folder,subject_name_arr{sbj_i});

    for subSim_index=sim_index_set(1):sim_index_set(2)
        subSim_dir=fullfile(subject_folder,[sim_subDir num2str(subSim_index)]);

        fprintf('Processing %s sim %d\n',subject_name_arr{sbj_i},subSim_index);

        %% load mua file
        mua_arr=[0.6 0.55 0.4 0.6 0.3 0];

        %% main calculation
        gate_arr=cell(1,25);
        
        % for .mat format
        filename=fullfile(subSim_dir,'cfg_1.mat');
        load(filename)
        filename=fullfile(subSim_dir,'PL_1.mat');
        load(filename)
        detp.ppath = 10*SDS_detpt_arr{num_SDS};
        photon_weight = each_photon_weight_arr(num_SDS);
        tof=mcxdettime(detp,cfg.prop);
        [tempcounts, idx]=histc(tof,0:cfg.tstep:cfg.tend);
        tempcounts = tempcounts';

        % pathlength array divide into 25 groups by time gate
        last_percent=0;
        count=0;
        for k=1:length(detp.ppath)
            gate = idx(k);
            if gate==0
                gate = gate+1;
            end
            gate_arr{1,gate} = [gate_arr{1,gate}; detp.ppath(k,:)];
            
            percent = fix(k/length(detp.ppath)*100);
            if last_percent~=percent
                
                fprintf(1, repmat('\b',1,count));
                count=fprintf(1,'Divide progress : %d %%',percent);
            end
            last_percent=percent;
        end
        fprintf(1, repmat('\b',1,count));

        for gate_index=1:num_gate
            fprintf('%d,',gate_index);
            reflectance_arr=[];
            Pathlength_arr=[];
            reflectance_changed=[];

            reflectance_arr=1/photon_weight*sum(exp(sum(-1*double(gate_arr{gate_index}(:,:)).*mua_arr,2)));
            
            for l=1:num_layer
                changed_mua=mua_arr;
                changed_mua(l)=changed_mua(l)*(1+mua_change_rate);
                reflectance_changed(1,l)=1/photon_weight*sum(exp(sum(-1*double(gate_arr{gate_index}(:,:)).*changed_mua,2)));
            end
            Pathlength_arr(1:num_layer)=log(reflectance_arr./reflectance_changed)./((mua_arr*mua_change_rate));
            
            temp_PL_arr=[];
            temp_PL_arr(end+1:end+length(to_output_layers))=Pathlength_arr(to_output_layers);
            Pathlength_arr=temp_PL_arr;

            %% save
            output_folder=fullfile(output_dir,[sim_subDir num2str(subSim_index)]);

            if exist(fullfile(input_folder,output_folder),'dir')==0
                mkdir(fullfile(input_folder,output_folder));
            end

            output_folder=fullfile(output_dir,[sim_subDir num2str(subSim_index)]);
            save(fullfile(input_folder,output_folder,['gate_' num2str(gate_index) '_average_pathlength.txt']),'Pathlength_arr','-ascii','-tabs');

        end
        fprintf('\n\n');
        subSim_index=subSim_index+1;

    end
end

disp('Done!');