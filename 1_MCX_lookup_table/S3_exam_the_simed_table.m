%{
Check the simulated result is completed or not

Benjamin Kao
Last update: 2020/06/10
%}

clc;clear;close all;

%% param
input_folder='CS'; % the folder to exam
sim_index=1:18; % how many simulation are there to check
num_SDS=1; % how many SDS are there
num_photon=2E8; % the number of photon

use_par_mode=1; % =1 to use parallel (multi-core) at the same time

%% main
% load the one config file
one_conf=load(fullfile(input_folder,'sim_1','cfg_1.mat'));
one_detpos=one_conf.cfg.detpos;

if use_par_mode==0
    f=waitbar(0);
    for i=sim_index
        waitbar(((i-min(sim_index)+1)/length(sim_index)),f,['Checking sim ' num2str(i)]);
        try
            temp_sim_sum=jsondecode(fileread(fullfile(input_folder,['sim_' num2str(i)],'sim_summary_1.json')));
        catch
            assert(false,['ERROR: sim ' num2str(i) ' is not completed! sim summary not exist!']);
        end
        try
            temp_PL=load(fullfile(input_folder,['sim_' num2str(i)],'PL_1.mat'));
        catch
            assert(false,['ERROR: sim ' num2str(i) ' is not completed! PL.mat not exist!']);
        end
        try
            temp_conf=load(fullfile(input_folder,['sim_' num2str(i)],'cfg_1.mat'));
        catch
            assert(false,['ERROR: sim ' num2str(i) ' is not completed! cfg.mat not exist!']);
        end
        
        assert(sum(sum(temp_conf.cfg.detpos~=one_detpos))==0,['ERROR: sim ' num2str(i) '''s detpos is not consistent with other!']);
        
        assert(temp_sim_sum.num_photon==num_photon,['ERROR: the photon number of sim ' num2str(i) ' is wrong!']);
        for s=1:num_SDS
            sum_npt=temp_sim_sum.SDS_detected_number(s);
            pl_npt=size(temp_PL.SDS_detpt_arr{s},1);
            assert(sum_npt==pl_npt,['EEEOE: the detected photon number in sim ' num2str(i) ' SDS ' num2str(s) ' is inconsistent!']);
            if pl_npt==0
                fprintf('Warning: the detected photon number in sim %d SDS %d is 0!\n',i,s);
            end
        end
    end
    delete(f);
elseif use_par_mode==1
    fprintf('Checking sim :\n');
    parfor i=sim_index
        fprintf('Checking %d\n',i);
        try
            temp_sim_sum=jsondecode(fileread(fullfile(input_folder,['sim_' num2str(i)],'sim_summary_1.json')));
        catch
            assert(false,['ERROR: sim ' num2str(i) ' is not completed! sim summary not exist!']);
        end
        try
            temp_PL=load(fullfile(input_folder,['sim_' num2str(i)],'PL_1.mat'));
        catch
            assert(false,['ERROR: sim ' num2str(i) ' is not completed! PL.mat not exist!']);
        end
        try
            temp_conf=load(fullfile(input_folder,['sim_' num2str(i)],'cfg_1.mat'));
        catch
            assert(false,['ERROR: sim ' num2str(i) ' is not completed! cfg.mat not exist!']);
        end
        
        assert(sum(sum(temp_conf.cfg.detpos~=one_detpos))==0,['ERROR: sim ' num2str(i) '''s detpos is not consistent with other!']);
        
        assert(temp_sim_sum.num_photon==num_photon,['ERROR: the photon number of sim ' num2str(i) ' is wrong!']);
        for s=1:num_SDS
%             sum_npt=temp_sim_sum.SDS_detected_number(s);
            pl_npt=size(temp_PL.SDS_detpt_arr{s},1);
%             assert(sum_npt==pl_npt,['EEEOE: the detected photon number in sim ' num2str(i) ' SDS ' num2str(s) ' is inconsistent!']);
            if pl_npt==0
                fprintf('Warning: the detected photon number in sim %d SDS %d is 0!\n',i,s);
            end
        end
    end
end
disp('Done! No error found.');