%{
Calculate the grey matter sensitivity of pathlength simulation

Clancy Lin
Last update: 2022/12/25
%}

clc;clear;close all;

input_folder='sim_2E10_n1457_diffNA_16';
input_subDir='cal_PL';
output_dir='gate_PL';
sim_subDir='sim_'; % the prefix of the subDir
to_plot_SDS_layer_index=4; % the index of the SDS and the layer to plot
num_gate=25; % number of timegates
sim_index_set=load('thisPC_sim_wl_index.txt');
load('find_sim_prop.txt')
find_wl=load('find_wl.txt');

GM_sen_arr=[];
sbj_PL_arr=[];
timegate_x=1:1:25;

for gate_index=1:num_gate
    for subSim_index=sim_index_set(1):sim_index_set(2)
        Pathlength_arr=load(fullfile(input_folder,input_subDir,['sim_' int2str(subSim_index)],['gate_' num2str(gate_index) '_average_pathlength.txt']));
        total_pathlength_arr(:,:,subSim_index)=Pathlength_arr;
    end
    mean_pathlength_arr=mean(total_pathlength_arr,3);

    if exist(fullfile(input_folder,output_dir),'dir')==0
        mkdir(fullfile(input_folder,output_dir));
    end

    save(fullfile(input_folder,output_dir,['gate_' num2str(gate_index) '_mean_pathlength.txt']),'mean_pathlength_arr','-ascii','-tabs');

    sbj_PL_arr(1,gate_index)=mean_pathlength_arr(:,to_plot_SDS_layer_index);
    sbj_PL_arr(2,gate_index)=sum(mean_pathlength_arr(:));
end

for gate_index=1:num_gate
    total_PL_arr=sum(sbj_PL_arr(2,:));
    GM_sen=sbj_PL_arr(1,gate_index)/total_PL_arr*100;
    GM_sen_arr(:,gate_index)=GM_sen;
end

figure
plot(timegate_x,GM_sen_arr)
title('Grey matter sensitivity of all time gate')
xlim([1 25])
xticks([1 5 10 15 20 25])
xlabel('time gate')
ylabel('sensitivity(%)')

print(fullfile(input_folder,'GM_sen_all.png'),'-dpng','-r200');

for gate_index=1:num_gate
    total_PL_arr=sum(sbj_PL_arr(2,:));
    GM_sen = sbj_PL_arr(1,gate_index)/sbj_PL_arr(2,gate_index)*100;
    GM_sen_arr(:,gate_index)=GM_sen;
end

figure
plot(timegate_x,GM_sen_arr)
title('Grey matter sensitivity of each time gate')

xlim([1 25])
xticks([1 5 10 15 20 25])
xlabel('time gate')
ylabel('sensitivity(%)')
print(fullfile(input_folder,'GM_sen_each.png'),'-dpng','-r200');

%% Find in different grey matter mus

f=figure;
t=tiledlayout(3,1);

for wl_index=1:length(find_wl)
    nexttile
    layer_mus_prop= find_sim_prop(wl_index,2:end);
    GM_sen_plot(layer_mus_prop,num_gate,input_folder,input_subDir,to_plot_SDS_layer_index,timegate_x)
    title(sprintf('%d nm',find_wl(wl_index)))
    legend('Location','northwest')
end

title(t,'Grey matter sensitivity in different GM layer mus')
xlim([1 25])
xticks([1 5 10 15 20 25])
xlabel(t,'time gate')
ylabel(t,'sensitivity^2(%)')
% line = legend();
% line.Layout.Tile = 'east';
f.Position = [250 250 600 600];

print(fullfile(input_folder,'GM_sen_mus.png'),'-dpng','-r200');

%% Average mus

f2=figure;
t2=tiledlayout(3,1);

for wl_index=1:length(find_wl)
    nexttile
    layer_mus_prop= find_sim_prop(wl_index,2:end);
    for gate_index=1:num_gate
        for subSim_index=layer_mus_prop(1):layer_mus_prop(length(layer_mus_prop))
            Pathlength_arr=load(fullfile(input_folder,input_subDir,['sim_' int2str(subSim_index)],['gate_' num2str(gate_index) '_average_pathlength.txt']));
            total_pathlength_arr(:,:,subSim_index)=Pathlength_arr;
        end
        mean_pathlength_arr=mean(total_pathlength_arr,3);
    
        if exist(fullfile(input_folder,output_dir),'dir')==0
            mkdir(fullfile(input_folder,output_dir));
        end
    
        sbj_PL_arr(1,gate_index)=mean_pathlength_arr(:,to_plot_SDS_layer_index);
        sbj_PL_arr(2,gate_index)=sum(mean_pathlength_arr(:));
    end
    
    for gate_index=1:num_gate
        total_PL_arr=sum(sbj_PL_arr(2,:));
%         GM_sen = sbj_PL_arr(1,gate_index)/sbj_PL_arr(2,gate_index)*100;
        GM_sen = sbj_PL_arr(1,gate_index)/total_PL_arr*100;
        GM_sen_arr(:,gate_index)=GM_sen;
    end
    plot(timegate_x,GM_sen_arr)
    title(sprintf('%d nm',find_wl(wl_index)))
end

title(t2,'Average grey matter sensitivity in different GM layer mus')
xlim([1 25])
xticks([1 5 10 15 20 25])
xlabel(t2,'time gate')
ylabel(t2,'sensitivity(%)')
f2.Position = [250 250 600 600];

print(fullfile(input_folder,'GM_sen_mus_avg.png'),'-dpng','-r200');


%% Calculate grey matter sensitivity and plot

function GM_sen_plot(layer_mus_prop,num_gate,input_folder,input_subDir,to_plot_SDS_layer_index,timegate_x)
    for subSim_index=layer_mus_prop(1):layer_mus_prop(length(layer_mus_prop))
        for gate_index=1:num_gate
            Pathlength_arr=load(fullfile(input_folder,input_subDir,['sim_' int2str(subSim_index)],['gate_' num2str(gate_index) '_average_pathlength.txt']));
            layer_PL_arr(1,gate_index)=Pathlength_arr(:,to_plot_SDS_layer_index);
            layer_PL_arr(2,gate_index)=sum(Pathlength_arr(:));
        end
        for gate_index=1:num_gate
            total_PL_arr=sum(layer_PL_arr(2,:));
            layer_GM_sen(:,gate_index)=power(layer_PL_arr(1,gate_index)/total_PL_arr*100,2);
        end
%         plot(timegate_x,layer_PL_arr(1,:),'DisplayName',sprintf('mus=%d',(subSim_index-layer_mus_prop(1))*60+25))
        plot(timegate_x,layer_GM_sen,'DisplayName',sprintf('mus=%d',(subSim_index-layer_mus_prop(1))*60+25))
        hold on
    end
    hold off
end