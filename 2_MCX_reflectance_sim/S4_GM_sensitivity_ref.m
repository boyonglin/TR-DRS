%{
Calculate the sensitivity for each SDS by averaging the calculated sensitivities for different variations

Ting-Yi Guo
Last update: 2022/12/25
%}

clear all;

dir = 'CS';
% load(fullfile('tissue_index.mat'));
tissue_index = [1 2 4 5];

sds=1;
to_output_layer=1:6;
base_reflectance_arr = zeros(1,25);
reflectance_arr = zeros(1,25);
sensitivity_arr = [];

% delta_OP = [0.33333 0.23077 0.47368 0.35294]*100; % need to be modified
% delta_OP = [20 8 8 20];
delta_OP = [-20 -8 8 20];

% Calculate base reflectance
filename = fullfile(dir,'sim_3','cfg_1.mat');
load(filename)
filename = fullfile(dir,'sim_3','PL_1.mat');
load(filename)
filename = fullfile(dir,'sim_3','mu.txt');
mu = load(filename);

detp.ppath = 10*SDS_detpt_arr{sds};
photon_weight = each_photon_weight_arr(sds);
tof=mcxdettime(detp,cfg.prop);
[tempcounts, idx]=histc(tof,0:cfg.tstep:cfg.tend);
tempcounts = tempcounts';
detp.ppath = SDS_detpt_arr{sds};

for gate = 1:25
    gate_arr = [];
    index = find(idx==gate);
    gate_arr = detp.ppath(index,:);
    base_reflectance_arr(1,gate)=1/photon_weight*sum(exp(-double(sum(gate_arr.*mu((2*to_output_layer)-1),2))));
end


% Calculate other reflectance
for tissue = 1 % 1:4
    reflectance_arr = zeros(2,25);
    for num = 1:4
        filename = fullfile(dir,['sim_' num2str(tissue_index(tissue,num))],'cfg_1.mat');
        load(filename)
        filename = fullfile(dir,['sim_' num2str(tissue_index(tissue,num))],'PL_1.mat');
        load(filename)
        filename = fullfile(dir,['sim_' num2str(tissue_index(tissue,num))],'mu.txt');
        mu = load(filename);
        
        detp.ppath = 10*SDS_detpt_arr{sds};
        photon_weight = each_photon_weight_arr(sds);
        tof=mcxdettime(detp,cfg.prop);
        [tempcounts, idx]=histc(tof,0:cfg.tstep:cfg.tend);
        tempcounts = tempcounts';
%         temp_arr(run,:) = tempcounts/photon_weight;
        detp.ppath = SDS_detpt_arr{sds};
        
        for gate = 1:25
            gate_arr = [];
            index = find(idx==gate);
            gate_arr = detp.ppath(index,:);
            reflectance_arr(num,gate)=1/photon_weight*sum(exp(-double(sum(gate_arr.*mu((2*to_output_layer)-1),2))));
        end
    end
    
    sensitivity_arr(tissue,:) = sum(((reflectance_arr./base_reflectance_arr)-1)/(delta_OP(num)))/2;
%     sensitivity_arr(tissue,:) = sum(abs((reflectance_arr./base_reflectance_arr)-1)/(delta_OP(num)))/2;
    
    figure('Units','pixels','position',[0 0 1920 1080]);

    for gate = 1:25
        subplot(5,5,gate);
        x = [-20 -8 0 8 20];
        y = [reflectance_arr(1:2,gate)' base_reflectance_arr(gate) reflectance_arr(3:4,gate)'];
        plot(x,y,'-o')
        title(['Gate' num2str(gate)]);
        xlabel('us change(%)');
        ylabel('reflectance');
    end
end
print(fullfile('different_gate.png'),'-dpng','-r200');
plot_gate = 25;

figure;
plot(1:1:plot_gate,sensitivity_arr(1,1:plot_gate))
% hold on;
% plot(1:1:plot_gate,sensitivity_arr(2,1:plot_gate))
% hold on;
% plot(1:1:plot_gate,sensitivity_arr(3,1:plot_gate))
% hold on;
% plot(1:1:plot_gate,sensitivity_arr(4,1:plot_gate))
% hold on;
title('Grey matter sensitivity of each time gate');
xlabel('time gate');
ylabel('sensitivity');

% legend('scalp','skull','cerebrospinal fluid','grey matter')

print(fullfile('reflectance.png'),'-dpng','-r200');
