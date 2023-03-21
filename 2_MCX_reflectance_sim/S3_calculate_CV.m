%{
Calculation of Coefficient of Variation (CV) at different total photon counts

Ting-Yi Guo
Last update: 2022/12/24
%}

dir = 'CS';
count = load(fullfile(dir,'Time_domain','final_counts.txt'));

scalp = [81 297 513 729 945 1161 1377 1593 1809 2025 2241 2457 2673];
skull = load(fullfile('skull_index.txt'));
fluid = load(fullfile('fluid_index.txt'));
gray = 943:948;

% for i = 1:6
%     figure;
%     plot(0:cfg.tstep:cfg.tend,new_count(i,:))
% end

% scalp
new_count=[];
for i=scalp
    new_count(end+1,:) = count(i,:);
end

scalp_reflentance_arr=[];

scalp_reflentance_arr(1,:)=mean(new_count);
scalp_reflentance_arr(2,:)=std(new_count);
scalp_reflentance_arr(3,:)=scalp_reflentance_arr(2,:)./scalp_reflentance_arr(1,:);

% skull
new_count=[];
for i=skull
    new_count(end+1,:) = count(i,:);
end

skull_reflentance_arr=[];

skull_reflentance_arr(1,:)=mean(new_count);
skull_reflentance_arr(2,:)=std(new_count);
skull_reflentance_arr(3,:)=skull_reflentance_arr(2,:)./skull_reflentance_arr(1,:);

% cerebrospinal fluid
new_count=[];
for i=fluid
    new_count(end+1,:) = count(i,:);
end

fluid_reflentance_arr=[];

fluid_reflentance_arr(1,:)=mean(new_count);
fluid_reflentance_arr(2,:)=std(new_count);
fluid_reflentance_arr(3,:)=fluid_reflentance_arr(2,:)./fluid_reflentance_arr(1,:);

% grey matter
new_count=[];
new_count = count(gray,:);

gray_reflentance_arr=[];

gray_reflentance_arr(1,:)=mean(new_count);
gray_reflentance_arr(2,:)=std(new_count);
gray_reflentance_arr(3,:)=gray_reflentance_arr(2,:)./gray_reflentance_arr(1,:);

figure;
plot(1:1:26,100*scalp_reflentance_arr(3,:))
hold on;
plot(1:1:26,100*skull_reflentance_arr(3,:))
hold on;
plot(1:1:26,100*fluid_reflentance_arr(3,:))
hold on;
plot(1:1:26,100*gray_reflentance_arr(3,:))
hold on;
title('The CV value of reflectance in different gates');
xlabel('time gate');
ylabel('CV(%)');
xlim([1 25])
legend('scalp','skull','cerebrospinal fluid','grey matter')

print(fullfile('reflectance.png'),'-dpng','-r200');
