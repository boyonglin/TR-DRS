%{
Find the position and direction of the probe before simulation.

Benjamin Kao
Last update: 2020/11/28
%}

clc;clear;close all;

%% param
model_folder='models';
subject_name_arr={'CS'};

to_plot_figure=1; % plot the probe position and direction or not

num_SDS=1; % number of detectors
%     SDS_x_arr=[0.8 1.5 1.5 3 3 4.5 4.5]; % the SDS x displacement (cm)
%     SDS_z_arr=[0 0 -1.5 0 -1.5 0 -1.5]; % the SDS z displacement (cm)
SDS_x_arr=[0 0.8]; % the SDS x displacement (cm) to Fp2, the 1st is the source
SDS_z_arr=[0.75 0.75]; % the SDS z displacement (cm) to Fp2, the 1st is the source
reference_index_arr=[2 3]; % the reference EEG point for each detector on the head, ['Fp2h','Fp2','AFp8','AF8','AFF8']

%% main
for sbj=1:length(subject_name_arr)
    subject_name=subject_name_arr{sbj};
    MRI_voxel_file=fullfile(model_folder,['headModel' subject_name '_EEG.mat']); % containing the MRI voxel model, also the EEG point and head surface mesh
    
    MRI_model=load(MRI_voxel_file);
    if isfield(MRI_model,'model_version')==0
        MRI_model.model_version=1;
    end
    
    %% find the position and direction of probes
    [p_pos,p_dir]=fun_cal_probe_position_and_direction_2D(MRI_model.EEG,SDS_x_arr,SDS_z_arr,MRI_model.vol,MRI_model.headsurf,MRI_model.voxel_size,reference_index_arr,MRI_model.model_version);
    
    save(fullfile(model_folder,[subject_name '_probe_pos.txt']),'p_pos','-ascii','-tabs');
    save(fullfile(model_folder,[subject_name '_probe_dir.txt']),'p_dir','-ascii','-tabs');
    
    
    %% plot the probes
    if to_plot_figure    
        model_pic_dir=fullfile(model_folder,'model_pic');
        if exist(model_pic_dir,'dir')==0
            mkdir(model_pic_dir);
        end
        plot_probe_position(MRI_model.vol,MRI_model.headsurf,num_SDS,p_pos,p_dir,[],[],[],[],model_pic_dir,subject_name,MRI_model.model_version);
    end
    close all;
    
end

disp('Done!');

%% functions

%% plot the 2D and 3D plot for the probe position to make sure the probe is on the subject's head surface
function plot_probe_position(vol,headsurf,num_SDS,p_pos,p_dir,src_pos,src_dir,set_r,set_angle,model_folder,subject_name,model_version)
%% param
do_plot_2D=1;
do_plot_3D=1;
do_plot_disk_source=0;
do_plot_cone_source=0;
do_SDS_text=1;
cone_hight=40;

dir_len=10;

%% main

image_ratio=size(vol,1)/size(vol,2);
if image_ratio>1
    image_ratio=1/image_ratio;
end


if do_plot_2D
    for s=1:num_SDS+1
        figure('Units','pixels','position',[0 0 round(image_ratio*800) 800]);
        imagesc(vol(:,:,round(p_pos(s,3))));
        hold on;
        if model_version==1
            plot(p_pos(s,1),p_pos(s,2),'xr','LineWidth',2,'MarkerSize',24);
            quiver(p_pos(s,1)-p_dir(s,1)*dir_len,p_pos(s,2)-p_dir(s,2)*dir_len,p_dir(s,1),p_dir(s,2),dir_len,'k','LineWidth',3,'MarkerSize',24);
        elseif model_version>=2
            plot(p_pos(s,2),p_pos(s,1),'xr','LineWidth',2,'MarkerSize',24); % swap the x and y
            quiver(p_pos(s,2)-p_dir(s,2)*dir_len,p_pos(s,1)-p_dir(s,1)*dir_len,p_dir(s,2),p_dir(s,1),dir_len,'k','LineWidth',3,'MarkerSize',24); % swap the x and y
        end
        if s==1
            title(['source at layer ' num2str(round(p_pos(s,3)))]);
        else
            title(['SDS ' num2str(s-1) ' at layer ' num2str(round(p_pos(s,3)))]);
        end
        saveas(gcf,fullfile(model_folder,[subject_name '_SDS_' num2str(s-1) '.png']));
        close all;
    end
end

if do_plot_3D
    %% plot the head and the SDSs
    figure('Units','pixels','position',[0 0 1400 1080]);
    patch(headsurf,'FaceColor',[1,.75,.65],'EdgeColor','none','FaceAlpha','0.9');
    lighting gouraud
    lightangle(0,30);
    lightangle(120,30);
    lightangle(-120,30);
    axis('image');
    hold on;
    plot3(p_pos(:,1),p_pos(:,2),p_pos(:,3),'ro','MarkerSize',10,'LineWidth',2);
    for s=1:size(p_pos,1)
        if s==1
            quiver3(p_pos(s,1)-p_dir(s,1)*dir_len,p_pos(s,2)-p_dir(s,2)*dir_len,p_pos(s,3)-p_dir(s,3)*dir_len,p_dir(s,1)*dir_len,p_dir(s,2)*dir_len,p_dir(s,3)*dir_len,'k','LineWidth',3,'MarkerSize',24);
            if do_SDS_text
                text(p_pos(s,1)-p_dir(s,1)*dir_len,p_pos(s,2)-p_dir(s,2)*dir_len,p_pos(s,3)-p_dir(s,3)*dir_len,'source','FontSize',20,'FontName','Times New Roman','Color','g');
            end
        else
            quiver3(p_pos(s,1),p_pos(s,2),p_pos(s,3),p_dir(s,1)*dir_len,p_dir(s,2)*dir_len,p_dir(s,3)*dir_len,'k','LineWidth',3,'MarkerSize',24);
            if do_SDS_text
                text(p_pos(s,1)+p_dir(s,1)*dir_len,p_pos(s,2)+p_dir(s,2)*dir_len,p_pos(s,3)+p_dir(s,3)*dir_len,['SDS ' num2str(s-1)],'FontSize',20,'FontName','Times New Roman','Color','g');
            end
        end
    end
    
    if do_plot_disk_source
        %% generate cylinder for disk source
        [xx,yy,zz]=cylinder(set_r,100);
        zz=zz.*cone_hight;
        
        % find the angle to turn
        orig_cyl_dir=[0 0 1];
        phy=acos(src_dir(3)); % altitude angle
        orig_cyl_dir=[cos(phy) 0 sin(phy);0 1 0;-sin(phy) 0 cos(phy)]*orig_cyl_dir'; % rotate on x-z plan
        theta=acos(src_dir(1)./sqrt(sum(src_dir(1:2).^2))); % azimuth angle
        if model_version>=2
            theta=-theta;
        end
        orig_cyl_dir=[cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]*orig_cyl_dir;  % rotate on x-y plan
        orig_cyl_dir=orig_cyl_dir'; % this should be the same as p_dir
        
        temp_cyl_arr=[reshape(xx,[],1) reshape(yy,[],1) reshape(zz,[],1)];
        temp_cyl_arr=[cos(phy) 0 sin(phy);0 1 0;-sin(phy) 0 cos(phy)]*temp_cyl_arr';
        temp_cyl_arr=[cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]*temp_cyl_arr;
        temp_cyl_arr=temp_cyl_arr';
        xx=reshape(temp_cyl_arr(:,1),size(xx))+src_pos(1);
        yy=reshape(temp_cyl_arr(:,2),size(yy))+src_pos(2);
        zz=reshape(temp_cyl_arr(:,3),size(zz))+src_pos(3);
        
        surf(xx,yy,zz);
    end
    
    if do_plot_cone_source
        [xx,yy,zz]=cylinder([0 cone_hight*tan(set_angle/180*pi)],100);
        zz=zz*cone_hight;
        
        % find the angle to turn
        orig_cyl_dir=[0 0 1];
        phy=acos(src_dir(3)); % altitude angle
        orig_cyl_dir=[cos(phy) 0 sin(phy);0 1 0;-sin(phy) 0 cos(phy)]*orig_cyl_dir'; % rotate on x-z plan
        theta=acos(src_dir(1)./sqrt(sum(src_dir(1:2).^2))); % azimuth angle
        if model_version>=2
            theta=-theta;
        end
        orig_cyl_dir=[cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]*orig_cyl_dir;  % rotate on x-y plan
        orig_cyl_dir=orig_cyl_dir'; % this should be the same as p_dir
        
        temp_cyl_arr=[reshape(xx,[],1) reshape(yy,[],1) reshape(zz,[],1)];
        temp_cyl_arr=[cos(phy) 0 sin(phy);0 1 0;-sin(phy) 0 cos(phy)]*temp_cyl_arr';
        temp_cyl_arr=[cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]*temp_cyl_arr;
        temp_cyl_arr=temp_cyl_arr';
        xx=reshape(temp_cyl_arr(:,1),size(xx))+src_pos(1);
        yy=reshape(temp_cyl_arr(:,2),size(yy))+src_pos(2);
        zz=reshape(temp_cyl_arr(:,3),size(zz))+src_pos(3);
        
        surf(xx,yy,zz);
    end
    
    view([230 0]);
    saveas(gcf,fullfile(model_folder,['headModel' subject_name '_fiber_plot.fig']));
end
end