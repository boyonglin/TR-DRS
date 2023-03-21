%{
Given some position of EEG points, and find the position and direction of optical probe
Also make sure the source probe is outside the head model

assert the source probe is at FP2 point
input: 
EEG: a structure, the 3-D point of some EEG points
SDS_x_arr: the array of SDS in x direction in cm
SDS_z_arr: the array of SDS in z direction in cm
vol: the 3D head voxel, each voxel is a number, representing the tissue type of this voxel.  0 for fiber.
voxel_size: the length of side of each voxel in the MRI voxel model, in mm
reference_index_arr: the index of reference point for each SDS
model_version: The version of the model.  The model orientation in version 1 and 2 are different.

output:
p_pos: the array of probe position
p_dir: the array of probe direction

Benjamin Kao
Last update: 2020/11/26
%}


function [p_pos, p_dir]=fun_cal_probe_position_and_direction_2D(EEG,SDS_x_arr,SDS_z_arr,vol,headsurf,voxel_size,reference_index_arr,model_version)
%% param
do_shift_source_outside=1;
do_plot_ref=0;
do_plot_circle=1; % plot the fitted circle and slice

%% init
SDS_x_arr=SDS_x_arr./voxel_size; % scale the distance to exactly how many voxel in length
SDS_z_arr=SDS_z_arr./voxel_size;

%% find the reference points
global Fp2_index;
Fp2_index=2; % which column of reference point is Fp2
reference_point_arr={'Fp2h','Fp2','AFp8','AF8','AFF8'; ...
                      'NFp2h','NFp2','AFp10h','AF10h','AFF10h'}; % the reference points to find the location of the probes
other_reference_point_arr={'Fp1','AFp7','AF7','AFF7','NFp1','AFp9h','AF9h','AFF9h'}; % the reference points used to find the sphere

% plot the reference points
if do_plot_ref
    figure('units','normalized','outerposition',[0 0 1 1]);
    hold on;
    axis([0,5,0,3]);
    plot([1 4],[1 1],'b--');
    plot([1 4],[2 2],'b--');
    for i=1:4
        plot([i i],[1 2],'b--');
    end
    for i=1:4
        for j=1:2
            text(i,3-j,reference_point_arr{j,i},'FontSize',20);
        end
    end

    for j=1:2
        for i=1:3
        %     point_1_pos=getfield(EEG,reference_point_arr{1,i});
        %     point_2_pos=getfield(EEG,reference_point_arr{1,i+1});
            point_1_pos=EEG.(reference_point_arr{j,i});
            point_2_pos=EEG.(reference_point_arr{j,i+1});
    %         fprintf('distance between %s and %s is %.2f mm\n',reference_point_arr{j,i},reference_point_arr{j,i+1},sqrt(sum((point_1_pos-point_2_pos).^2)));
            text(i+0.5,3-j+0.05,num2str(sqrt(sum((point_1_pos-point_2_pos).^2))),'FontSize',12);
        end
    end
    for i=1:4
        point_1_pos=EEG.(reference_point_arr{1,i});
        point_2_pos=EEG.(reference_point_arr{2,i});
        text(i+0.1,1.5,num2str(sqrt(sum((point_1_pos-point_2_pos).^2))),'FontSize',12);
    end
    axis off;
    title('reference points and distance');
    hold off;
end


%% use fitting to find the center and radius of the circle
global circle_points;
circle_points=[];
for j=1:size(reference_point_arr,1)
    for i=1:size(reference_point_arr,2)
        circle_points=[circle_points; EEG.(reference_point_arr{j,i})];
    end
end
% add other reference points
for i=1:size(other_reference_point_arr,2)
    circle_points=[circle_points; EEG.(other_reference_point_arr{1,i})];
end
if model_version==1
    x0=[circle_points(Fp2_index,1)-85 circle_points(Fp2_index,2)+5 circle_points(Fp2_index,3)-5 90]; % [x,y,z] and r of the circle
elseif model_version>=2
    x0=[circle_points(Fp2_index,1)+60 circle_points(Fp2_index,2)-25 circle_points(Fp2_index,3)-2 90]; % [x,y,z] and r of the circle
end
[x]=fminsearch(@circle_distance,x0);

fprintf('the radius of fitting circle is %.2f mm\n',x(4));
fprintf('the center is at %.2f, %.2f, %.2f\n',x(1),x(2),x(3));
fprintf('distance of center to reference points:\n');
for j=1:size(reference_point_arr,1)
    for i=1:size(reference_point_arr,2)
        fprintf('\t%.2f,',sqrt(sum((x(1:3)-circle_points(i+(j-1)*size(reference_point_arr,2),:)).^2)));
    end
end
fprintf('\n');

%% plot the fitted circle points
if do_plot_circle
    figure();
    imagesc(vol(:,:,round(x(3))));
    hold on;
    if model_version==1
        plot(circle_points(:,1),circle_points(:,2),'xr','MarkerSize',20); % switch x and y
        plot(x(1),x(2),'xr','MarkerSize',20);
    elseif model_version>=2
        plot(circle_points(:,2),circle_points(:,1),'xr','MarkerSize',20); % switch x and y
        plot(x(2),x(1),'xr','MarkerSize',20);
    end
end

%% find the direction of probe
% because the SDS is not always near the reference points, so I didn't use
% for loop to process each SDS
global center radius SDS_arr reference_index SDS_index
center=x(1:3);
radius=x(4);
SDS_arr=SDS_x_arr;

%% for source
% p_dir(1,:)=center-circle_points(Fp2_index,:);
% p_dir(1,:)=p_dir(1,:)./sqrt(sum(p_dir(1,:).^2));
% p_pos(1,:)=circle_points(Fp2_index,:);
% 
% % make sure the source probe is out of the head
% if do_shift_source_outside==1
%     shift_length=0.1;
%     shift_time=0;
%     old_p_pos1=p_pos(1,:);
%     while vol(ceil(p_pos(1,2)),ceil(p_pos(1,1)),ceil(p_pos(1,3)))~=0 % swap the x and y
%         p_pos(1,:)=p_pos(1,:)-p_dir(1,:)*shift_length;
%         shift_time=shift_time+1;
%     end
%     if shift_time~=0
%         fprintf('The Source Probe is not in air, shift it %d * %.3f, that is %.3f, to the air.\n',shift_time,shift_length,sqrt(sum((p_pos(1,:)-old_p_pos1).^2)));
%         fprintf('\tFrom %.2f, %.2f, %.2f to %.2f, %.2f, %.2f\n',old_p_pos1(1),old_p_pos1(2),old_p_pos1(3),p_pos(1,1),p_pos(1,2),p_pos(1,3));
%     end
% end


%% for detectors

for s=1:length(SDS_arr)
    reference_index=reference_index_arr(s);
    SDS_index=s;
    if SDS_x_arr(s)~=0
        alpha=1; % init value, the ratio of probe to Fp2 compare to reference point to Fp2
        alpha=fminsearch(@find_ratio,alpha);
        p_pos_prime=circle_points(Fp2_index,:)+alpha*(circle_points(reference_index,:)-circle_points(Fp2_index,:)); % calculate the approximate point on the line of reference EEG points
        p_dir(s,:)=p_pos_prime-center; % calculate the direction from the center to the approximate point
        p_dir(s,:)=p_dir(s,:)./sqrt(sum(p_dir(s,:).^2)); % normalize the direction vector
    else
        p_dir(s,:)=circle_points(Fp2_index,:)-center;
        p_dir(s,:)=p_dir(s,:)./sqrt(sum(p_dir(s,:).^2)); % normalize the direction vector
    end
    if SDS_z_arr(s)~=0
        angle=SDS_z_arr(s)*10/radius; % calculate the angle change in z direction, also turn cm into mm
        temp_point=center+radius*p_dir(s,:); % the temp point (z shift = 0)
        FR=circle_points(1,:)-circle_points(Fp2_index,:); % the vector of Fp2 to Fp2h
        cross_CP_FR=[p_dir(s,2)*FR(3)-p_dir(s,3)*FR(2) p_dir(s,3)*FR(1)-p_dir(s,1)*FR(3) p_dir(s,1)*FR(2)-p_dir(s,2)*FR(1)]; % the corss of [center to point] and [Fp2 to Fp2h];
        cross_CP_FR=cross_CP_FR./sqrt(sum(cross_CP_FR.^2));
        temp_point_to_true_point=cross_CP_FR*tan(angle)*radius; % the vector from the temp point (z shift = 0) to the true point
        temp_point=temp_point+temp_point_to_true_point;
        p_dir(s,:)=temp_point-center;
        p_dir(s,:)=p_dir(s,:)./sqrt(sum(p_dir(s,:).^2)); % normalize the direction vector
    end
    p_pos(s,:)=center+radius*p_dir(s,:); % extend the direction to the radius of the circle
    if s==1 % flip the direction of the source probe
        p_dir(s,:)=-1*p_dir(s,:);
    end
    
    % adjust the probe position to she surface
    dist_2_vertices=sum((headsurf.vertices-p_pos(s,:)).^2,2);
    [dist_2_vertices,closest_index]=sort(dist_2_vertices);
    to_surf_threshold=1;
    if min(dist_2_vertices)>to_surf_threshold
        % move it toward closest point
        num_pt_to_interp=3;
        v_point=headsurf.vertices(closest_index(1:num_pt_to_interp),:);
        dist_2_vertices=dist_2_vertices(1:num_pt_to_interp);
        new_pos=[0 0 0];
        for pt=1:num_pt_to_interp
            new_pos=new_pos+v_point(pt,:)./dist_2_vertices(pt);
        end
        new_pos=new_pos./sum(1./dist_2_vertices);
        if s>1
            fprintf('\tMove SDS %d %2f voxel to surface.\n',s-1,sqrt(sum((p_pos(s,:)-new_pos).^2)));
        else
            fprintf('\tMove source %2f voxel to surface.\n',sqrt(sum((p_pos(s,:)-new_pos).^2)));
        end
        p_pos(s,:)=new_pos;
    end
    
    if do_shift_source_outside==1 && s==1
        shift_length=0.1;
        shift_time=0;
        old_p_pos1=p_pos(1,:);
        if model_version==1
            while vol(ceil(p_pos(1,2)),ceil(p_pos(1,1)),ceil(p_pos(1,3)))~=0 % swap the x and y
                p_pos(1,:)=p_pos(1,:)-p_dir(1,:)*shift_length;
                shift_time=shift_time+1;
            end
        elseif model_version==2
            while vol(ceil(p_pos(1,1)),ceil(p_pos(1,2)),ceil(p_pos(1,3)))~=0
                p_pos(1,:)=p_pos(1,:)-p_dir(1,:)*shift_length;
                shift_time=shift_time+1;
            end
        end
        if shift_time~=0
            fprintf('The Source Probe is not in air, shift it %d * %.3f, that is %.3f, to the air.\n',shift_time,shift_length,sqrt(sum((p_pos(1,:)-old_p_pos1).^2)));
            fprintf('\tFrom ( %.2f, %.2f, %.2f ) to ( %.2f, %.2f, %.2f )\n',old_p_pos1(1),old_p_pos1(2),old_p_pos1(3),p_pos(1,1),p_pos(1,2),p_pos(1,3));
        end
    end
end

%% report the distance
for s=2:length(SDS_x_arr)
    fprintf('Distance of SDS %d to source = %f mm\n',s-1,sqrt(sum((p_pos(1,:)-p_pos(s,:)).^2,2))*voxel_size);
    if SDS_z_arr(s)~=SDS_z_arr(1)
        fprintf('Distance of SDS %d to SDS %d = %f mm\n',s-1,s-2,sqrt(sum((p_pos(s-1,:)-p_pos(s,:)).^2,2))*voxel_size);
    end
end

end

function error=circle_distance(x)
global circle_points;
error=0;
for i=1:size(circle_points,1)
    error=error+(sqrt(sum((circle_points(i,:)-x(1:3)).^2))-x(4))^2;
end
end

function error=find_ratio(alpha)
global circle_points center radius SDS_arr reference_index SDS_index Fp2_index
dir_to_fp2=circle_points(Fp2_index,:)-center; % the dir from the center of sphere to Fp2
probe_pos=circle_points(Fp2_index,:)+alpha*(circle_points(reference_index,:)-circle_points(Fp2_index,:)); % calculate the probe position
dir_to_probe=probe_pos-center; % the dir from the center of sphere to the probe position

target_angle=SDS_arr(SDS_index)*10/(2*pi*radius)*2*pi; % turn the SDS cm into mm
cos_theta=sum(dir_to_probe.*dir_to_fp2)/sqrt(sum(dir_to_probe.^2))/sqrt(sum(dir_to_fp2.^2));
error=abs(abs(target_angle)-acos(cos_theta));
end