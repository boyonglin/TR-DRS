%{
Run MCX to find the fluence rate or reflectance
dist2axis: use the distance to the axis of the detector to determine if a
photon is detected or not.

Benjamin Kao
last update: 2020/06/21
Version: 4.41

input:
tissue_param: the optical parameter of tissue in 1/cm, a n*m matrix, with n is the wavelength number, m = 2*(number of tissue types), [mua_1 mus_1 mua_2 mus_2......]
head_vol: the 3D head voxel, each voxel is a number, representing the tissue type of this voxel.  0 for fiber.
voxel_size: the length of side of each voxel in the MRI voxel model, in mm
p_pos: the position of the probes, including the source probe, a (s+1)*3 matrix
p_dir: the direction of the probes, including the source probe, a (s+1)*3 matrix
sim_set: the simulation setting structure
output_folder: the folder to store the simulation output
output_fluence_rate: =1 to output the fluence rate, which is get from all the incident photon.
output_pathlength: =1 to output the pathlength of detected photons
output_jacobian: =1 to output the jacobian of detected photons
do_checkpoint: =1 to store checkpoint, each checkpoint store the current detected photon infomation. =0 for no checkpoint.
max_CV: if the CV of detected photon is smaller than this value, reduce the detected photon number to reduce size. =0 to turn off this feature.
min_detpt: if the photon is larger than this value, try do reduce them.  Work with max_CV.
%}

function fun_MCX_sim_dist2axis(tissue_param,head_vol,voxel_size,p_pos,p_dir,sim_set,output_folder,output_fluence_rate,output_pathlength,output_jacobian,do_checkpoint,max_CV,min_detpt)

%% param
% about simulation
show_MCX_info=0; % show the MCX default run information

num_shuffle_time=10; % if check CV, do this times to average the CV
num_CV_group=5; % if check CV, how many groups should the photon be devided.

% about tissue and instrument
num_photon=sim_set.num_photon;
photon_per_simulation=sim_set.photon_per_simulation; % howmany photon to run each time, so not cost too much memory simultaneously
mcx_max_detpt=sim_set.mcx_max_detpt;
num_SDS=sim_set.num_SDS;
num_layer=sim_set.num_layer;
to_output_layer=sim_set.to_output_layer; % the layers to output pathlength
detector_r=sim_set.detector_r./voxel_size; %mm
detector_larger_r=sim_set.detector_larger_r./voxel_size; %mm, use the larger detector to find the detected photon, and calculate the distance of photon to the center of detector
detector_NA=sim_set.detector_NA; % the numerical aperture of the fibers
source_type=sim_set.source_type; % the type of source, 'cone' or 'disk' or 'pencil'
source_r=sim_set.source_r./voxel_size; % mm, the radius of the source
source_NA=sim_set.source_NA; % the numerical aperture of the source
n=sim_set.n; % the refractive index of tissue and outer air/(probe surface)
g=sim_set.g; % the anisotropic factor of tissue and outer air/(probe surface)
fiber_n=sim_set.fiber_n; % the refractive index of fibers

%% initialize
global sim_version;
assert(sim_version==4.41,'ERROR: The version of main and function is not consistent!');

assert(sum([output_fluence_rate output_pathlength])>0,'ERROR: There will be no output!');
assert(output_jacobian<=output_pathlength,'ERROR: cannot output jacobian without replay!');

assert(sum(size(p_pos)==size(p_dir))==2,'ERROR: size of p_pos and p_dir not match!');
assert(size(p_pos,1)==num_SDS+1,'ERROR: num_SDS not match with probe setting!');
assert(size(p_pos,1)==size(detector_r,1)+1,'ERROR: detector_r not match with probe setting!')
assert(size(p_pos,1)==size(detector_larger_r,1)+1,'ERROR: detector_larger_r not match with probe setting!')
assert(size(p_pos,1)==size(detector_NA,1)+1,'ERROR: detector_NA not match with probe setting!')
assert(max(max(max(head_vol)))<=num_layer,'ERROR: num_layer and volume not match!');
assert(max(max(max(head_vol)))*2<=size(tissue_param,2),'ERROR: tissue_param and volume not match!');

assert(max_CV>=0,'ERROR: max_CV setting error!');
if max_CV>0
    assert(min_detpt>0,'ERROR: min_detpt not setting!');
end


if num_photon<photon_per_simulation
    photon_per_simulation=num_photon;
    num_need_to_run=1;
else
    assert(rem(num_photon,photon_per_simulation)==0,'num photon is not the multiple of num per simulation!');
    num_need_to_run=ceil(num_photon/photon_per_simulation); % how many times should run
end

if exist(output_folder,'dir')==0
    mkdir(output_folder);
end
if exist(fullfile(output_folder,'sim_set.mat'),'file')==0
    save(fullfile(output_folder,'sim_set.mat'),'sim_set');
end

%% prepare probe information for MCX simulation
p_pos(:,[1 2])=p_pos(:,[2 1]); % exchange the x, y pos
p_dir(:,[1 2])=p_dir(:,[2 1]); % exchange the x, y pos
SDS_normal_vector=p_dir(2:end,:); % the normal vector of the detectors

%% set the simulation
cfg.nphoton=photon_per_simulation;
cfg.vol=uint8(head_vol);
cfg.srcpos=p_pos(1,:);
cfg.srcdir=p_dir(1,:);
if strcmp(source_type,'pencil')
    assert(source_NA==0,'source angle setting ERROR!');
    assert(source_r==0,'source radius setting ERROR!');
    cfg.srctype='pencil'; % source type
elseif strcmp(source_type,'cone')
    assert(source_r==0,'source radius setting ERROR!');
    cfg.srctype='cone';
    sin_angle=source_NA/fiber_n;
    assert(sin_angle>=0,'source NA setting ERROR!');
    assert(sin_angle<=1,'source NA setting ERROR!');
    source_angle=asin(sin_angle);
    cfg.srcparam1=[source_angle 0 0 0];
    cfg.srcparam2=[0 0 0 0];
elseif strcmp(source_type,'disk')
    assert(source_NA==0,'source angle setting ERROR!');
    cfg.srctype='disk';
    assert(source_r>=0,'source radius setting ERROR!');
    cfg.srcparam1=[source_r 0 0 0];
    cfg.srcparam2=[0 0 0 0];
end

cfg.detpos=p_pos(2:end,:);
cfg.detpos(:,4)=detector_larger_r;


% make a reference point on the axis of detector
reference_point_arr=zeros(size(cfg.detpos,1),3);
for s=1:size(cfg.detpos,1)
    reference_point_arr(s,:)=cfg.detpos(s,1:3)-SDS_normal_vector(s,:).*detector_larger_r(s);
end

%% for GPU using
sim_summary.sim_GPU='';

[~,gpuinfo]=evalc('mcxlab(''gpuinfo'')');

if length(gpuinfo)==1
    cfg.gpuid=1; % =1 use the first GPU, =2, use the 2nd, ='11' use both GPUs together
    sim_summary.sim_GPU=gpuinfo(1).name;
else
    GPU_setting=load('GPU_setting.txt');
    assert(sum(GPU_setting(1,:))>=1,'ERROR: No GPU selection!');
    assert(length(gpuinfo)>=max(find(GPU_setting(1,:)==1)),'ERROR: Select GPU not exist!');
    if sum(GPU_setting(1,:))==1
        cfg.gpuid=find(GPU_setting(1,:)==1);
        sim_summary.sim_GPU=gpuinfo(cfg.gpuid).name;
    else
        gpuid_arr='';
        workload_arr=[];
        for gp=1:length(gpuinfo)
            gpuid_arr(1,gp)=num2str(GPU_setting(1,gp),'%d');
            workload_arr(1,gp)=GPU_setting(2,gp);
            if GPU_setting(1,gp)==1
                sim_summary.sim_GPU=[sim_summary.sim_GPU gpuinfo(gp).name ','];
            end
        end
        cfg.gpuid=gpuid_arr;
        cfg.workload=workload_arr;
    end
end


%% other setting for MCX
cfg.tstart=0;                        % start of the simulation time window (in second)
cfg.tend=5e-9;                       % end of the simulation time window (in second)
cfg.tstep=2e-10;                     % time gate width (in second), here we asks mcxlab to a "videos" at 50 time gates
cfg.isnormalized=1; % normalize the output fluence to unitary source
cfg.unitinmm=voxel_size; % the length unit for a grid edge length

% cfg.savedetflag='dspxvw';
cfg.savedetflag='dpxv';
% 1 d  output detector ID (1)
% 2 s  output partial scat. even counts (#media)
% 4 p  output partial path-lengths (#media)
% 8 m  output momentum transfer (#media)
% 16 x  output exit position (3)
% 32 v  output exit direction (3)
% 64 w  output initial weight (1)
cfg.outputtype='energy';
cfg.maxdetphoton=mcx_max_detpt; % maximum number of photons saved by the detectors
cfg.isgpuinfo=0; % 1-print GPU info, [0]-do not print
cfg.autopilot=1;  % let mcxlab to automatically decide the threads/blocks
cfg.issrcfrom0=1; % first voxel is [1 1 1]
cfg.issaveexit=1; % save the exiting photons
cfg.issaveref=1;
cfg.isreflect=0; % consider index mismatch reflect
cfg.isspecular=1; % calculate specular reflection if source is outside

reflectance_arr=[];

wl_start_index=1;
from_chk_flag=0;
if do_checkpoint
    if exist(fullfile(output_folder,'checkpioint_index.txt'),'file')~=0
        from_chk_flag=1;
        wl_runin_chk=load(fullfile(output_folder,'checkpioint_index.txt'));
        fprintf('\nStart from checkpoint wl: %d, subsim: %d \n',wl_runin_chk(1),wl_runin_chk(2));
        if wl_runin_chk(2)>0 % if some subsim for this wl had been simed
            old_simtime_chk=load(fullfile(output_folder,'simtime_chk.txt'));
        else
            old_simtime_chk=0;
        end
        % check the runned chackpoint file
        if output_pathlength || output_jacobian
            for chkidx=1:wl_runin_chk(2)
                assert(exist(fullfile(output_folder,['detpt_chk_' num2str(chkidx) '.mat']),'file')~=0,['ERROR: checkpoint ' num2str(chkidx) ' has loss!']);
            end
            if wl_runin_chk(2)>0 % if some subsim for this wl had been simed
                old_numpt_chk=load(fullfile(output_folder,'numpt_chk.txt'));
            else
                old_numpt_chk=zeros(1,num_SDS);
            end
        end
        if output_fluence_rate
            if wl_runin_chk(2)>0 % if some subsim for this wl had been simed
                old_flu_chk=load(fullfile(output_folder,'fluencerate_chk.mat'));
                old_flu_chk=old_flu_chk.fluencerate_arr;
            else
                old_flu_chk=[];
            end
        end
        wl_start_index=wl_runin_chk(1);
    end
end

for wl=wl_start_index:size(tissue_param,1)
    fprintf('sim wl %d/%d subsim: ,',wl,size(tissue_param,1));
    
    %% start timer
    sim_summary.num_photon=num_photon;
    
    rng('default');
    
    %% set the optical parameter
    cfg.prop=[0 0 1 fiber_n];
    for l=1:num_layer
        cfg.prop(l+1,:)=[tissue_param(wl,2*l-1) tissue_param(wl,2*l) g(l) n(l)];
    end
    cfg.prop(:,1:2)=cfg.prop(:,1:2).*0.1; % turn 1/cm into 1/mm
    
    %% init array to store the detected photon from each detector
    if output_jacobian
        seed_arr=cell(1,num_SDS);
        detpt_arr=cell(1,num_SDS);
        for d=1:num_SDS
            detpt_arr{d}=[];
        end
    end
    if output_pathlength
        detPL_arr=cell(1,num_SDS);
        for d=1:num_SDS
            detPL_arr{d}=[];
        end
    end
    
    if max_CV>0
        if exist(fullfile(output_folder,'cheCV_chk.mat'),'file') && from_chk_flag
            CV_chk=load(fullfile(output_folder,'cheCV_chk.mat'));
            detpt_weight_arr=CV_chk.detpt_weight_arr;
            SDS_CV_arr=CV_chk.SDS_CV_arr;
            each_photon_weight_arr=CV_chk.each_photon_weight_arr;
            SDS_done_flag=CV_chk.SDS_done_flag;
        else
            detpt_weight_arr=cell(1,num_SDS); % the weight of each detected photon.
            SDS_CV_arr=ones(1,num_SDS)*inf; % the CV if this SDS had been check
            each_photon_weight_arr=zeros(1,num_SDS); % the number of simulated photon for this SDS.
            SDS_done_flag=zeros(1,num_SDS); % whether the SDS is done due to small enough CV
        end
    else
        SDS_done_flag=zeros(1,num_SDS); % whether the SDS is done due to small enough CV
    end
    
    
    if from_chk_flag==0
        sim_time_calc=0;
        fluencerate_arr=[];
        if output_jacobian || output_pathlength
            SDS_detected_numpt=zeros(1,num_SDS);
        end
        run_index=0;
    else
        sim_time_calc=old_simtime_chk;
        if output_fluence_rate
            fluencerate_arr=old_flu_chk;
        end
        if output_jacobian || output_pathlength
            SDS_detected_numpt=old_numpt_chk;
        end
        run_index=wl_runin_chk(2);
    end
        
    %% run first time to get the random seed for detected photon
    
    
    while run_index<num_need_to_run
        timer_sim=tic;
        fprintf(' %d,',run_index+1);
        cfg.seed=randi([1,1000000]);
        save(fullfile(output_folder,['cfg_' num2str(wl) '.mat']),'cfg');
        if show_MCX_info
            [fluencerate,detpt,~,seeds]=mcxlab(cfg);
        else
            [~,fluencerate,detpt,~,seeds]=evalc('mcxlab(cfg);');
        end
        
        if output_pathlength || output_jacobian
            if isfield(detpt,'detid') % to avoid that there are no photon detected in this simulation
                assert(length(detpt.detid)<cfg.maxdetphoton,'Error: detected too much photon!'); % to avoid the loss due to too many detected photon
                for s=1:num_SDS
                    if SDS_done_flag(s)==1
                        continue;
                    end
                    %% find the photon in the true detector r
                    dist2reference_square=sum((detpt.p-reference_point_arr(s,:)).^2,2); % the distance of Photon and Reference
                    dot_RP_RD=sum((detpt.p-reference_point_arr(s,:)).*SDS_normal_vector(s,:),2); % the cosine value of RP X RD (D for Detector). RD is equal to the normal vetector * dist(RD)
                    % because the length of normal vector is 1, so the dot value is equal to the RP * cos(theta)
                    dist2axis_square=dist2reference_square-dot_RP_RD.^2;
                    in_dist_detpt=find(dist2axis_square<=detector_r(s)^2);
                    
                    %% find the photon in the NA
                    if length(in_dist_detpt)>0
                        det_vector=detpt.v(in_dist_detpt,:);
                        dot_value=sum(SDS_normal_vector(s,:).*det_vector,2);
                        dot_value(dot_value>1)=1;
                        det_angle=acos(dot_value);
                        in_NA_index=in_dist_detpt(dot_value>0 & sin(det_angle)<=(detector_NA(s)/fiber_n));
                        SDS_detected_numpt(s)=SDS_detected_numpt(s)+length(find(in_NA_index));
                        if output_pathlength
                            detPL_arr{s}=[detPL_arr{s}; detpt.ppath(in_NA_index,to_output_layer)./10.*voxel_size]; % turn voxel length (mm) into cm
                        end
                        if output_jacobian
                            seed_arr{s}=[seed_arr{s} seeds.data(:,in_NA_index)]; % store the detected photon random seed
                            detpt_arr{s}=[detpt_arr{s} detpt.data(:,in_NA_index)]; % store the detected photon data
                        end
                    end
                end
            end
        end
        
        clear detpt seeds;
        run_index=run_index+1;
        if output_fluence_rate
            if run_index==1
                fluencerate_arr=fluencerate.data;
            else
                fluencerate_arr=fluencerate_arr+fluencerate.data;
            end
        end
        
        sim_time_calc=sim_time_calc+toc(timer_sim);
        
        %% check CV
        if max_CV>0 && output_pathlength
            temp_mua_arr=tissue_param(wl,1:2:end); % in 1/cm
            temp_mua_arr=temp_mua_arr(to_output_layer);
            for s=1:num_SDS
                if SDS_done_flag(s)==0
                    if size(detPL_arr{s},1)>0 %% add the detpt from this subsim to weight arr
                        if do_checkpoint~=0
                            temp_weight_arr=exp(-1*sum(detPL_arr{s}.*temp_mua_arr,2));
                            detpt_weight_arr{s}(end+1:end+length(temp_weight_arr),1)=temp_weight_arr;
                        elseif do_checkpoint==0
                            detpt_weight_arr{s}=exp(-1*sum(detPL_arr{s}.*temp_mua_arr,2));
                        end
                    end
                    
                    if SDS_detected_numpt(s)>min_detpt % only check CV for large enough photon number
                        photon_in_subGroup=ones(1,num_CV_group)*floor(size(detpt_weight_arr{s},1)/num_CV_group);
                        photon_in_subGroup(1:rem(size(detpt_weight_arr{s},1),num_CV_group))=photon_in_subGroup(1:rem(size(detpt_weight_arr{s},1),num_CV_group))+1;
                        assert(sum(photon_in_subGroup)==size(detpt_weight_arr{s},1));
                        photon_in_subGroup=cumsum([0 photon_in_subGroup]);
                        CV_ref_subGroup=zeros(num_CV_group,num_shuffle_time);
                        for sh=1:num_shuffle_time
                            shuffled_index=randperm(size(detpt_weight_arr{s},1));
                            for sg=1:num_CV_group
                                CV_ref_subGroup(sg,sh)=sum(detpt_weight_arr{s}(shuffled_index(photon_in_subGroup(sg)+1:photon_in_subGroup(sg+1))));
                            end
                        end
                        CV_ref_subGroup=std(CV_ref_subGroup)./mean(CV_ref_subGroup,1);
                        temp_CV=sqrt(mean(CV_ref_subGroup.^2));
                        if temp_CV<=max_CV
                            detpt_weight_arr{s}=[]; % clear this array to save ckeck point time.
                            SDS_CV_arr(s)=temp_CV;
                            SDS_done_flag(s)=1;
                            each_photon_weight_arr(s)=run_index*photon_per_simulation;
                            fprintf('\tSDS %d has CV %f%%, early stop!\n',s,temp_CV*100);
                        end
                    end
                end
            end
            save(fullfile(output_folder,'cheCV_chk.mat'),'SDS_done_flag','each_photon_weight_arr','detpt_weight_arr','SDS_CV_arr');
        end
                
        %% checkpoint
        if do_checkpoint~=0
            % save this sub simulation
            to_save=[wl, run_index];
            save(fullfile(output_folder,'checkpioint_index.txt'),'to_save','-ascii','-tabs');
            save(fullfile(output_folder,'simtime_chk.txt'),'sim_time_calc','-ascii','-tabs');
            
            if output_fluence_rate
                save(fullfile(output_folder,'fluencerate_chk.mat'),'fluencerate_arr');
            end
            
            if output_jacobian || output_pathlength
                save(fullfile(output_folder,'numpt_chk.txt'),'SDS_detected_numpt','-ascii','-tabs');
            end
            if output_jacobian && output_pathlength
                save(fullfile(output_folder,['detpt_chk_' num2str(run_index) '.mat']),'seed_arr','detpt_arr','detPL_arr');
            elseif output_jacobian
                save(fullfile(output_folder,['detpt_chk_' num2str(run_index) '.mat']),'seed_arr','detpt_arr');
            elseif output_pathlength
                save(fullfile(output_folder,['detpt_chk_' num2str(run_index) '.mat']),'detPL_arr');
            end
            
            % clear the var for save RAM
            if output_jacobian
                seed_arr=cell(1,num_SDS);
                for d=1:num_SDS
                    detpt_arr{d}=[];
                end
            end
            if output_pathlength
                for d=1:num_SDS
                    detPL_arr{d}=[];
                end
            end
        end
        if rem(run_index,200)==0
            fprintf('\n');
        end
        
        %% if terminate manually
        if exist('stop_flag.txt','file')>0
            stop_flat=load('stop_flag.txt');
            if stop_flat==1
                assert(false,'User terminate the script!');
            end
        end
        
        %% if CV is small enough
        if max_CV~=0 && sum(SDS_done_flag==0)==0
            break;
        end
    end
    fprintf('\n');
    
    timer_sim=tic;
    
    if max_CV==0
        assert(run_index==num_need_to_run,'ERROR: run_index not correct!');
    end
    
    %% average the fluence rate
    if output_fluence_rate
        average_fluence_rate=fluencerate_arr./run_index;
        save(fullfile(output_folder,['average_fluence_' num2str(wl) '.mat']),'average_fluence_rate');
    end
    
    %% calculate PL
    if output_pathlength
        
        % if do checkpoint, then load the simulated seed array
        if do_checkpoint~=0
            % init
            SDS_detpt_arr=cell(1,num_SDS);
            for s=1:num_SDS
                SDS_detpt_arr{s}=single(zeros(SDS_detected_numpt(s),length(to_output_layer)));
            end
            % load the checkpoints
            SDS_load_numpt=ones(1,num_SDS);
            for ridx=1:run_index
                temp_detPL=load(fullfile(output_folder,['detpt_chk_' num2str(ridx) '.mat']));
                for s=1:num_SDS
                    old_SDS_numpt=SDS_load_numpt(s);
                    SDS_load_numpt(s)=SDS_load_numpt(s)+size(temp_detPL.detPL_arr{s},1);
                    SDS_detpt_arr{s}(old_SDS_numpt:SDS_load_numpt(s)-1,:)=temp_detPL.detPL_arr{s};
                end
            end
        else
            SDS_detpt_arr=detPL_arr;
        end
        
        sim_summary.SDS_detected_number=zeros(1,num_SDS);
        for s=1:num_SDS
            sim_summary.SDS_detected_number(1,s)=size(SDS_detpt_arr{s},1);
        end
        
        if max_CV>0
            for s=1:num_SDS
                if each_photon_weight_arr(s)==0
                    assert(SDS_done_flag(s)==0,'ERRROR: SDS had done but don''t know when!');
                    each_photon_weight_arr(s)=num_photon;
                end
            end
        else
            each_photon_weight_arr=ones(1,num_SDS)*num_photon;
        end
                
        %% save output file
        detPL_state=whos('SDS_detpt_arr');
        if detPL_state.bytes/1024^2>500 % if the file will be too large
            save(fullfile(output_folder,['PL_' num2str(wl) '.mat']),'SDS_detpt_arr','each_photon_weight_arr','-v7.3');
        else
            save(fullfile(output_folder,['PL_' num2str(wl) '.mat']),'SDS_detpt_arr','each_photon_weight_arr');
            if exist(fullfile(output_folder,['PL_' num2str(wl) '.mat']),'file')==0
                save(fullfile(output_folder,['PL_' num2str(wl) '.mat']),'SDS_detpt_arr','each_photon_weight_arr','-v7.3');
            end
        end
        
        if max_CV>0
            save(fullfile(output_folder,['SDS_CV_arr_' num2str(wl) '.txt']),'SDS_CV_arr','-ascii','-tabs');
            sim_summary.had_check_CV=1;
        end
        
        %% calculate reflectance
        for s=1:num_SDS
            if size(SDS_detpt_arr{s},1)==0
                reflectance_arr(wl,s)=0;
            else
                reflectance_arr(wl,s)=1/each_photon_weight_arr(s)*sum(exp(-double(sum(SDS_detpt_arr{s}.*tissue_param(wl,(2*to_output_layer)-1),2))));
            end
        end
        to_save=reflectance_arr;
        save(fullfile(output_folder,'reflectance.txt'),'to_save','-ascii','-tabs');
    end
    
    %% stop timer
    sim_summary.sim_time=sim_time_calc+toc(timer_sim);
    
    if output_jacobian
        %% start replay timer
        timer_replay=tic;
        
        %% if do checkpoint, then load the simulated seed array
        if do_checkpoint~=0
            % init
            seed_arr=cell(1,num_SDS);
            detpt_arr=cell(1,num_SDS);
            temp_detpt=load(fullfile(output_folder,['detpt_chk_1.mat']));
            for s=1:num_SDS
                seed_arr{s}=uint8(zeros(size(temp_detpt.seed_arr{1},1),SDS_detected_numpt(s)));
                detpt_arr{s}=single(zeros(size(temp_detpt.detpt_arr{1},1),SDS_detected_numpt(s)));
            end
            % load the checkpoints
            SDS_load_numpt=ones(1,num_SDS);
            for ridx=1:run_index
                temp_detpt=load(fullfile(output_folder,['detpt_chk_' num2str(ridx) '.mat']));
                for s=1:num_SDS
                    old_SDS_numpt=SDS_load_numpt(s);
                    SDS_load_numpt(s)=SDS_load_numpt(s)+size(temp_detpt.seed_arr{s},2);
                    seed_arr{s}(:,old_SDS_numpt:SDS_load_numpt(s)-1)=temp_detpt.seed_arr{s};
                    detpt_arr{s}(:,old_SDS_numpt:SDS_load_numpt(s)-1)=temp_detpt.detpt_arr{s};
                end
            end
        end
        
        %% init array to store the detected photon
        jacobian_arr=cell(1,num_SDS);
        
        %% replay the detected photon
        fprintf('\tReplaying SDS:');
        for d=1:num_SDS
            fprintf('\t%d,',d);
            if size(seed_arr{d},2)>0
                replay_run_index=0;
                replay_num_need_to_run=1;
                if size(seed_arr{d},2)>mcx_max_detpt
                    replay_num_need_to_run=ceil(size(seed_arr{d},2)/mcx_max_detpt);
                end
                while replay_run_index<replay_num_need_to_run
                    newcfg=cfg;
                    newcfg.gpuid=1;
                    newcfg.workload=100;
                    newcfg.outputtype='jacobian';
                    if replay_run_index<replay_num_need_to_run-1
                        newcfg.seed=seed_arr{d}(:,replay_run_index*mcx_max_detpt+1:(replay_run_index+1)*mcx_max_detpt);
                        newcfg.detphotons=detpt_arr{d}(:,replay_run_index*mcx_max_detpt+1:(replay_run_index+1)*mcx_max_detpt);
                    else
                        newcfg.seed=seed_arr{d}(:,replay_run_index*mcx_max_detpt+1:end);
                        newcfg.detphotons=detpt_arr{d}(:,replay_run_index*mcx_max_detpt+1:end);
                    end
                    if show_MCX_info
                        [new_fluencerate,new_detpt,~,~]=mcxlab(newcfg);
                    else
                        [T,new_fluencerate,new_detpt,~,~]=evalc('mcxlab(newcfg);');
                    end
                    if isfield(new_detpt,'detid')
                        assert(length(new_detpt.detid)<=newcfg.maxdetphoton,'Error: detected too much photon!'); % to avoid the loss due to too many detected photon
                        if replay_run_index==0
                            jacobian_arr{d}=new_fluencerate.data;
                        else
                            jacobian_arr{d}=jacobian_arr{d}+new_fluencerate.data;
                        end
                    end
                    replay_run_index=replay_run_index+1;
                end
                
            else % if there are no detected photon for this detector
                jacobian_arr{d}=0;
            end
        end
        fprintf('\n');
        
        %% save output file
        if output_jacobian
            save(fullfile(output_folder,['jacobian_' num2str(wl) '.mat']),'jacobian_arr');
        end
        
        sim_summary.replay_time=toc(timer_replay);
        sim_summary.sim_speed=sim_summary.num_photon/(sim_summary.sim_time+sim_summary.replay_time);
    else
        sim_summary.sim_speed=sim_summary.num_photon/sim_summary.sim_time;
    end
    
    %% output the summary file
    try
        fid=fopen(fullfile(output_folder,['sim_summary_' num2str(wl) '.json']),'w');
        fprintf(fid,jsonencode(sim_summary));
        fclose(fid);
    end
    
    %% edit and delete the checkpoint file
    to_save=[wl+1, 0];
    save(fullfile(output_folder,'checkpioint_index.txt'),'to_save','-ascii','-tabs');
    
    if exist(fullfile(output_folder,'fluencerate_chk.mat'),'file')~=0
        delete(fullfile(output_folder,'fluencerate_chk.mat'));
    end
    
    if exist(fullfile(output_folder,'numpt_chk.txt'),'file')~=0
        delete(fullfile(output_folder,'numpt_chk.txt'));
    end
    
    if exist(fullfile(output_folder,'simtime_chk.txt'),'file')~=0
        delete(fullfile(output_folder,'simtime_chk.txt'));
    end
    
    for chkidx=1:num_need_to_run
        if exist(fullfile(output_folder,['detpt_chk_' num2str(chkidx) '.mat']),'file')~=0
            delete(fullfile(output_folder,['detpt_chk_' num2str(chkidx) '.mat']));
        end 
    end
    from_chk_flag=0;
    
    if exist(fullfile(output_folder,'cheCV_chk.mat'),'file')~=0
        delete(fullfile(output_folder,'cheCV_chk.mat'));
    end
end

%% delete the checkpoint file
if exist(fullfile(output_folder,'checkpioint_index.txt'),'file')~=0
    delete(fullfile(output_folder,'checkpioint_index.txt'));
end
end
