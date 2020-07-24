%Copyright © 2018, Sampsa Pursiainen
function [void] = print_meshes(void);

void = [];

length_reconstruction_cell = 1;

aux_wm_ind = -1;


if ismember(evalin('base','zef.on_screen'), [0,1]) && not(evalin('base','zef.visualization_type')==3) 
    
h_axes_text = [];    
    
number_of_frames = evalin('base','zef.number_of_frames');
file_index = evalin('base','zef.file_index');
file_name = evalin('base','zef.file');
file_path = evalin('base','zef.file_path');
snapshot_resolution = round([evalin('base','zef.snapshot_vertical_resolution') evalin('base','zef.snapshot_horizontal_resolution')]);
    
c_va = camva(evalin('base','zef.h_axes1'));
c_pos = campos(evalin('base','zef.h_axes1'));
c_ta = camtarget(evalin('base','zef.h_axes1'));
c_p = camproj(evalin('base','zef.h_axes1')); 
c_u = camup(evalin('base','zef.h_axes1'));

h_fig_aux = figure;set(h_fig_aux,'visible','on'); 
clf;
set(h_fig_aux,'renderer','opengl');
set(h_fig_aux,'paperunits','inches');
set(h_fig_aux,'papersize',snapshot_resolution);
set(h_fig_aux,'paperposition',[0 0 fliplr(snapshot_resolution)]);
light('Position',[0 0 1],'Style','infinite');
light('Position',[0 0 -1],'Style','infinite');
h_axes_image = get(h_fig_aux,'currentaxes');
hold on;

cp_a = evalin('base','zef.cp_a');
cp_b = evalin('base','zef.cp_b');
cp_c = evalin('base','zef.cp_c');
cp_d = evalin('base','zef.cp_d');
cp2_a = evalin('base','zef.cp2_a');
cp2_b = evalin('base','zef.cp2_b');
cp2_c = evalin('base','zef.cp2_c');
cp2_d = evalin('base','zef.cp2_d');
cp3_a = evalin('base','zef.cp3_a');
cp3_b = evalin('base','zef.cp3_b');
cp3_c = evalin('base','zef.cp3_c');
cp3_d = evalin('base','zef.cp3_d');

if iscell(evalin('base','zef.reconstruction'))
length_reconstruction_cell = length(evalin('base','zef.reconstruction'));
else
length_reconstruction_cell = 1;
end
frame_start = evalin('base','zef.frame_start');
frame_stop = evalin('base','zef.frame_stop');
frame_step = evalin('base','zef.frame_step');
frame_step = max(frame_step,1);
if frame_start == 0
    frame_start = 1;
end
if frame_stop == 0
    frame_stop = length_reconstruction_cell;
end
frame_start = max(frame_start,1);
frame_start = min(length_reconstruction_cell,frame_start);
frame_stop = max(frame_stop,1);
frame_stop = min(length_reconstruction_cell,frame_stop);
if length([frame_start : frame_step : frame_stop]) > 1 
is_video = 1;
else 
is_video = 0;
end

if iscell(evalin('base','zef.reconstruction')) &  evalin('base','zef.visualization_type') == 2 & is_video & file_index == 4
avi_file_temp = [file_path file_name(1:end-4) '_temp.avi'];
avi_file = [file_path file_name];
video_quality = str2num(evalin('base','zef.video_codec'));
h_aviobj = VideoWriter(avi_file_temp);
h_aviobj.Quality = video_quality; 
open(h_aviobj);
end

%T
sensors = evalin('base','zef.sensors');
surface_triangles = evalin('base','zef.surface_triangles');
nodes = evalin('base','zef.nodes');

if size(sensors,2) == 6 & ismember(evalin('base','zef.imaging_method'), [1 4 5])
    electrode_model = 2;
else
    electrode_model = 1;
end

[X_s, Y_s, Z_s] = sphere(20);  
aux_scale_ind = evalin('base','zef.location_unit');
switch aux_scale_ind 
    case 1
        aux_scale_val = 1;
    case 2 
        aux_scale_val = 0.01;
    case 3 
        aux_scale_val = 0.001;
end


aux_ind = []; 
if evalin('base','zef.s_visible')
sphere_scale = 3.2*aux_scale_val;    
X_s = sphere_scale*X_s;
Y_s = sphere_scale*Y_s;
Z_s = sphere_scale*Z_s;
if evalin('base','zef.cp_on');
if not(isempty(aux_ind))
aux_ind = intersect(aux_ind,find(sum(sensors(:,1:3).*repmat([cp_a cp_b cp_c],size(sensors,1),1),2) >= cp_d));
else
aux_ind = [find(sum(sensors(:,1:3).*repmat([cp_a cp_b cp_c],size(sensors,1),1),2) >= cp_d)];
end
end
if evalin('base','zef.cp2_on');
if not(isempty(aux_ind))
aux_ind = intersect(aux_ind,find(sum(sensors(:,1:3).*repmat([cp2_a cp2_b cp2_c],size(sensors,1),1),2) >= cp2_d));
else
aux_ind = [find(sum(sensors(:,1:3).*repmat([cp2_a cp2_b cp2_c],size(sensors,1),1),2) >= cp2_d)];
end
end
if evalin('base','zef.cp3_on');
if not(isempty(aux_ind))
aux_ind = intersect(aux_ind,find(sum(sensors(:,1:3).*repmat([cp3_a cp3_b cp3_c],size(sensors,1),1),2) >= cp3_d));
else
aux_ind = [find(sum(sensors(:,1:3).*repmat([cp3_a cp3_b cp3_c],size(sensors,1),1),2) >= cp3_d)];
end
end
if evalin('base','zef.cp_on') || evalin('base','zef.cp2_on') || evalin('base','zef.cp3_on')
if evalin('base','zef.cp_mode') == 1
sensors = sensors(aux_ind,:);
elseif evalin('base','zef.cp_mode') == 2
aux_ind = setdiff([1:size(sensors,1)]',aux_ind);
sensors = sensors(aux_ind,:);   
end
end
aux_ind = [];

if electrode_model == 1 & evalin('base','zef.attach_electrodes') & ismember(evalin('base','zef.imaging_method'),[1 4 5]) 
sensors = attach_sensors_volume(sensors); 
elseif electrode_model==2  & ismember(evalin('base','zef.imaging_method'),[1 4 5]) 
sensors = attach_sensors_volume(sensors);
end

if electrode_model == 1 | not(ismember(evalin('base','zef.imaging_method'),[1,4,5]))
for i = 1 : size(sensors,1)
h = surf(sensors(i,1) + X_s, sensors(i,2) + Y_s, sensors(i,3) + Z_s);
set(h,'facecolor',evalin('base','zef.s_color'));
set(h,'edgecolor','none');
set(h,'specularstrength',0.3);
set(h,'diffusestrength',0.7);
set(h,'ambientstrength',0.7);
set(h,'facealpha',evalin('base','zef.layer_transparency'));
end
else
h = trisurf(sensors(:,2:4),nodes(:,1),nodes(:,2),nodes(:,3));
set(h,'facecolor',evalin('base','zef.s_color'));
set(h,'edgecolor',evalin('base','zef.s_color')); 
set(h,'specularstrength',0.3);
set(h,'diffusestrength',0.7);
set(h,'ambientstrength',0.7);
set(h,'facealpha',evalin('base','zef.layer_transparency'));
set(h,'edgealpha',evalin('base','zef.layer_transparency'));    
end
if ismember(evalin('base','zef.imaging_method'),[2,3])
sensors(:,4:6) = sensors(:,4:6)./repmat(sqrt(sum(sensors(:,4:6).^2,2)),1,3);
h=coneplot(sensors(:,1) + 4.5*sensors(:,4),sensors(:,2) + 4.5*sensors(:,5),sensors(:,3) + 4.5*sensors(:,6),8*sensors(:,4),8*sensors(:,5),8*sensors(:,6),0,'nointerp');
set(h,'facecolor',evalin('base','zef.s_color'));
set(h,'edgecolor','none'); 
set(h,'specularstrength',0.3);
set(h,'diffusestrength',0.7);
set(h,'ambientstrength',0.7);
set(h,'facealpha',evalin('base','zef.layer_transparency'));
if size(sensors,2) == 9
sensors(:,7:9) = sensors(:,7:9)./repmat(sqrt(sum(sensors(:,7:9).^2,2)),1,3);
h=coneplot(sensors(:,1) + 4.5*sensors(:,7),sensors(:,2) + 4.5*sensors(:,8),sensors(:,3) + 4.5*sensors(:,9),8*sensors(:,7),8*sensors(:,8),8*sensors(:,9),0,'nointerp');
set(h,'facecolor',0.9*[0 1 1]);
set(h,'edgecolor','none'); 
set(h,'specularstrength',0.3);
set(h,'diffusestrength',0.7);
set(h,'ambientstrength',0.7);
set(h,'facealpha',evalin('base','zef.layer_transparency'));
end
end
end

i = 0;
length_reuna = 0;
sigma_vec = [];
priority_vec = [];
visible_vec = [];
color_cell = cell(0);
aux_brain_ind = [];
for k = 1 : 18   
switch k
    case 1
        var_0 = 'zef.d1_on';
        var_1 = 'zef.d1_sigma';
        var_2 = 'zef.d1_priority';
        var_3 = 'zef.d1_visible';
    color_str = evalin('base','zef.d1_color');
     case 2
        var_0 = 'zef.d2_on';
        var_1 = 'zef.d2_sigma';   
        var_2 = 'zef.d2_priority';
        var_3 = 'zef.d2_visible';
        color_str = evalin('base','zef.d2_color');
     case 3
        var_0 = 'zef.d3_on';
        var_1 = 'zef.d3_sigma';   
        var_2 = 'zef.d3_priority';
        var_3 = 'zef.d3_visible';
        color_str = evalin('base','zef.d3_color');
     case 4
        var_0 = 'zef.d4_on';
        var_1 = 'zef.d4_sigma';   
        var_2 = 'zef.d4_priority';
        var_3 = 'zef.d4_visible';
        color_str = evalin('base','zef.d4_color');
     case 5
        var_0 = 'zef.d5_on';
        var_1 = 'zef.d5_sigma';
        var_2 = 'zef.d5_priority';
        var_3 = 'zef.d5_visible';
        color_str = evalin('base','zef.d5_color');
     case 6
        var_0 = 'zef.d6_on';
        var_1 = 'zef.d6_sigma';   
        var_2 = 'zef.d6_priority';
        var_3 = 'zef.d6_visible';
        color_str = evalin('base','zef.d6_color');
     case 7
        var_0 = 'zef.d7_on';
        var_1 = 'zef.d7_sigma';   
        var_2 = 'zef.d7_priority';
        var_3 = 'zef.d7_visible';
        color_str = evalin('base','zef.d7_color');
     case 8
        var_0 = 'zef.d8_on';
        var_1 = 'zef.d8_sigma';   
        var_2 = 'zef.d8_priority';
        var_3 = 'zef.d8_visible';
        color_str = evalin('base','zef.d8_color');
     case 9
        var_0 = 'zef.d9_on';
        var_1 = 'zef.d9_sigma';
        var_2 = 'zef.d9_priority';
        var_3 = 'zef.d9_visible';
       color_str = evalin('base','zef.d9_color');
     case 10
        var_0 = 'zef.d10_on';
        var_1 = 'zef.d10_sigma';   
        var_2 = 'zef.d10_priority';
        var_3 = 'zef.d10_visible';
        color_str = evalin('base','zef.d10_color');
     case 11
        var_0 = 'zef.d11_on';
        var_1 = 'zef.d11_sigma';   
        var_2 = 'zef.d11_priority';
        var_3 = 'zef.d11_visible';
        color_str = evalin('base','zef.d11_color');
     case 12
        var_0 = 'zef.d12_on';
        var_1 = 'zef.d12_sigma';   
        var_2 = 'zef.d12_priority';
        var_3 = 'zef.d12_visible';
        color_str = evalin('base','zef.d12_color');
    case 13
        var_0 = 'zef.d13_on';
        var_1 = 'zef.d13_sigma';   
        var_2 = 'zef.d13_priority';
        var_3 = 'zef.d13_visible';
        color_str = evalin('base','zef.d13_color');
    case 14
        var_0 = 'zef.w_on';
        var_1 = 'zef.w_sigma';    
        var_2 = 'zef.w_priority';
        var_3 = 'zef.w_visible';
        color_str = evalin('base','zef.w_color');
    case 15
        var_0 = 'zef.g_on';
        var_1 = 'zef.g_sigma';
        var_2 = 'zef.g_priority';
        var_3 = 'zef.g_visible';
        color_str = evalin('base','zef.g_color');
    case 16
        var_0 = 'zef.c_on';
        var_1 = 'zef.c_sigma';
        var_2 = 'zef.c_priority';
        var_3 = 'zef.c_visible';
        color_str = evalin('base','zef.c_color');
     case 17
        var_0 = 'zef.sk_on';
        var_1 = 'zef.sk_sigma';
        var_2 = 'zef.sk_priority';
        var_3 = 'zef.sk_visible';
        color_str = evalin('base','zef.sk_color');
     case 18
        var_0 = 'zef.sc_on';
        var_1 = 'zef.sc_sigma';
        var_2 = 'zef.sc_priority';
        var_3 = 'zef.sc_visible';
        color_str = evalin('base','zef.sc_color');
     end
on_val = evalin('base',var_0);      
sigma_val = evalin('base',var_1);  
priority_val = evalin('base',var_2);
visible_val = evalin('base',var_3);
if on_val
i = i + 1;
sigma_vec(i,1) = sigma_val;
priority_vec(i,1) = priority_val;
color_cell{i} = color_str;
visible_vec(i,1) = i*visible_val;

if k == 1 && evalin('base','zef.d1_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 2 && evalin('base','zef.d2_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 3 && evalin('base','zef.d3_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 4 && evalin('base','zef.d4_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 5 && evalin('base','zef.d5_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 6 && evalin('base','zef.d6_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 7 && evalin('base','zef.d7_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 8 && evalin('base','zef.d8_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 9 && evalin('base','zef.d9_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 10 && evalin('base','zef.d10_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 11 && evalin('base','zef.d11_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 12 && evalin('base','zef.d12_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 13 && evalin('base','zef.d13_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 14 && evalin('base','zef.wm_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 15 && evalin('base','zef.g_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 16 && evalin('base','zef.c_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 17 && evalin('base','zef.sk_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 18 && evalin('base','zef.sc_sources');
    aux_brain_ind = [aux_brain_ind i];
end
end
end

johtavuus = evalin('base','zef.sigma');
johtavuus = johtavuus(:,2);
 
I = find(ismember(johtavuus,visible_vec));
johtavuus = johtavuus(I);
tetra = evalin('base','zef.tetra');

tetra = tetra(I,:);
tetra_c = (1/4)*(nodes(tetra(:,1),:) + nodes(tetra(:,2),:) + nodes(tetra(:,3),:) + nodes(tetra(:,4),:));

aux_ind = [];
if evalin('base','zef.cp_on');

aux_ind = [find(sum(tetra_c.*repmat([cp_a cp_b cp_c],size(tetra_c,1),1),2) >= cp_d)];
end
if evalin('base','zef.cp2_on');

if not(isempty(aux_ind))
aux_ind = intersect(aux_ind,find(sum(tetra_c.*repmat([cp2_a cp2_b cp2_c],size(tetra_c,1),1),2) >= cp2_d));
else
aux_ind = [find(sum(tetra_c.*repmat([cp2_a cp2_b cp2_c],size(tetra_c,1),1),2) >= cp2_d)];
end
end
if evalin('base','zef.cp3_on');

if not(isempty(aux_ind))
aux_ind = intersect(aux_ind,find(sum(tetra_c.*repmat([cp3_a cp3_b cp3_c],size(tetra_c,1),1),2) >= cp3_d));
else
aux_ind = [find(sum(tetra_c.*repmat([cp3_a cp3_b cp3_c],size(tetra_c,1),1),2) >= cp3_d)];
end
end

if evalin('base','zef.cp_on') || evalin('base','zef.cp2_on') || evalin('base','zef.cp3_on')
if evalin('base','zef.cp_mode') == 1
tetra = tetra(aux_ind,:);
elseif evalin('base','zef.cp_mode') == 2
aux_ind = setdiff([1:size(tetra,1)]',aux_ind);
tetra = tetra(aux_ind,:);   
elseif evalin('base','zef.cp_mode') == 3
aux_ind = union(aux_ind,find(johtavuus==aux_brain_ind));
tetra = tetra(aux_ind,:);  
elseif evalin('base','zef.cp_mode') == 4
aux_ind = setdiff([1:size(tetra,1)]',aux_ind);
aux_ind = union(aux_ind,find(johtavuus==aux_brain_ind));
tetra = tetra(aux_ind,:);  
end
else
aux_ind = [1:size(tetra,1)]';
end;
I_aux = I(aux_ind);

 ind_m = [ 2 3 4;
           1 4 3;
           1 2 4; 
           1 3 2]; 

tetra_sort = [tetra(:,[2 3 4]) ones(size(tetra,1),1) [1:size(tetra,1)]'; 
              tetra(:,[1 4 3]) 2*ones(size(tetra,1),1) [1:size(tetra,1)]'; 
              tetra(:,[1 2 4]) 3*ones(size(tetra,1),1) [1:size(tetra,1)]'; 
              tetra(:,[1 3 2]) 4*ones(size(tetra,1),1) [1:size(tetra,1)]';];
tetra_sort(:,1:3) = sort(tetra_sort(:,1:3),2);

if evalin('base','zef.use_gpu') == 1 & gpuDeviceCount > 0
tetra_sort = gpuArray(uint32(tetra_sort));
tetra_sort = gather(sortrows(tetra_sort,[1 2 3])); 
else
tetra_sort = sortrows(tetra_sort,[1 2 3]);   
end
tetra_ind = zeros(size(tetra_sort,1),1);
I = find(sum(abs(tetra_sort(2:end,1:3)-tetra_sort(1:end-1,1:3)),2)==0);
tetra_ind(I) = 1;
tetra_ind(I+1) = 1;
I = find(tetra_ind == 0);
tetra_ind = sub2ind(size(tetra),repmat(tetra_sort(I,5),1,3),ind_m(tetra_sort(I,4),:));
surface_triangles = tetra(tetra_ind);
tetra_ind = tetra_sort(I,5);
clear tetra_sort;

n_compartments = i;

max_abs_reconstruction = -Inf;
min_rec = Inf;
max_rec = -Inf;
frame_start = 1;
frame_stop = 1;
frame_step = 1;

if ismember(evalin('base','zef.visualization_type'), [2,4])
s_i_ind = evalin('base','zef.source_interpolation_ind{1}');
end

if evalin('base','zef.use_parcellation')
selected_list = evalin('base','zef.parcellation_selected');
p_i_ind = evalin('base','zef.parcellation_interp_ind');
end

if iscell(evalin('base','zef.reconstruction')) &&  evalin('base','zef.visualization_type') == 2
length_reconstruction_cell = length(evalin('base','zef.reconstruction'));
frame_start = evalin('base','zef.frame_start');
frame_stop = evalin('base','zef.frame_stop');
frame_step = evalin('base','zef.frame_step');
if frame_start == 0
    frame_start = 1;
end
if frame_stop == 0
    frame_stop = length_reconstruction_cell;
end
frame_start = max(frame_start,1);
frame_start = min(length_reconstruction_cell,frame_start);
frame_stop = max(frame_stop,1);
frame_stop = min(length_reconstruction_cell,frame_stop);
number_of_frames = length([frame_start : frame_step : frame_stop]);
for f_ind = frame_start : frame_step : frame_stop
reconstruction = evalin('base',['zef.reconstruction{' int2str(f_ind) '}']);
reconstruction = reconstruction(:);  
reconstruction = reshape(reconstruction,3,length(reconstruction)/3);
if ismember(evalin('base','zef.reconstruction_type'), 6)
reconstruction = (1/sqrt(3))*sum(reconstruction)';
else    
reconstruction = sqrt(sum(reconstruction.^2))';
end
reconstruction = sum(reconstruction(s_i_ind),2)/4;
max_abs_reconstruction = max([max_abs_reconstruction ; (reconstruction(:))]);
min_rec = min([min_rec ; (reconstruction(:))]);
max_rec = max_abs_reconstruction;
end
if not(ismember(evalin('base','zef.reconstruction_type'), [6]))
if evalin('base','zef.inv_scale') == 1
min_rec = 10*log10(max(min_rec/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
max_rec = 10*log10(max(max_rec/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
elseif evalin('base','zef.inv_scale') == 2
min_rec = (max(min_rec/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
max_rec = (max(max_rec/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
elseif evalin('base','zef.inv_scale') == 3
min_rec = sqrt(max(min_rec/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
max_rec = sqrt(max(max_rec/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
end
end
elseif  evalin('base','zef.visualization_type') == 2
%s_i_ind = evalin('base','zef.source_interpolation_ind{1}');
reconstruction = evalin('base','zef.reconstruction');
reconstruction = reconstruction(:);  
reconstruction = reshape(reconstruction,3,length(reconstruction)/3);
if ismember(evalin('base','zef.reconstruction_type'), 6)
reconstruction = (1/sqrt(3))*sum(reconstruction)';
else    
reconstruction = sqrt(sum(reconstruction.^2))';
end
reconstruction = sum(reconstruction(s_i_ind),2)/4;
max_abs_reconstruction = max([max_abs_reconstruction ; (reconstruction(:))]);
min_rec = min([min_rec ; (reconstruction(:))]);
max_rec = max_abs_reconstruction;
if not(ismember(evalin('base','zef.reconstruction_type'), [6]))
if evalin('base','zef.inv_scale') == 1
min_rec = 10*log10(max(min_rec/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
max_rec = 10*log10(max(max_rec/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
elseif evalin('base','zef.inv_scale') == 2
min_rec = (max(min_rec/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
max_rec = (max(max_rec/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
elseif evalin('base','zef.inv_scale') == 3
min_rec = sqrt(max(min_rec/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
max_rec = sqrt(max(max_rec/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
end
end
end

if ismember(evalin('base','zef.visualization_type'),[2])
if  iscell(evalin('base','zef.reconstruction')) & evalin('base','zef.visualization_type') == 2
h_waitbar = waitbar(1/number_of_frames,['Frame ' int2str(1) ' of ' int2str(number_of_frames) '.']);    
set(h_waitbar,'handlevisibility','off');
end
end

f_ind_aux = 1;
for f_ind = frame_start : frame_start
if  iscell(evalin('base','zef.reconstruction')) & evalin('base','zef.visualization_type') == 2    
waitbar(f_ind_aux/number_of_frames,h_waitbar,['Frame ' int2str(f_ind_aux) ' of ' int2str(number_of_frames) '.']);    
set(h_waitbar,'handlevisibility','off');
end
%axes(evalin('base','zef.h_axes1'));
%if not(isempty(h_colorbar))
%colorbar(h_colorbar,'delete'); 
%h_colorbar = [];
%end
hold on;
%**************************************************************************
if ismember(evalin('base','zef.visualization_type'),[2,4])
brain_ind_aux = evalin('base','zef.brain_ind');
brain_ind = brain_ind_aux;
[aux_vec, brain_ind, I_2] = intersect(I_aux,brain_ind);
clear aux_vec;
johtavuus(aux_ind(brain_ind))=0;
I_3 = find(ismember(tetra_ind,brain_ind));

if evalin('base','zef.use_parcellation')
reconstruction_p_1 = ones(size(I_3,1),1);
reconstruction_p_2 = zeros(size(I_3,1),1);
for p_ind = selected_list
p_ind_aux = brain_ind_aux(p_i_ind{p_ind}{1});
[p_ind_aux,p_ind_aux_1,p_ind_aux_2] = intersect(I_aux, p_ind_aux);
[p_ind_aux] = find(ismember(tetra_ind(I_3),p_ind_aux_1));
reconstruction_p_1(p_ind_aux) = p_ind+1;
reconstruction_p_2(p_ind_aux) = 1;
end
end
end

if ismember(evalin('base','zef.visualization_type'),[4])
reconstruction = reconstruction_p_1;
min_rec = 1; 
max_rec = size(evalin('base','zef.parcellation_colormap'),1);
end

if ismember(evalin('base','zef.visualization_type'),[2])

%******************************************************
if iscell(evalin('base','zef.reconstruction')) 
reconstruction = evalin('base',['zef.reconstruction{' int2str(frame_start) '}']);
else
reconstruction = evalin('base','zef.reconstruction');
end
reconstruction = reconstruction(:);  
reconstruction = reshape(reconstruction,3,length(reconstruction)/3);
reconstruction = sqrt(sum(reconstruction.^2))';
reconstruction = sum(reconstruction(s_i_ind),2)/4;
if evalin('base','zef.inv_scale') == 1
reconstruction = 10*log10(max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
elseif evalin('base','zef.inv_scale') == 2
reconstruction = (max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
elseif evalin('base','zef.inv_scale') == 3
reconstruction = sqrt(max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
end
[a_hist, b_hist] = hist(reconstruction,min_rec:(max_rec-min_rec)/50:max_rec - (max_rec - min_rec)/50);
a_hist = max(0,real(log10(a_hist)));
length_reconstruction = length(reconstruction);
%******************************************************


%******************************************************
if iscell(evalin('base','zef.reconstruction')) 
reconstruction = evalin('base',['zef.reconstruction{' int2str(frame_start) '}']);
else
reconstruction = evalin('base','zef.reconstruction');
end
reconstruction = reconstruction(:);  
reconstruction = reshape(reconstruction,3,length(reconstruction)/3);
reconstruction = sqrt(sum(reconstruction.^2))';
reconstruction = sum(reconstruction(s_i_ind),2)/4;
if evalin('base','zef.inv_scale') == 1
reconstruction = 10*log10(max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
elseif evalin('base','zef.inv_scale') == 2
reconstruction = (max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
elseif evalin('base','zef.inv_scale') == 3
reconstruction = sqrt(max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
end
[a_hist, b_hist] = hist(reconstruction,min_rec:(max_rec-min_rec)/50:max_rec - (max_rec - min_rec)/50);
a_hist = max(0,real(log10(a_hist)));
length_reconstruction = length(reconstruction);
%******************************************************

if iscell(evalin('base','zef.reconstruction'))
reconstruction = evalin('base',['zef.reconstruction{' int2str(f_ind) '}']);
else
reconstruction = evalin('base','zef.reconstruction');  
end
reconstruction = reconstruction(:);  
reconstruction = reshape(reconstruction,3,length(reconstruction)/3);

if ismember(evalin('base','zef.reconstruction_type'),[1 7])
reconstruction = sqrt(sum(reconstruction.^2))';
elseif evalin('base','zef.reconstruction_type') == 6
reconstruction = (1/sqrt(3))*sum(reconstruction)';
end
if ismember(evalin('base','zef.reconstruction_type'), [1 6 7])
reconstruction = sum(reconstruction(s_i_ind),2)/4;
reconstruction = reconstruction(I_2);
I_2_b_rec = I_2;
I_3_rec = I_3;
I_2 = zeros(length(aux_ind),1);
I_2(brain_ind) = [1:length(brain_ind)]';
I_2_rec = I_2;
I_1 = tetra_ind(I_3);
I_1_rec = I_1;
reconstruction = reconstruction(I_2(I_1));
end

if ismember(evalin('base','zef.reconstruction_type'), [2 3 4 5])
rec_x = reconstruction(1,:)';
rec_y = reconstruction(2,:)';
rec_z = reconstruction(3,:)';
rec_x = sum(rec_x(s_i_ind),2)/4;
rec_y = sum(rec_y(s_i_ind),2)/4;
rec_z = sum(rec_z(s_i_ind),2)/4;
rec_x = rec_x(I_2);
rec_y = rec_y(I_2);
rec_z = rec_z(I_2);
I_2_b_rec = I_2;
I_3_rec = I_3;
I_2 = zeros(length(aux_ind),1);
I_2(brain_ind) = [1:length(brain_ind)]';
I_2_rec = I_2;
I_1 = tetra_ind(I_3);
I_1_rec = I_1;
rec_x = rec_x(I_2(I_1));
rec_y = rec_y(I_2(I_1));
rec_z = rec_z(I_2(I_1));
n_vec_aux = cross(nodes(surface_triangles(I_3,2),:)' - nodes(surface_triangles(I_3,1),:)',...
nodes(surface_triangles(I_3,3),:)' - nodes(surface_triangles(I_3,1),:)')';
n_vec_aux = n_vec_aux./repmat(sqrt(sum(n_vec_aux.^2,2)),1,3);
end 

if ismember(evalin('base','zef.reconstruction_type'), [2 3 4 5])
reconstruction = sqrt((rec_x.*n_vec_aux(:,1)).^2 + (rec_y.*n_vec_aux(:,2)).^2 + (rec_z.*n_vec_aux(:,3)).^2);
end

if evalin('base','zef.reconstruction_type') == 3
reconstruction = sqrt((rec_x - rec_x.*abs(n_vec_aux(:,1))).^2 + (rec_y - rec_y.*abs(n_vec_aux(:,2))).^2 + (rec_z - rec_z.*abs(n_vec_aux(:,3))).^2);
end

if evalin('base','zef.reconstruction_type') == 4
aux_rec = rec_x.*n_vec_aux(:,1) + rec_y.*n_vec_aux(:,2) + rec_z.*n_vec_aux(:,3);
I_aux_rec = find(aux_rec > 0);
reconstruction(I_aux_rec) = 0;
%reconstruction = reconstruction./max(abs(reconstruction(:)));
end

if evalin('base','zef.reconstruction_type') == 5 
aux_rec = rec_x.*n_vec_aux(:,1) + rec_y.*n_vec_aux(:,2) + rec_z.*n_vec_aux(:,3);
I_aux_rec = find(aux_rec <= 0);
reconstruction(I_aux_rec) = 0;
%reconstruction = reconstruction./max(abs(reconstruction(:)));
end

if ismember(evalin('base','zef.reconstruction_type'), [2 3 4 5 7])
reconstruction = smooth_field(surface_triangles(I_3,:), reconstruction, size(nodes,1),3);
end

if not(ismember(evalin('base','zef.reconstruction_type'),[6]))
if evalin('base','zef.inv_scale') == 1
reconstruction = 10*log10(max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
elseif evalin('base','zef.inv_scale') == 2
reconstruction = (max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
elseif evalin('base','zef.inv_scale') == 3
reconstruction = sqrt(max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
end
end
end

if ismember(evalin('base','zef.visualization_type'),[2,4])

if evalin('base','zef.use_parcellation')
reconstruction = reconstruction.*reconstruction_p_2;
end

colormap_size = 4096;
if evalin('base','zef.inv_colormap') == 1
c_aux_1 = floor(colormap_size/3);
c_aux_2 = floor(2*colormap_size/3);
colormap_vec_aux = [([20*[c_aux_1:-1:1] zeros(1,colormap_size-c_aux_1)]); ([15*[1: c_aux_1] 15*[c_aux_2-c_aux_1:-1:1] zeros(1,colormap_size-c_aux_2)]) ; ([zeros(1,c_aux_1) 6*[1:c_aux_2-c_aux_1] 6*[colormap_size-c_aux_2:-1:1]]);([zeros(1,c_aux_2) 7.5*[1:colormap_size-c_aux_2]])];
colormap_vec = zeros(3,size(colormap_vec_aux,2));
colormap_vec = colormap_vec + 0.52*[50*colormap_vec_aux(1,:) ; 50*colormap_vec_aux(1,:) ; 50*colormap_vec_aux(1,:)];
colormap_vec = colormap_vec + 0.5*[85*colormap_vec_aux(3,:) ; 197*colormap_vec_aux(3,:) ; 217*colormap_vec_aux(3,:)];
colormap_vec = colormap_vec + 0.1*[2*colormap_vec_aux(2,:) ; 118*colormap_vec_aux(2,:) ; 132*colormap_vec_aux(2,:)];
colormap_vec = colormap_vec + [203*colormap_vec_aux(4,:) ; 203*colormap_vec_aux(4,:) ; 100*colormap_vec_aux(4,:)];
clear colormap_vec_aux;
colormap_vec = colormap_vec'/max(colormap_vec(:));
colormap_vec = colormap_vec(:,1:3);
set(h_fig_aux,'colormap',colormap_vec);
elseif evalin('base','zef.inv_colormap') == 2
c_aux_1 = floor(colormap_size/3);
c_aux_2 = floor(2*colormap_size/3);
colormap_vec = zeros(3,colormap_size);
colormap_vec(1,:) =10*([colormap_size:-1:1]/colormap_size);
colormap_vec(2,:) = [10*( (3/2)*[c_aux_2:-1:1]/colormap_size) zeros(1,colormap_size-c_aux_2)];
colormap_vec(3,:) = [10*((3)*[c_aux_1:-1:1]/colormap_size) zeros(1,colormap_size-c_aux_1)];  
colormap_vec = colormap_vec'/max(colormap_vec(:));
colormap_vec = flipud(colormap_vec);
colormap_vec = colormap_vec(:,1:3);
colormap_vec = colormap_vec + repmat(0.2*([colormap_size:-1:1]'/colormap_size),1,3);
colormap_vec = colormap_vec/max(colormap_vec(:));
set(h_fig_aux,'colormap',colormap_vec);
elseif evalin('base','zef.inv_colormap') == 3
c_aux_1 = floor(colormap_size/3);
c_aux_2 = floor(2*colormap_size/3);
colormap_vec = zeros(3,colormap_size);
colormap_vec(2,:) = 10*([colormap_size:-1:1]/colormap_size);
colormap_vec(1,:) = [10*(3/2)*[c_aux_2:-1:1]/colormap_size zeros(1,colormap_size-c_aux_2)];
colormap_vec(3,:) = [10*(3*[c_aux_1:-1:1]/colormap_size) zeros(1,colormap_size-c_aux_1)];
%colormap_vec([1 3 2],:) = 0.75*colormap_vec([1 3 2],:) + 0.25*colormap_vec([1 2 3],:);
colormap_vec = colormap_vec'/max(colormap_vec(:));
colormap_vec = flipud(colormap_vec);
colormap_vec = colormap_vec(:,1:3);
colormap_vec = colormap_vec + repmat(0.2*([colormap_size:-1:1]'/colormap_size),1,3);
colormap_vec = colormap_vec/max(colormap_vec(:));
set(h_fig_aux,'colormap',colormap_vec);
elseif evalin('base','zef.inv_colormap') == 4
c_aux_1 = floor(colormap_size/3);
c_aux_2 = floor(2*colormap_size/3);
colormap_vec = zeros(3,colormap_size);
colormap_vec(3,:) = 10*([colormap_size:-1:1]/colormap_size);
colormap_vec(2,:) = [10*(3/2)*[c_aux_2:-1:1]/colormap_size zeros(1,colormap_size-c_aux_2)];
colormap_vec(1,:) = [10*(3*[c_aux_1:-1:1]/colormap_size) zeros(1,colormap_size-c_aux_1)];
%colormap_vec([1 3 2],:) = 0.75*colormap_vec([1 3 2],:) + 0.25*colormap_vec([1 2 3],:);
colormap_vec = colormap_vec'/max(colormap_vec(:));
colormap_vec = flipud(colormap_vec);
colormap_vec = colormap_vec(:,1:3);
colormap_vec = colormap_vec + repmat(0.2*([colormap_size:-1:1]'/colormap_size),1,3);
colormap_vec = colormap_vec/max(colormap_vec(:));
set(h_fig_aux,'colormap',colormap_vec);
elseif evalin('base','zef.inv_colormap') == 5
c_aux_1 = floor(colormap_size/3);
c_aux_2 = floor(2*colormap_size/3);
colormap_vec = [([20*[c_aux_1:-1:1] zeros(1,colormap_size-c_aux_1)]); ([15*[1: c_aux_1] 15*[c_aux_2-c_aux_1:-1:1] zeros(1,colormap_size-c_aux_2)]) ; ([zeros(1,c_aux_1) 6*[1:c_aux_2-c_aux_1] 6*[colormap_size-c_aux_2:-1:1]]);([zeros(1,c_aux_2) 7.5*[1:colormap_size-c_aux_2]])];
colormap_vec([1 2],:) = colormap_vec([2 1],:);
colormap_vec(1,:) = colormap_vec(1,:) + colormap_vec(2,:);
colormap_vec(3,:) = colormap_vec(4,:) + colormap_vec(3,:);
colormap_vec(2,:) = colormap_vec(4,:) + colormap_vec(2,:);
colormap_vec(1,:) = colormap_vec(4,:) + colormap_vec(1,:);
colormap_vec = colormap_vec'/max(colormap_vec(:));
colormap_vec = flipud(colormap_vec);
colormap_vec = colormap_vec(:,1:3);
set(h_fig_aux,'colormap',colormap_vec);
elseif evalin('base','zef.inv_colormap') == 6
c_aux_1 = floor(colormap_size/3);
c_aux_2 = floor(2*colormap_size/3);
c_aux_3 = floor(colormap_size/2);
colormap_vec = [([20*[c_aux_3:-1:1] zeros(1,colormap_size-c_aux_3)]); ([15*[1: c_aux_1] 15*[c_aux_2-c_aux_1:-1:1] zeros(1,colormap_size-c_aux_2)]) ; ([zeros(1,c_aux_1) 7*[1:c_aux_2-c_aux_1] 7*[colormap_size-c_aux_2:-1:1]]);([zeros(1,c_aux_2) 10.5*[1:colormap_size-c_aux_2]])];
colormap_vec(3,:) = colormap_vec(4,:) + colormap_vec(3,:);
colormap_vec(2,:) = colormap_vec(4,:) + colormap_vec(2,:);
colormap_vec(1,:) = colormap_vec(4,:) + colormap_vec(1,:);
colormap_vec = colormap_vec'/max(colormap_vec(:));
colormap_vec = flipud(colormap_vec);
colormap_vec = colormap_vec(:,1:3);
set(h_fig_aux,'colormap',colormap_vec);
elseif evalin('base','zef.inv_colormap') == 7
c_aux_1 = floor(colormap_size/2);    
colormap_vec = [10*[c_aux_1:-1:1] zeros(1,colormap_size-c_aux_1); 3*[1: c_aux_1] 3*[colormap_size-c_aux_1:-1:1]; zeros(1,c_aux_1) 3.8*[1:colormap_size-c_aux_1]];
colormap_vec([1 2],:) = colormap_vec([2 1],:);
colormap_vec(1,:) = colormap_vec(1,:) + colormap_vec(2,:);
colormap_vec(1,:) = colormap_vec(3,:) + colormap_vec(1,:);
colormap_vec(2,:) = colormap_vec(3,:) + colormap_vec(2,:);
colormap_vec = colormap_vec'/max(colormap_vec(:));
colormap_vec = flipud(colormap_vec);
set(h_fig_aux,'colormap',colormap_vec);
elseif evalin('base','zef.inv_colormap') == 8
c_aux_1 = floor(colormap_size/2);    
colormap_vec = [10*[c_aux_1:-1:1] zeros(1,colormap_size-c_aux_1); 3*[1: c_aux_1] 3*[colormap_size-c_aux_1:-1:1]; zeros(1,c_aux_1) 3.8*[1:colormap_size-c_aux_1]];
colormap_vec([2 3],:) = colormap_vec([3 2],:);
colormap_vec(1,:) = colormap_vec(2,:) + colormap_vec(1,:);
colormap_vec(3,:) = colormap_vec(2,:) + colormap_vec(3,:);
colormap_vec = colormap_vec+100;
colormap_vec = colormap_vec'/max(colormap_vec(:));
colormap_vec = flipud(colormap_vec);
set(h_fig_aux,'colormap',colormap_vec);
elseif evalin('base','zef.inv_colormap') == 9
c_aux_1 = floor(colormap_size/2);    
colormap_vec = [10*[c_aux_1:-1:1] zeros(1,colormap_size-c_aux_1); 2*[1: c_aux_1] 2*[colormap_size-c_aux_1:-1:1]; zeros(1,c_aux_1) 3.8*[1:colormap_size-c_aux_1]];
colormap_vec(1,:) = colormap_vec(3,:) + colormap_vec(1,:);
colormap_vec(2,:) = colormap_vec(3,:) + colormap_vec(2,:);
colormap_vec = colormap_vec+100;
colormap_vec = colormap_vec'/max(colormap_vec(:));
colormap_vec = flipud(colormap_vec);
set(h_fig_aux,'colormap',colormap_vec);
elseif evalin('base','zef.inv_colormap') == 10
c_aux_1 = floor(colormap_size/2);    
colormap_vec = [10*[c_aux_1:-1:1] zeros(1,colormap_size-c_aux_1); 8*[1: c_aux_1] 8*[colormap_size-c_aux_1:-1:1]; zeros(1,c_aux_1) 5*[1:colormap_size-c_aux_1]];
colormap_vec = colormap_vec+100;
colormap_vec = colormap_vec'/max(colormap_vec(:));
colormap_vec = flipud(colormap_vec);
set(h_fig_aux,'colormap',colormap_vec);
elseif evalin('base','zef.inv_colormap') == 11
colormap_vec = [(colormap_size/5)^3 + colormap_size^2*[1 : colormap_size] ; (colormap_size/2)^3 + ((colormap_size)/2)*[1:colormap_size].^2 ; ...
    (0.7*colormap_size)^3+(0.5*colormap_size)^2*[1:colormap_size]];
colormap_vec = colormap_vec'/max(colormap_vec(:));
set(h_fig_aux,'colormap',colormap_vec);
elseif evalin('base','zef.inv_colormap') == 12
colormap_vec = [[1:colormap_size] ; 0.5*[1:colormap_size] ; 0.5*[colormap_size:-1:1] ];
colormap_vec = colormap_vec'/max(colormap_vec(:));
set(h_fig_aux,'colormap',colormap_vec);
elseif evalin('base','zef.inv_colormap') == 13
set(h_fig_aux,'colormap',evalin('base','zef.parcellation_colormap'));
end

axes(h_axes_image); set(h_fig_aux,'visible','on');
h_surf_2 = trimesh(surface_triangles(I_3,:),nodes(:,1),nodes(:,2),nodes(:,3),reconstruction);
set(h_surf_2,'edgecolor','none','facecolor','flat','facelighting','flat','CDataMapping','scaled');
set(h_axes_image,'CLim',[min_rec max_rec]);
set(h_surf_2,'specularstrength',0.2);
set(h_surf_2,'specularexponent',0.8);
set(h_surf_2,'SpecularColorReflectance',0.8);
set(h_surf_2,'diffusestrength',1);
set(h_surf_2,'ambientstrength',1);
if evalin('base','zef.brain_transparency') < 1
f_alpha_aux = zeros(size(nodes,1),1);
I_tr = I_3;
if evalin('base','zef.inv_scale') == 1
r_alpha_aux = (reconstruction-min(reconstruction))/(max(reconstruction)-min(reconstruction));
else
r_alpha_aux = abs(reconstruction)/max(abs(reconstruction));
end
r_alpha_aux= max(0,r_alpha_aux-min(r_alpha_aux));
r_alpha_aux = r_alpha_aux/max(r_alpha_aux);
f_alpha_aux(surface_triangles(I_tr,1)) = r_alpha_aux/3;
f_alpha_aux(surface_triangles(I_tr,2)) = f_alpha_aux(surface_triangles(I_tr,2)) + r_alpha_aux/3;
f_alpha_aux(surface_triangles(I_tr,3)) = f_alpha_aux(surface_triangles(I_tr,3)) + r_alpha_aux/3; 
set(h_surf_2,'FaceVertexAlpha',max(evalin('base','zef.brain_transparency'),f_alpha_aux));
set(h_surf_2,'FaceAlpha','interp');
set(h_surf_2,'AlphaDataMapping','none'); 
end

if ismember(evalin('base','zef.visualization_type'),[2])
h_colorbar = colorbar('EastOutside','Position',[0.92 0.647 0.01 0.29]);
end

lighting phong;

if evalin('base','zef.visualization_type') == 2
%h_colorbar = colorbar('EastOutside','Position',[0.94 0.675 0.01 0.29],'fontsize',1500);
h_axes_hist = axes('position',[0.03 0.035 0.2 0.1],'visible','off');
h_bar = bar(h_axes_hist,b_hist+(max_rec-min_rec)/(2*50),a_hist,'hist');
set(h_bar,'facecolor',[0.5 0.5 0.5]);
set(h_axes_hist,'ytick',[]);
hist_ylim =  [0 log10(length_reconstruction)];
set(h_axes_hist,'ylim',hist_ylim);
set(h_axes_hist,'fontsize',1500);
set(h_axes_hist,'xlim',[min_rec max_rec]);
set(h_axes_hist,'linewidth',200);
set(h_axes_hist,'ticklength',[0 0]);
end
end
axes(h_axes_image);set(h_fig_aux,'visible','on');

for i = 1 : n_compartments

if visible_vec(i)    
I_2 = find(johtavuus(aux_ind) == i);
I_3 = find(ismember(tetra_ind,I_2));
color_str = color_cell{i};
if not(isempty(I_3))
h_surf = trimesh(surface_triangles(I_3,:),nodes(:,1),nodes(:,2),nodes(:,3),'edgecolor','none','facecolor',color_str,'facelighting','flat');
set(h_surf,'specularstrength',0.1);
set(h_surf,'diffusestrength',0.5);
set(h_surf,'ambientstrength',0.85);
if not(evalin('base','zef.visualization_type') == 2) || not(i == aux_wm_ind )
set(h_surf,'facealpha',evalin('base','zef.layer_transparency'));
end
lighting phong;
end
end
end


view(evalin('base','zef.azimuth'),evalin('base','zef.elevation'));
axis('image');
camva(evalin('base','zef.cam_va'));
if evalin('base','zef.axes_visible')
set(h_axes_image,'visible','on');
set(h_axes_image,'xGrid','on');
set(h_axes_image,'yGrid','on');
set(h_axes_image,'zGrid','on');
else
set(h_axes_image,'visible','off');
set(h_axes_image,'xGrid','off');
set(h_axes_image,'yGrid','off');
set(h_axes_image,'zGrid','off');
end

  if evalin('base','zef.visualization_type') == 2
  h_axes_text = axes('position',[0.03 0.94 0.01 0.05],'visible','off');
  h_text = text(0, 0.5, ['Time: ' num2str(evalin('base','zef.inv_time_1') + evalin('base','zef.inv_time_2')/2 + frame_step*(f_ind - 1)*evalin('base','zef.inv_time_3'),'%0.6f') ' s, Frame: ' num2str(f_ind) ' / ' num2str(length_reconstruction_cell) '.']);
  set(h_text,'visible','on','fontsize',1500);
  end
  
  drawnow;
      
if iscell(evalin('base','zef.reconstruction')) &  evalin('base','zef.visualization_type') == 2

  if is_video  
  if file_index == 1; 
  print(h_fig_aux,'-djpeg95','-r1',[file_path  file_name(1:end-4) '_' int2str(f_ind) file_name(end-3:end)]); 
  elseif file_index ==2; 
  print(h_fig_aux,'-dtiff','-r1',[file_path  file_name(1:end-4) '_' int2str(f_ind) file_name(end-3:end)]); 
  elseif file_index ==3;
  print(h_fig_aux,'-dpng','-r1',[file_path  file_name(1:end-4) '_' int2str(f_ind) file_name(end-3:end)]); 
  elseif file_index ==4;
  bmp_file_temp = [file_path  file_name(1:end-4) '_temp.bmp'];
  print(h_fig_aux,'-dbmp','-r1',bmp_file_temp);
  [movie_frame] = imread(bmp_file_temp,'bmp');
  h_frame = im2frame(movie_frame);
  writeVideo(h_aviobj,h_frame);
  delete(bmp_file_temp);
  end;
  else
  if file_index == 1; 
  print(h_fig_aux,'-djpeg95','-r1',[file_path  file_name]); 
  elseif file_index ==2; 
  print(h_fig_aux,'-dtiff','-r1',[file_path  file_name]); 
  elseif file_index==3; 
  print(h_fig_aux,'-dpng','-r1',[file_path  file_name]); 
  end;
  end;
 
else
camva(c_va);
campos(c_pos);
camtarget(c_ta);
camproj(c_p); 
camup(c_u);
drawnow;
if file_index == 1; 
print(h_fig_aux,[file_path file_name],'-djpeg95','-r1'); 
elseif file_index ==2; 
print(h_fig_aux,[file_path file_name],'-dtiff','-r1'); 
elseif file_index ==3; 
print(h_fig_aux,[file_path file_name],'-dpng','-r1'); 
end;
end

end

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
tic;
for f_ind = frame_start + frame_step : frame_step : frame_stop
pause(1/30);
f_ind_aux = f_ind_aux + 1;
if  iscell(evalin('base','zef.reconstruction')) & evalin('base','zef.visualization_type') == 2    
waitbar(f_ind_aux/number_of_frames,h_waitbar,['Frame ' int2str(f_ind_aux) ' of ' int2str(number_of_frames) '.']); 
set(h_waitbar,'handlevisibility','off');
end
delete(h_bar);
delete(h_text);
axes(h_axes_image);set(h_fig_aux,'visible','on');
delete(h_surf_2);
hold on;


%******************************************************
if iscell(evalin('base','zef.reconstruction'))
reconstruction = evalin('base',['zef.reconstruction{' int2str(f_ind) '}']);
else
reconstruction = evalin('base','zef.reconstruction');  
end
reconstruction = reconstruction(:);  
reconstruction = reshape(reconstruction,3,length(reconstruction)/3);
reconstruction = sqrt(sum(reconstruction.^2))';
reconstruction = sum(reconstruction(s_i_ind),2)/4;
if evalin('base','zef.inv_scale') == 1
reconstruction = 10*log10(max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
elseif evalin('base','zef.inv_scale') == 2
reconstruction = (max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
elseif evalin('base','zef.inv_scale') == 3
reconstruction = sqrt(max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
end
[a_hist, b_hist] = hist(reconstruction,min_rec:(max_rec-min_rec)/50:max_rec - (max_rec - min_rec)/50);
a_hist = max(0,real(log10(a_hist)));
length_reconstruction = length(reconstruction);
%******************************************************

if iscell(evalin('base','zef.reconstruction'))
reconstruction = evalin('base',['zef.reconstruction{' int2str(f_ind) '}']);
else
reconstruction = evalin('base','zef.reconstruction');  
end
reconstruction = reconstruction(:);  
reconstruction = reshape(reconstruction,3,length(reconstruction)/3);

if ismember(evalin('base','zef.reconstruction_type'),[1 7])
reconstruction = sqrt(sum(reconstruction.^2))';
elseif evalin('base','zef.reconstruction_type') == 6
reconstruction = (1/sqrt(3))*sum(reconstruction)';
end
if ismember(evalin('base','zef.reconstruction_type'), [1 6 7])
reconstruction = sum(reconstruction(s_i_ind),2)/4;
reconstruction = reconstruction(I_2_b_rec);
reconstruction = reconstruction(I_2_rec(I_1_rec));
end

if ismember(evalin('base','zef.reconstruction_type'), [2 3 4 5])
rec_x = reconstruction(1,:)';
rec_y = reconstruction(2,:)';
rec_z = reconstruction(3,:)';
rec_x = sum(rec_x(s_i_ind),2)/4;
rec_y = sum(rec_y(s_i_ind),2)/4;
rec_z = sum(rec_z(s_i_ind),2)/4;
rec_x = rec_x(I_2_b_rec);
rec_y = rec_y(I_2_b_rec);
rec_z = rec_z(I_2_b_rec);
rec_x = rec_x(I_2_rec(I_1_rec));
rec_y = rec_y(I_2_rec(I_1_rec));
rec_z = rec_z(I_2_rec(I_1_rec));
end 

if ismember(evalin('base','zef.reconstruction_type'), [2 3 4 5])
reconstruction = sqrt((rec_x.*n_vec_aux(:,1)).^2 + (rec_y.*n_vec_aux(:,2)).^2 + (rec_z.*n_vec_aux(:,3)).^2);
end

if evalin('base','zef.reconstruction_type') == 3
reconstruction = sqrt((rec_x - rec_x.*abs(n_vec_aux(:,1))).^2 + (rec_y - rec_y.*abs(n_vec_aux(:,2))).^2 + (rec_z - rec_z.*abs(n_vec_aux(:,3))).^2);
end

if evalin('base','zef.reconstruction_type') == 4
aux_rec = rec_x.*n_vec_aux(:,1) + rec_y.*n_vec_aux(:,2) + rec_z.*n_vec_aux(:,3);
I_aux_rec = find(aux_rec > 0);
reconstruction(I_aux_rec) = 0;
%reconstruction = reconstruction./max(abs(reconstruction(:)));
end

if evalin('base','zef.reconstruction_type') == 5 
aux_rec = rec_x.*n_vec_aux(:,1) + rec_y.*n_vec_aux(:,2) + rec_z.*n_vec_aux(:,3);
I_aux_rec = find(aux_rec <= 0);
reconstruction(I_aux_rec) = 0;
%reconstruction = reconstruction./max(abs(reconstruction(:)));
end

if ismember(evalin('base','zef.reconstruction_type'), [2 3 4 5 7])
reconstruction = smooth_field(surface_triangles(I_3_rec,:), reconstruction, size(nodes,1),3);
end

if not(ismember(evalin('base','zef.reconstruction_type'), [6]))
if evalin('base','zef.inv_scale') == 1
reconstruction = 10*log10(max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
elseif evalin('base','zef.inv_scale') == 2
reconstruction = (max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
elseif evalin('base','zef.inv_scale') == 3
reconstruction = sqrt(max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
end
end

if evalin('base','zef.use_parcellation')
reconstruction = reconstruction.*reconstruction_p_2;
end

h_surf_2 = trimesh(surface_triangles(I_3_rec,:),nodes(:,1),nodes(:,2),nodes(:,3),reconstruction);
set(h_surf_2,'edgecolor','none','facecolor','flat','facelighting','flat','CDataMapping','scaled');
set(gca,'CLim',[min_rec max_rec]); 
set(h_surf_2,'specularstrength',0.2);
set(h_surf_2,'specularexponent',0.8);
set(h_surf_2,'SpecularColorReflectance',0.8);
set(h_surf_2,'diffusestrength',1);
set(h_surf_2,'ambientstrength',1);
if evalin('base','zef.brain_transparency') < 1
%f_alpha_aux = zeros(size(nodes,1),1);
if evalin('base','zef.inv_scale') == 1
r_alpha_aux = (reconstruction-min(reconstruction))/(max(reconstruction)-min(reconstruction));
else
r_alpha_aux = abs(reconstruction)/max(abs(reconstruction));
end
r_alpha_aux= max(0,r_alpha_aux-min(r_alpha_aux));
r_alpha_aux = r_alpha_aux/max(r_alpha_aux);
f_alpha_aux(surface_triangles(I_tr,1)) = r_alpha_aux/3;
f_alpha_aux(surface_triangles(I_tr,2)) = f_alpha_aux(surface_triangles(I_tr,2)) + r_alpha_aux/3;
f_alpha_aux(surface_triangles(I_tr,3)) = f_alpha_aux(surface_triangles(I_tr,3)) + r_alpha_aux/3; 
set(h_surf_2,'FaceVertexAlpha',max(evalin('base','zef.brain_transparency'),f_alpha_aux));
set(h_surf_2,'FaceAlpha','interp');
set(h_surf_2,'AlphaDataMapping','none'); 
end
lighting phong;
camorbit(frame_step*evalin('base','zef.orbit_1')/15,frame_step*evalin('base','zef.orbit_2')/15);

if evalin('base','zef.visualization_type') == 2
h_bar = bar(h_axes_hist,b_hist+(max_rec-min_rec)/(2*50),a_hist,'hist');
set(h_bar,'facecolor',[0.5 0.5 0.5]);
set(h_axes_hist,'ytick',[]);
hist_ylim =  [0 log10(length_reconstruction)];
set(h_axes_hist,'ylim',hist_ylim);
set(h_axes_hist,'fontsize',1500);
set(h_axes_hist,'xlim',[min_rec max_rec]);
set(h_axes_hist,'linewidth',200);
set(h_axes_hist,'ticklength',[0 0]);

  axes(h_axes_text);set(h_fig_aux,'visible','on');
  h_text = text(0, 0.5, ['Time: ' num2str(evalin('base','zef.inv_time_1') + evalin('base','zef.inv_time_2')/2 + frame_step*(f_ind - 1)*evalin('base','zef.inv_time_3'),'%0.6f') ' s, Frame: ' num2str(f_ind) ' / ' num2str(length_reconstruction_cell) '.']);
  set(h_text,'visible','on','fontsize',1500);
  
end
 
  if file_index == 1; 
  print(h_fig_aux,'-djpeg95','-r1',[file_path  file_name(1:end-4) '_' int2str(f_ind) file_name(end-3:end)]); 
  elseif file_index ==2; 
  print(h_fig_aux,'-dtiff','-r1',[file_path  file_name(1:end-4) '_' int2str(f_ind) file_name(end-3:end)]); 
  elseif file_index ==3; 
  print(h_fig_aux,'-dpng','-r1',[file_path  file_name(1:end-4) '_' int2str(f_ind) file_name(end-3:end)]); 
  elseif file_index == 4;
  print(h_fig_aux,'-dbmp','-r1',bmp_file_temp);
  [movie_frame] = imread(bmp_file_temp,'bmp');
  h_frame = im2frame(movie_frame);
  writeVideo(h_aviobj,h_frame); 
  delete(bmp_file_temp);
  end;

end


%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



if iscell(evalin('base','zef.reconstruction')) &  evalin('base','zef.visualization_type') == 2 & is_video & file_index == 4
close(h_aviobj);
warning off; delete('avi_file'); warning on
movefile(avi_file_temp, avi_file, 'f');
end
close(h_fig_aux);
if  iscell(evalin('base','zef.reconstruction')) &  evalin('base','zef.visualization_type') == 2    
close(h_waitbar);     
end
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************

else

number_of_frames = evalin('base','zef.number_of_frames');    
file_index = evalin('base','zef.file_index');
file_name = evalin('base','zef.file');
file_path = evalin('base','zef.file_path');    
    
c_va = camva(evalin('base','zef.h_axes1'));
c_pos = campos(evalin('base','zef.h_axes1'));
c_ta = camtarget(evalin('base','zef.h_axes1'));
c_p = camproj(evalin('base','zef.h_axes1')); 
c_u = camup(evalin('base','zef.h_axes1'));


if iscell(evalin('base','zef.reconstruction'))
length_reconstruction_cell = length(evalin('base','zef.reconstruction'));
else
length_reconstruction_cell = 1;
end
frame_start = evalin('base','zef.frame_start');
frame_stop = evalin('base','zef.frame_stop');
frame_step = evalin('base','zef.frame_step');
if frame_start == 0
    frame_start = 1;
end
if frame_stop == 0
    frame_stop = length_reconstruction_cell;
end
frame_start = max(frame_start,1);
frame_start = min(length_reconstruction_cell,frame_start);
frame_stop = max(frame_stop,1);
frame_stop = min(length_reconstruction_cell,frame_stop);
if length([frame_start : frame_step : frame_stop]) > 1 
is_video = 1;
else 
is_video = 0;
end

if iscell(evalin('base','zef.reconstruction')) &  evalin('base','zef.visualization_type') == 3 & is_video & file_index == 4
avi_file_temp = [file_path file_name(1:end-4) '_temp.avi'];
avi_file = [file_path file_name];
video_quality = str2num(evalin('base','zef.video_codec'));
h_aviobj = VideoWriter(avi_file_temp); 
h_aviobj.Quality = video_quality;
open(h_aviobj);
end

if ismember(evalin('base','zef.visualization_type'),[3, 4])
s_i_ind = evalin('base','zef.source_interpolation_ind{2}');
s_i_ind_2 =  evalin('base','zef.source_interpolation_ind{1}');
end

if evalin('base','zef.use_parcellation')
selected_list = evalin('base','zef.parcellation_selected');
p_i_ind = evalin('base','zef.parcellation_interp_ind');
end

if ismember(evalin('base','zef.visualization_type'),[3])
max_abs_reconstruction = 0;
min_rec = Inf;
max_rec = -Inf;
if iscell(evalin('base','zef.reconstruction'))
length_reconstruction_cell = length(evalin('base','zef.reconstruction'));
frame_start = evalin('base','zef.frame_start');
frame_stop = evalin('base','zef.frame_stop');
frame_step = evalin('base','zef.frame_step');
if frame_start == 0
frame_start = 1;
end
if frame_stop == 0
frame_stop = length_reconstruction_cell;
end
frame_start = max(frame_start,1);
frame_start = min(length_reconstruction_cell,frame_start);
frame_stop = max(frame_stop,1);
frame_stop = min(length_reconstruction_cell,frame_stop);
number_of_frames = length([frame_start : frame_step : frame_stop]);
for f_ind = frame_start : frame_step : frame_stop
reconstruction = single(evalin('base',['zef.reconstruction{' int2str(f_ind) '}']));
reconstruction = reconstruction(:);  
reconstruction = reshape(reconstruction,3,length(reconstruction)/3);
if ismember(evalin('base','zef.reconstruction_type'), 6)
reconstruction = (1/sqrt(3))*sum(reconstruction)';
else    
reconstruction = sqrt(sum(reconstruction.^2))';
end
reconstruction = sum(reconstruction(s_i_ind_2),2)/4;
max_abs_reconstruction = max([max_abs_reconstruction ; (reconstruction(:))]);
min_rec = min([min_rec ; (reconstruction(:))]);
max_rec = max_abs_reconstruction;
end
if not(ismember(evalin('base','zef.reconstruction_type'), [6]))
if evalin('base','zef.inv_scale') == 1
min_rec = 10*log10(max(min_rec/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
max_rec = 10*log10(max(max_rec/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
elseif evalin('base','zef.inv_scale') == 2
min_rec = (max(min_rec/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
max_rec = (max(max_rec/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
elseif evalin('base','zef.inv_scale') == 3
min_rec = sqrt(max(min_rec/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
max_rec = sqrt(max(max_rec/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
end
end
else
frame_start = 1;
frame_stop = 1;
frame_step = 1;
number_of_frames = 1;    
reconstruction = evalin('base','zef.reconstruction');
reconstruction = reconstruction(:);  
reconstruction = reshape(reconstruction,3,length(reconstruction)/3);
if ismember(evalin('base','zef.reconstruction_type'), 6)
reconstruction = (1/sqrt(3))*sum(reconstruction)';
else    
reconstruction = sqrt(sum(reconstruction.^2))';
end
reconstruction = sum(reconstruction(s_i_ind_2),2)/4;
max_abs_reconstruction = max([max_abs_reconstruction ; (reconstruction(:))]);
min_rec = min([min_rec ; (reconstruction(:))]);
max_rec = max_abs_reconstruction;
if not(ismember(evalin('base','zef.reconstruction_type'), [6]))
if evalin('base','zef.inv_scale') == 1
min_rec = 10*log10(max(min_rec/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
max_rec = 10*log10(max(max_rec/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
elseif evalin('base','zef.inv_scale') == 2
min_rec = (max(min_rec/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
max_rec = (max(max_rec/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
elseif evalin('base','zef.inv_scale') == 3
min_rec = sqrt(max(min_rec/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
max_rec = sqrt(max(max_rec/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
end   
end 
end
end
   

cb_done = 0;
i = 0;
length_reuna = 0;
sigma_vec = [];
priority_vec = [];
visible_vec = [];
color_cell = cell(0);
aux_brain_ind = [];
for k = 1 : 18   
switch k
    case 1
        var_0 = 'zef.d1_on';
        var_1 = 'zef.d1_sigma';
        var_2 = 'zef.d1_priority';
        var_3 = 'zef.d1_visible';
    color_str = evalin('base','zef.d1_color');
     case 2
        var_0 = 'zef.d2_on';
        var_1 = 'zef.d2_sigma';   
        var_2 = 'zef.d2_priority';
        var_3 = 'zef.d2_visible';
        color_str = evalin('base','zef.d2_color');
     case 3
        var_0 = 'zef.d3_on';
        var_1 = 'zef.d3_sigma';   
        var_2 = 'zef.d3_priority';
        var_3 = 'zef.d3_visible';
        color_str = evalin('base','zef.d3_color');
     case 4
        var_0 = 'zef.d4_on';
        var_1 = 'zef.d4_sigma';   
        var_2 = 'zef.d4_priority';
        var_3 = 'zef.d4_visible';
        color_str = evalin('base','zef.d4_color');
     case 5
        var_0 = 'zef.d5_on';
        var_1 = 'zef.d5_sigma';
        var_2 = 'zef.d5_priority';
        var_3 = 'zef.d5_visible';
        color_str = evalin('base','zef.d5_color');
     case 6
        var_0 = 'zef.d6_on';
        var_1 = 'zef.d6_sigma';   
        var_2 = 'zef.d6_priority';
        var_3 = 'zef.d6_visible';
        color_str = evalin('base','zef.d6_color');
     case 7
        var_0 = 'zef.d7_on';
        var_1 = 'zef.d7_sigma';   
        var_2 = 'zef.d7_priority';
        var_3 = 'zef.d7_visible';
        color_str = evalin('base','zef.d7_color');
     case 8
        var_0 = 'zef.d8_on';
        var_1 = 'zef.d8_sigma';   
        var_2 = 'zef.d8_priority';
        var_3 = 'zef.d8_visible';
        color_str = evalin('base','zef.d8_color');
     case 9
        var_0 = 'zef.d9_on';
        var_1 = 'zef.d9_sigma';   
        var_2 = 'zef.d9_priority';
        var_3 = 'zef.d9_visible';
        color_str = evalin('base','zef.d9_color');
             case 10
        var_0 = 'zef.d10_on';
        var_1 = 'zef.d10_sigma';   
        var_2 = 'zef.d10_priority';
        var_3 = 'zef.d10_visible';
        color_str = evalin('base','zef.d10_color');
     case 11
        var_0 = 'zef.d11_on';
        var_1 = 'zef.d11_sigma';   
        var_2 = 'zef.d11_priority';
        var_3 = 'zef.d11_visible';
        color_str = evalin('base','zef.d11_color');
     case 12
        var_0 = 'zef.d12_on';
        var_1 = 'zef.d12_sigma';   
        var_2 = 'zef.d12_priority';
        var_3 = 'zef.d12_visible';
        color_str = evalin('base','zef.d12_color');
     case 13
        var_0 = 'zef.d13_on';
        var_1 = 'zef.d13_sigma';   
        var_2 = 'zef.d13_priority';
        var_3 = 'zef.d13_visible';
        color_str = evalin('base','zef.d13_color');
    case 14
        var_0 = 'zef.w_on';
        var_1 = 'zef.w_sigma';    
        var_2 = 'zef.w_priority';
        var_3 = 'zef.w_visible';
        color_str = evalin('base','zef.w_color');
    case 15
        var_0 = 'zef.g_on';
        var_1 = 'zef.g_sigma';
        var_2 = 'zef.g_priority';
        var_3 = 'zef.g_visible';
        color_str = evalin('base','zef.g_color');
    case 16
        var_0 = 'zef.c_on';
        var_1 = 'zef.c_sigma';
        var_2 = 'zef.c_priority';
        var_3 = 'zef.c_visible';
        color_str = evalin('base','zef.c_color');
     case 17
        var_0 = 'zef.sk_on';
        var_1 = 'zef.sk_sigma';
        var_2 = 'zef.sk_priority';
        var_3 = 'zef.sk_visible';
        color_str = evalin('base','zef.sk_color');
     case 18
        var_0 = 'zef.sc_on';
        var_1 = 'zef.sc_sigma';
        var_2 = 'zef.sc_priority';
        var_3 = 'zef.sc_visible';
        color_str = evalin('base','zef.sc_color');
     end
on_val = evalin('base',var_0);      
sigma_val = evalin('base',var_1);  
priority_val = evalin('base',var_2);
visible_val = evalin('base',var_3);
if on_val
i = i + 1;
sigma_vec(i,1) = sigma_val;
priority_vec(i,1) = priority_val;
color_cell{i} = color_str;
visible_vec(i,1) = i*visible_val;

if k == 1 && evalin('base','zef.d1_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 2 && evalin('base','zef.d2_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 3 && evalin('base','zef.d3_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 4 && evalin('base','zef.d4_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 5 && evalin('base','zef.d5_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 6 && evalin('base','zef.d6_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 7 && evalin('base','zef.d7_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 8 && evalin('base','zef.d8_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 9 && evalin('base','zef.d9_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 10 && evalin('base','zef.d10_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 11 && evalin('base','zef.d11_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 12 && evalin('base','zef.d12_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 13 && evalin('base','zef.d13_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 14 && evalin('base','zef.wm_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 15 && evalin('base','zef.g_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 16 && evalin('base','zef.c_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 17 && evalin('base','zef.sk_sources');
    aux_brain_ind = [aux_brain_ind i];
end
if k == 18 && evalin('base','zef.sc_sources');
    aux_brain_ind = [aux_brain_ind i];
end
end
end


snapshot_resolution = round([evalin('base','zef.snapshot_vertical_resolution') evalin('base','zef.snapshot_horizontal_resolution')]);

h_fig_aux = figure;set(h_fig_aux,'visible','on'); 
clf;
set(h_fig_aux,'renderer','opengl');
set(h_fig_aux,'paperunits','inches');
set(h_fig_aux,'papersize',snapshot_resolution);
set(h_fig_aux,'paperposition',[0 0 fliplr(snapshot_resolution)]);
light('Position',[0 0 1],'Style','infinite');
light('Position',[0 0 -1],'Style','infinite');
h_axes_image = get(h_fig_aux,'currentaxes');
hold on;
sensors = evalin('base','zef.sensors');
reuna_p = evalin('base','zef.reuna_p');
reuna_t = evalin('base','zef.reuna_t');  
if evalin('base','zef.cp_on') || evalin('base','zef.cp2_on') || evalin('base','zef.cp3_on')
    for i = 1 : length(reuna_t)
        reuna_t{i} = uint32(reuna_t{i});
triangle_c{i} = (1/3)*(reuna_p{i}(reuna_t{i}(:,1),:) + reuna_p{i}(reuna_t{i}(:,2),:) + reuna_p{i}(reuna_t{i}(:,3),:));
    end
end

aux_ind_1 = [];
aux_ind_2 = cell(1,length(reuna_t));

if evalin('base','zef.cp_on')
cp_a = evalin('base','zef.cp_a');
cp_b = evalin('base','zef.cp_b');
cp_c = evalin('base','zef.cp_c');
cp_d = evalin('base','zef.cp_d');
% if cp_a ~= 0 | cp_b ~=0    
% light('Position',[-cp_a -cp_b -cp_b],'Style','infinite');
% end
if not(isempty(aux_ind_1))
aux_ind_1 = intersect(aux_ind_1,find(sum(sensors(:,1:3).*repmat([cp_a cp_b cp_c],size(sensors,1),1),2) >= cp_d));
else
aux_ind_1 = find(sum(sensors(:,1:3).*repmat([cp_a cp_b cp_c],size(sensors,1),1),2) >= cp_d);
end
for i = 1 : length(reuna_t)
    if not(isempty(aux_ind_2{i}))
aux_ind_2{i} = intersect(aux_ind_2{i},find(sum(triangle_c{i}.*repmat([cp_a cp_b cp_c],size(triangle_c{i},1),1),2) >= cp_d));
else
aux_ind_2{i} = find(sum(triangle_c{i}.*repmat([cp_a cp_b cp_c],size(triangle_c{i},1),1),2) >= cp_d);
    end
end    
end


if evalin('base','zef.cp2_on')
cp2_a = evalin('base','zef.cp2_a');
cp2_b = evalin('base','zef.cp2_b');
cp2_c = evalin('base','zef.cp2_c');
cp2_d = evalin('base','zef.cp2_d');
% if cp2_a ~= 0 | cp2_b ~=0    
% light('Position',[-cp2_a -cp2_b -cp2_b],'Style','infinite');
% end
if not(isempty(aux_ind_1))
aux_ind_1 = intersect(aux_ind_1,find(sum(sensors(:,1:3).*repmat([cp2_a cp2_b cp2_c],size(sensors,1),1),2) >= cp2_d));
else
aux_ind_1 = find(sum(sensors(:,1:3).*repmat([cp2_a cp2_b cp2_c],size(sensors,1),1),2) >= cp2_d);
end
for i = 1 : length(reuna_t)
    if not(isempty(aux_ind_2{i}))
aux_ind_2{i} = intersect(aux_ind_2{i},find(sum(triangle_c{i}.*repmat([cp2_a cp2_b cp2_c],size(triangle_c{i},1),1),2) >= cp2_d));
else
aux_ind_2{i} = find(sum(triangle_c{i}.*repmat([cp2_a cp2_b cp2_c],size(triangle_c{i},1),1),2) >= cp2_d);
    end
end    
end

if evalin('base','zef.cp3_on')
cp3_a = evalin('base','zef.cp3_a');
cp3_b = evalin('base','zef.cp3_b');
cp3_c = evalin('base','zef.cp3_c');
cp3_d = evalin('base','zef.cp3_d');
% if cp3_a ~= 0 | cp3_b ~=0    
% light('Position',[-cp3_a -cp3_b -cp3_b],'Style','infinite');
% end
if not(isempty(aux_ind_1))
aux_ind_1 = intersect(aux_ind_1,find(sum(sensors(:,1:3).*repmat([cp3_a cp3_b cp3_c],size(sensors,1),1),2) >= cp3_d));
else
aux_ind_1 = find(sum(sensors(:,1:3).*repmat([cp3_a cp3_b cp3_c],size(sensors,1),1),2) >= cp3_d);
end
for i = 1 : length(reuna_t)
if not(isempty(aux_ind_2{i}))
aux_ind_2{i} = intersect(aux_ind_2{i},find(sum(triangle_c{i}.*repmat([cp3_a cp3_b cp3_c],size(triangle_c{i},1),1),2) >= cp3_d));
else
aux_ind_2{i} = find(sum(triangle_c{i}.*repmat([cp3_a cp3_b cp3_c],size(triangle_c{i},1),1),2) >= cp3_d);
end
end    
end


if evalin('base','zef.cp_on') || evalin('base','zef.cp2_on') || evalin('base','zef.cp3_on')
if evalin('base','zef.cp_mode') == 1
sensors = sensors(aux_ind_1,:);
elseif evalin('base','zef.cp_mode') == 2
aux_ind_1 = setdiff([1:size(sensors,1)]',aux_ind_1);
sensors = sensors(aux_ind_1,:);   
end
for i = 1 : length(reuna_t)
if evalin('base','zef.cp_mode') == 1
reuna_t{i} = reuna_t{i}(aux_ind_2{i},:);
if ismember(evalin('base','zef.visualization_type'), [3,4])
if ismember(i, aux_brain_ind)
ab_ind = find(aux_brain_ind==i);
s_i_ind{ab_ind} = s_i_ind{ab_ind}(aux_ind_2{i},:);
if evalin('base','zef.use_parcellation')
for p_ind = selected_list
[aux_is_1, aux_is_2, aux_is_3] = intersect(p_i_ind{p_ind}{2}{ab_ind},aux_ind_2{i});
p_i_ind{p_ind}{2}{ab_ind}= aux_is_3;
end
end
end;
end
elseif evalin('base','zef.cp_mode') == 2
aux_ind_2{i} = setdiff([1:size(reuna_t{i},1)]',aux_ind_2{i});
reuna_t{i} = reuna_t{i}(aux_ind_2{i},:);   
if ismember(evalin('base','zef.visualization_type'), [3,4])
if ismember(i, aux_brain_ind)
ab_ind = find(aux_brain_ind==i);
s_i_ind{ab_ind} = s_i_ind{ab_ind}(aux_ind_2{i},:);
if evalin('base','zef.use_parcellation')
for p_ind = selected_list
[aux_is_1, aux_is_2, aux_is_3] = intersect(p_i_ind{p_ind}{2}{ab_ind},aux_ind_2{i});
p_i_ind{p_ind}{2}{ab_ind}= aux_is_3;
end
end
end;
end
elseif evalin('base','zef.cp_mode') == 3
if ismember(i, aux_brain_ind)
aux_ind_2{i} = reuna_t{i};
else
reuna_t{i} = reuna_t{i}(aux_ind_2{i},:);
end
elseif evalin('base','zef.cp_mode') == 4
if ismember(i, aux_brain_ind)
aux_ind_2{i} = reuna_t{i};
else  
aux_ind_2{i} = setdiff([1:size(reuna_t{i},1)]',aux_ind_2{i});
reuna_t{i} = reuna_t{i}(aux_ind_2{i},:); 
end
end
end
end

aux_ind_1 = [];
aux_ind_2 = cell(1,length(reuna_t));
triangle_c = cell(1,length(reuna_t));


if ismember(evalin('base','zef.imaging_method'), [1 4 5])  & size(sensors,2) == 6 
    electrode_model = 2;
elseif ismember(evalin('base','zef.imaging_method'), [1 4 5])
    electrode_model = 1;
else
    electrode_model = 0;
end

if evalin('base','zef.attach_electrodes') & electrode_model == 1
for i = 1 : size(sensors,1)
[min_val, min_ind] = min(sqrt(sum((reuna_p{end} - repmat(sensors(i,1:3),size(reuna_p{end},1),1)).^2,2)));
sensors(i,1:3) = reuna_p{end}(min_ind,:);
end
end

if evalin('base','zef.attach_electrodes') & electrode_model == 2
   center_points_aux = (1/3)*(reuna_p{end}(reuna_t{end}(:,1),:) + ...
                       reuna_p{end}(reuna_t{end}(:,2),:) + ...
                       reuna_p{end}(reuna_t{end}(:,3),:));
for i = 1 : size(sensors,1)
[min_val, min_ind] = min(sqrt(sum((reuna_p{end} - repmat(sensors(i,1:3),size(reuna_p{end},1),1)).^2,2)));
sensors(i,1:3) = reuna_p{end}(min_ind,:);
end
 sensors_aux = []; 
 for i = 1 : size(sensors,1)
 [dist_val] = (sqrt(sum((center_points_aux - repmat(sensors(i,1:3),size(center_points_aux,1),1)).^2,2)));
 dist_ind = find(dist_val <= sensors(i,4) & dist_val >= sensors(i,5)); 
 sensors_aux = [sensors_aux ; i*ones(length(dist_ind),1) reuna_t{end}(dist_ind,:)];
 end
 sensors = sensors_aux;
end

[X_s, Y_s, Z_s] = sphere(20);  
aux_scale_ind = evalin('base','zef.location_unit');
switch aux_scale_ind 
    case 1
        aux_scale_val = 1;
    case 2 
        aux_scale_val = 0.01;
    case 3 
        aux_scale_val = 0.001;
end



if evalin('base','zef.s_visible')
if electrode_model == 1 | not(ismember(evalin('base','zef.imaging_method'),[1,4,5]))
sphere_scale = 3.2*aux_scale_val;    
X_s = sphere_scale*X_s;
Y_s = sphere_scale*Y_s;
Z_s = sphere_scale*Z_s;
for i = 1 : size(sensors,1)
h = surf(sensors(i,1) + X_s, sensors(i,2) + Y_s, sensors(i,3) + Z_s);
set(h,'facecolor',evalin('base','zef.s_color'));
set(h,'edgecolor','none'); 
set(h,'specularstrength',0.3);
set(h,'diffusestrength',0.7);
set(h,'ambientstrength',0.7);
set(h,'facealpha',evalin('base','zef.layer_transparency'));
end
elseif electrode_model == 2  
h = trisurf(sensors(:,2:4),reuna_p{end}(:,1),reuna_p{end}(:,2),reuna_p{end}(:,3));
set(h,'facecolor',evalin('base','zef.s_color'));
set(h,'edgecolor',evalin('base','zef.s_color')); 
set(h,'specularstrength',0.3);
set(h,'diffusestrength',0.7);
set(h,'ambientstrength',0.7);
set(h,'facealpha',evalin('base','zef.layer_transparency'));
set(h,'edgealpha',evalin('base','zef.layer_transparency'));
end

if ismember(evalin('base','zef.imaging_method'),[2 3])
sensors(:,4:6) = sensors(:,4:6)./repmat(sqrt(sum(sensors(:,4:6).^2,2)),1,3);
h=coneplot(sensors(:,1) + 4.5*sensors(:,4),sensors(:,2) + 4.5*sensors(:,5),sensors(:,3) + 4.5*sensors(:,6),8*sensors(:,4),8*sensors(:,5),8*sensors(:,6),0,'nointerp');
set(h,'facecolor',evalin('base','zef.s_color'));
set(h,'edgecolor','none'); 
set(h,'specularstrength',0.3);
set(h,'diffusestrength',0.7);
set(h,'ambientstrength',0.7);
set(h,'facealpha',evalin('base','zef.layer_transparency'));
if size(sensors,2) == 9
sensors(:,7:9) = sensors(:,7:9)./repmat(sqrt(sum(sensors(:,7:9).^2,2)),1,3);
h=coneplot(sensors(:,1) + 4.5*sensors(:,7),sensors(:,2) + 4.5*sensors(:,8),sensors(:,3) + 4.5*sensors(:,9),8*sensors(:,7),8*sensors(:,8),8*sensors(:,9),0,'nointerp');
set(h,'facecolor', 0.9*[1 1 1]);
set(h,'edgecolor','none'); 
set(h,'specularstrength',0.3);
set(h,'diffusestrength',0.7);
set(h,'ambientstrength',0.7);
set(h,'facealpha',evalin('base','zef.layer_transparency'));
end
end
end


if ismember(evalin('base','zef.visualization_type'),[3,4]) 

    if ismember(evalin('base','zef.visualization_type'),[3])
f_ind = frame_start;
end

i = 0;

aux_brain_visible_ind = [];

for k = 1 : 18
switch k
    case 1
        on_val = evalin('base','zef.d1_on');
        visible_val = evalin('base','zef.d1_visible');
        color_str =  evalin('base','zef.d1_color');
    case 2
        on_val = evalin('base','zef.d2_on');
        visible_val = evalin('base','zef.d2_visible');
        color_str = evalin('base','zef.d2_color');
    case 3
        on_val = evalin('base','zef.d3_on');
        visible_val = evalin('base','zef.d3_visible');
        color_str = evalin('base','zef.d3_color');
    case 4
        on_val = evalin('base','zef.d4_on');
        visible_val = evalin('base','zef.d4_visible');
        color_str = evalin('base','zef.d4_color');
    case 5
        on_val = evalin('base','zef.d5_on');
        visible_val = evalin('base','zef.d5_visible');
        color_str =  evalin('base','zef.d5_color');
    case 6
        on_val = evalin('base','zef.d6_on');
        visible_val = evalin('base','zef.d6_visible');
        color_str = evalin('base','zef.d6_color');
    case 7
        on_val = evalin('base','zef.d7_on');
        visible_val = evalin('base','zef.d7_visible');
        color_str = evalin('base','zef.d7_color');
    case 8
        on_val = evalin('base','zef.d8_on');
        visible_val = evalin('base','zef.d8_visible');
        color_str = evalin('base','zef.d8_color');
    case 9
        on_val = evalin('base','zef.d9_on');
        visible_val = evalin('base','zef.d9_visible');
        color_str =  evalin('base','zef.d9_color');
    case 10
        on_val = evalin('base','zef.d10_on');
        visible_val = evalin('base','zef.d10_visible');
        color_str = evalin('base','zef.d10_color');
    case 11
        on_val = evalin('base','zef.d11_on');
        visible_val = evalin('base','zef.d11_visible');
        color_str = evalin('base','zef.d11_color');
    case 12
        on_val = evalin('base','zef.d12_on');
        visible_val = evalin('base','zef.d12_visible');
        color_str = evalin('base','zef.d12_color');
    case 13
        on_val = evalin('base','zef.d13_on');
        visible_val = evalin('base','zef.d13_visible');
        color_str = evalin('base','zef.d13_color');
    case 14
        on_val = evalin('base','zef.w_on');
        visible_val = evalin('base','zef.w_visible');
        color_str = evalin('base','zef.w_color');
    case 15
        on_val = evalin('base','zef.g_on');
        visible_val = evalin('base','zef.g_visible');
        color_str = evalin('base','zef.g_color');
    case 16
        on_val = evalin('base','zef.c_on');
        visible_val = evalin('base','zef.c_visible');
        color_str = evalin('base','zef.c_color');
    case 17
        on_val = evalin('base','zef.sk_on');
        visible_val = evalin('base','zef.sk_visible');
        color_str = evalin('base','zef.sk_color');
    case 18
        on_val = evalin('base','zef.sc_on');
        visible_val = evalin('base','zef.sc_visible');
        color_str = evalin('base','zef.sc_color');
    end
if on_val  
i = i + 1;    
if visible_val
if ismember(i, aux_brain_ind)
ab_ind = find(aux_brain_ind==i);
aux_brain_visible_ind = [aux_brain_visible_ind i];

if ismember(evalin('base','zef.visualization_type'),[3])
if  i == aux_brain_visible_ind
if  iscell(evalin('base','zef.reconstruction')) 
h_waitbar = waitbar(1/number_of_frames,['Frame ' int2str(1) ' of ' int2str(number_of_frames) '.']);    
set(h_waitbar,'handlevisibility','off');
end    
end  
end
    
colormap_size = 4096;
if evalin('base','zef.inv_colormap') == 1
c_aux_1 = floor(colormap_size/3);
c_aux_2 = floor(2*colormap_size/3);
colormap_vec_aux = [([20*[c_aux_1:-1:1] zeros(1,colormap_size-c_aux_1)]); ([15*[1: c_aux_1] 15*[c_aux_2-c_aux_1:-1:1] zeros(1,colormap_size-c_aux_2)]) ; ([zeros(1,c_aux_1) 6*[1:c_aux_2-c_aux_1] 6*[colormap_size-c_aux_2:-1:1]]);([zeros(1,c_aux_2) 7.5*[1:colormap_size-c_aux_2]])];
colormap_vec = zeros(3,size(colormap_vec_aux,2));
colormap_vec = colormap_vec + 0.52*[50*colormap_vec_aux(1,:) ; 50*colormap_vec_aux(1,:) ; 50*colormap_vec_aux(1,:)];
colormap_vec = colormap_vec + 0.5*[85*colormap_vec_aux(3,:) ; 197*colormap_vec_aux(3,:) ; 217*colormap_vec_aux(3,:)];
colormap_vec = colormap_vec + 0.1*[2*colormap_vec_aux(2,:) ; 118*colormap_vec_aux(2,:) ; 132*colormap_vec_aux(2,:)];
colormap_vec = colormap_vec + [203*colormap_vec_aux(4,:) ; 203*colormap_vec_aux(4,:) ; 100*colormap_vec_aux(4,:)];
clear colormap_vec_aux;
colormap_vec = colormap_vec'/max(colormap_vec(:));
colormap_vec = colormap_vec(:,1:3);
set(h_fig_aux,'colormap',colormap_vec);
elseif evalin('base','zef.inv_colormap') == 2
c_aux_1 = floor(colormap_size/3);
c_aux_2 = floor(2*colormap_size/3);
colormap_vec = zeros(3,colormap_size);
colormap_vec(1,:) =10*([colormap_size:-1:1]/colormap_size);
colormap_vec(2,:) = [10*( (3/2)*[c_aux_2:-1:1]/colormap_size) zeros(1,colormap_size-c_aux_2)];
colormap_vec(3,:) = [10*((3)*[c_aux_1:-1:1]/colormap_size) zeros(1,colormap_size-c_aux_1)];  
colormap_vec = colormap_vec'/max(colormap_vec(:));
colormap_vec = flipud(colormap_vec);
colormap_vec = colormap_vec(:,1:3);
colormap_vec = colormap_vec + repmat(0.2*([colormap_size:-1:1]'/colormap_size),1,3);
colormap_vec = colormap_vec/max(colormap_vec(:));
set(h_fig_aux,'colormap',colormap_vec);
elseif evalin('base','zef.inv_colormap') == 3
c_aux_1 = floor(colormap_size/3);
c_aux_2 = floor(2*colormap_size/3);
colormap_vec = zeros(3,colormap_size);
colormap_vec(2,:) = 10*([colormap_size:-1:1]/colormap_size);
colormap_vec(1,:) = [10*(3/2)*[c_aux_2:-1:1]/colormap_size zeros(1,colormap_size-c_aux_2)];
colormap_vec(3,:) = [10*(3*[c_aux_1:-1:1]/colormap_size) zeros(1,colormap_size-c_aux_1)];
%colormap_vec([1 3 2],:) = 0.75*colormap_vec([1 3 2],:) + 0.25*colormap_vec([1 2 3],:);
colormap_vec = colormap_vec'/max(colormap_vec(:));
colormap_vec = flipud(colormap_vec);
colormap_vec = colormap_vec(:,1:3);
colormap_vec = colormap_vec + repmat(0.2*([colormap_size:-1:1]'/colormap_size),1,3);
colormap_vec = colormap_vec/max(colormap_vec(:));
set(h_fig_aux,'colormap',colormap_vec);
elseif evalin('base','zef.inv_colormap') == 4
c_aux_1 = floor(colormap_size/3);
c_aux_2 = floor(2*colormap_size/3);
colormap_vec = zeros(3,colormap_size);
colormap_vec(3,:) = 10*([colormap_size:-1:1]/colormap_size);
colormap_vec(2,:) = [10*(3/2)*[c_aux_2:-1:1]/colormap_size zeros(1,colormap_size-c_aux_2)];
colormap_vec(1,:) = [10*(3*[c_aux_1:-1:1]/colormap_size) zeros(1,colormap_size-c_aux_1)];
%colormap_vec([1 3 2],:) = 0.75*colormap_vec([1 3 2],:) + 0.25*colormap_vec([1 2 3],:);
colormap_vec = colormap_vec'/max(colormap_vec(:));
colormap_vec = flipud(colormap_vec);
colormap_vec = colormap_vec(:,1:3);
colormap_vec = colormap_vec + repmat(0.2*([colormap_size:-1:1]'/colormap_size),1,3);
colormap_vec = colormap_vec/max(colormap_vec(:));
set(h_fig_aux,'colormap',colormap_vec);
elseif evalin('base','zef.inv_colormap') == 5
c_aux_1 = floor(colormap_size/3);
c_aux_2 = floor(2*colormap_size/3);
colormap_vec = [([20*[c_aux_1:-1:1] zeros(1,colormap_size-c_aux_1)]); ([15*[1: c_aux_1] 15*[c_aux_2-c_aux_1:-1:1] zeros(1,colormap_size-c_aux_2)]) ; ([zeros(1,c_aux_1) 6*[1:c_aux_2-c_aux_1] 6*[colormap_size-c_aux_2:-1:1]]);([zeros(1,c_aux_2) 7.5*[1:colormap_size-c_aux_2]])];
colormap_vec([1 2],:) = colormap_vec([2 1],:);
colormap_vec(1,:) = colormap_vec(1,:) + colormap_vec(2,:);
colormap_vec(3,:) = colormap_vec(4,:) + colormap_vec(3,:);
colormap_vec(2,:) = colormap_vec(4,:) + colormap_vec(2,:);
colormap_vec(1,:) = colormap_vec(4,:) + colormap_vec(1,:);
colormap_vec = colormap_vec'/max(colormap_vec(:));
colormap_vec = flipud(colormap_vec);
colormap_vec = colormap_vec(:,1:3);
set(h_fig_aux,'colormap',colormap_vec);
elseif evalin('base','zef.inv_colormap') == 6
c_aux_1 = floor(colormap_size/3);
c_aux_2 = floor(2*colormap_size/3);
c_aux_3 = floor(colormap_size/2);
colormap_vec = [([20*[c_aux_3:-1:1] zeros(1,colormap_size-c_aux_3)]); ([15*[1: c_aux_1] 15*[c_aux_2-c_aux_1:-1:1] zeros(1,colormap_size-c_aux_2)]) ; ([zeros(1,c_aux_1) 7*[1:c_aux_2-c_aux_1] 7*[colormap_size-c_aux_2:-1:1]]);([zeros(1,c_aux_2) 10.5*[1:colormap_size-c_aux_2]])];
colormap_vec(3,:) = colormap_vec(4,:) + colormap_vec(3,:);
colormap_vec(2,:) = colormap_vec(4,:) + colormap_vec(2,:);
colormap_vec(1,:) = colormap_vec(4,:) + colormap_vec(1,:);
colormap_vec = colormap_vec'/max(colormap_vec(:));
colormap_vec = flipud(colormap_vec);
colormap_vec = colormap_vec(:,1:3);
set(h_fig_aux,'colormap',colormap_vec);
elseif evalin('base','zef.inv_colormap') == 7
c_aux_1 = floor(colormap_size/2);    
colormap_vec = [10*[c_aux_1:-1:1] zeros(1,colormap_size-c_aux_1); 3*[1: c_aux_1] 3*[colormap_size-c_aux_1:-1:1]; zeros(1,c_aux_1) 3.8*[1:colormap_size-c_aux_1]];
colormap_vec([1 2],:) = colormap_vec([2 1],:);
colormap_vec(1,:) = colormap_vec(1,:) + colormap_vec(2,:);
colormap_vec(1,:) = colormap_vec(3,:) + colormap_vec(1,:);
colormap_vec(2,:) = colormap_vec(3,:) + colormap_vec(2,:);
colormap_vec = colormap_vec'/max(colormap_vec(:));
colormap_vec = flipud(colormap_vec);
set(h_fig_aux,'colormap',colormap_vec);
elseif evalin('base','zef.inv_colormap') == 8
c_aux_1 = floor(colormap_size/2);    
colormap_vec = [10*[c_aux_1:-1:1] zeros(1,colormap_size-c_aux_1); 3*[1: c_aux_1] 3*[colormap_size-c_aux_1:-1:1]; zeros(1,c_aux_1) 3.8*[1:colormap_size-c_aux_1]];
colormap_vec([2 3],:) = colormap_vec([3 2],:);
colormap_vec(1,:) = colormap_vec(2,:) + colormap_vec(1,:);
colormap_vec(3,:) = colormap_vec(2,:) + colormap_vec(3,:);
colormap_vec = colormap_vec+100;
colormap_vec = colormap_vec'/max(colormap_vec(:));
colormap_vec = flipud(colormap_vec);
set(h_fig_aux,'colormap',colormap_vec);
elseif evalin('base','zef.inv_colormap') == 9
c_aux_1 = floor(colormap_size/2);    
colormap_vec = [10*[c_aux_1:-1:1] zeros(1,colormap_size-c_aux_1); 2*[1: c_aux_1] 2*[colormap_size-c_aux_1:-1:1]; zeros(1,c_aux_1) 3.8*[1:colormap_size-c_aux_1]];
colormap_vec(1,:) = colormap_vec(3,:) + colormap_vec(1,:);
colormap_vec(2,:) = colormap_vec(3,:) + colormap_vec(2,:);
colormap_vec = colormap_vec+100;
colormap_vec = colormap_vec'/max(colormap_vec(:));
colormap_vec = flipud(colormap_vec);
set(h_fig_aux,'colormap',colormap_vec);
elseif evalin('base','zef.inv_colormap') == 10
c_aux_1 = floor(colormap_size/2);    
colormap_vec = [10*[c_aux_1:-1:1] zeros(1,colormap_size-c_aux_1); 8*[1: c_aux_1] 8*[colormap_size-c_aux_1:-1:1]; zeros(1,c_aux_1) 5*[1:colormap_size-c_aux_1]];
colormap_vec = colormap_vec+100;
colormap_vec = colormap_vec'/max(colormap_vec(:));
colormap_vec = flipud(colormap_vec);
set(h_fig_aux,'colormap',colormap_vec);
elseif evalin('base','zef.inv_colormap') == 11
colormap_vec = [(colormap_size/5)^3 + colormap_size^2*[1 : colormap_size] ; (colormap_size/2)^3 + ((colormap_size)/2)*[1:colormap_size].^2 ; ...
    (0.7*colormap_size)^3+(0.5*colormap_size)^2*[1:colormap_size]];
colormap_vec = colormap_vec'/max(colormap_vec(:));
set(h_fig_aux,'colormap',colormap_vec);
elseif evalin('base','zef.inv_colormap') == 12
colormap_vec = [[1:colormap_size] ; 0.5*[1:colormap_size] ; 0.5*[colormap_size:-1:1] ];
colormap_vec = colormap_vec'/max(colormap_vec(:));
set(h_fig_aux,'colormap',colormap_vec);
elseif evalin('base','zef.inv_colormap') == 13
set(h_fig_aux,'colormap',evalin('base','zef.parcellation_colormap'));
end

if ismember(evalin('base','zef.visualization_type'),[4])
    

if evalin('base','zef.use_parcellation')
reconstruction = ones(size(reuna_t{i},1),1);
for p_ind = selected_list
reconstruction(p_i_ind{p_ind}{2}{ab_ind}) = p_ind+1;
end
end
min_rec = 1; 
max_rec = size(evalin('base','zef.parcellation_colormap'),1);


else

%******************************************************
if iscell(evalin('base','zef.reconstruction')) 
reconstruction = evalin('base',['zef.reconstruction{' int2str(frame_start) '}']);
else
reconstruction = evalin('base','zef.reconstruction');
end
reconstruction = reconstruction(:);  
reconstruction = reshape(reconstruction,3,length(reconstruction)/3);
reconstruction = sqrt(sum(reconstruction.^2))';
reconstruction = sum(reconstruction(s_i_ind_2),2)/4;
if evalin('base','zef.inv_scale') == 1
reconstruction = 10*log10(max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
elseif evalin('base','zef.inv_scale') == 2
reconstruction = (max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
elseif evalin('base','zef.inv_scale') == 3
reconstruction = sqrt(max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
end
[a_hist, b_hist] = hist(reconstruction,min_rec:(max_rec-min_rec)/50:max_rec - (max_rec - min_rec)/50);
a_hist = max(0,real(log10(a_hist)));
length_reconstruction = length(reconstruction);
%******************************************************
    

if iscell(evalin('base','zef.reconstruction')) 
reconstruction = single(evalin('base',['zef.reconstruction{' int2str(frame_start) '}']));
else
reconstruction = evalin('base','zef.reconstruction');
end
reconstruction = reconstruction(:);
reconstruction = reshape(reconstruction,3,length(reconstruction)/3);

if ismember(evalin('base','zef.reconstruction_type'),[1 7])
reconstruction = sqrt(sum(reconstruction.^2))';
elseif evalin('base','zef.reconstruction_type') == 6
reconstruction = (1/sqrt(3))*sum(reconstruction)';
end
if ismember(evalin('base','zef.reconstruction_type'), [1 6 7])
reconstruction = sum(reconstruction(s_i_ind{ab_ind}),2)/3;
end

if ismember(evalin('base','zef.reconstruction_type'), [2 3 4 5])
rec_x = reconstruction(1,:)';
rec_y = reconstruction(2,:)';
rec_z = reconstruction(3,:)';
rec_x = sum(rec_x(s_i_ind{ab_ind}),2)/3;
rec_y = sum(rec_y(s_i_ind{ab_ind}),2)/3;
rec_z = sum(rec_z(s_i_ind{ab_ind}),2)/3;
n_vec_aux = cross(reuna_p{i}(reuna_t{i}(:,2),:)' - reuna_p{i}(reuna_t{i}(:,1),:)',...
reuna_p{i}(reuna_t{i}(:,3),:)' - reuna_p{i}(reuna_t{i}(:,1),:)')';
n_vec_aux = n_vec_aux./repmat(sqrt(sum(n_vec_aux.^2,2)),1,3);
end

if ismember(evalin('base','zef.reconstruction_type'), [2 3 4 5])
reconstruction = sqrt((rec_x.*n_vec_aux(:,1)).^2 + (rec_y.*n_vec_aux(:,2)).^2 + (rec_z.*n_vec_aux(:,3)).^2);
end

if evalin('base','zef.reconstruction_type') == 3
reconstruction = sqrt((rec_x - rec_x.*abs(n_vec_aux(:,1))).^2 + (rec_y - rec_y.*abs(n_vec_aux(:,2))).^2 + (rec_z - rec_z.*abs(n_vec_aux(:,3))).^2);
end

if evalin('base','zef.reconstruction_type') == 4 
aux_rec = rec_x.*n_vec_aux(:,1) + rec_y.*n_vec_aux(:,2) + rec_z.*n_vec_aux(:,3);
I_aux_rec = find(aux_rec > 0);
reconstruction(I_aux_rec) = 0;
%reconstruction = reconstruction./max(abs(reconstruction(:)));
end

if evalin('base','zef.reconstruction_type') == 5 
aux_rec = rec_x.*n_vec_aux(:,1) + rec_y.*n_vec_aux(:,2) + rec_z.*n_vec_aux(:,3);
I_aux_rec = find(aux_rec <= 0);
reconstruction(I_aux_rec) = 0;
%reconstruction = reconstruction./max(abs(reconstruction(:)));
end

if ismember(evalin('base','zef.reconstruction_type'), [2 3 4 5 7])
reconstruction = smooth_field(reuna_t{i}, reconstruction, size(reuna_p{i}(:,1),1),3);
end

if not(ismember(evalin('base','zef.reconstruction_type'), [6]))
if evalin('base','zef.inv_scale') == 1
reconstruction = 10*log10(max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
elseif evalin('base','zef.inv_scale') == 2
reconstruction = (max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
elseif evalin('base','zef.inv_scale') == 3
reconstruction = sqrt(max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
end
end
end

if evalin('base','zef.use_parcellation')
reconstruction_aux = zeros(size(reconstruction));
for p_ind = selected_list
reconstruction_aux(p_i_ind{p_ind}{2}{ab_ind}) = reconstruction(p_i_ind{p_ind}{2}{ab_ind});
end
reconstruction = reconstruction_aux;
end

axes(h_axes_image); set(h_fig_aux,'visible','on');
h_surf_2{ab_ind} = trisurf(reuna_t{i},reuna_p{i}(:,1),reuna_p{i}(:,2),reuna_p{i}(:,3),reconstruction,'edgecolor','none');
set(h_surf_2{ab_ind},'edgecolor','none','facecolor','flat','facelighting','flat','CDataMapping','scaled');
set(gca,'CLim',[min_rec max_rec]); 
set(h_surf_2{ab_ind},'specularstrength',0.2);
set(h_surf_2{ab_ind},'specularexponent',0.8);
set(h_surf_2{ab_ind},'SpecularColorReflectance',0.8);
set(h_surf_2{ab_ind},'diffusestrength',1);
set(h_surf_2{ab_ind},'ambientstrength',1);
if evalin('base','zef.brain_transparency') < 1
f_alpha_aux = zeros(size(reuna_p{i},1),1);
if evalin('base','zef.inv_scale') == 1
r_alpha_aux = (reconstruction-min(reconstruction))/(max(reconstruction)-min(reconstruction));
else
r_alpha_aux = abs(reconstruction)/max(abs(reconstruction));
end
r_alpha_aux= max(0,r_alpha_aux-min(r_alpha_aux));
r_alpha_aux = r_alpha_aux/max(r_alpha_aux);
f_alpha_aux(reuna_t{i}(:,1)) = f_alpha_aux(reuna_t{i}(:,1)) + r_alpha_aux/3;
f_alpha_aux(reuna_t{i}(:,2)) = f_alpha_aux(reuna_t{i}(:,2)) + r_alpha_aux/3;
f_alpha_aux(reuna_t{i}(:,3)) = f_alpha_aux(reuna_t{i}(:,3)) + r_alpha_aux/3; 
set(h_surf_2{ab_ind},'FaceVertexAlpha',max(evalin('base','zef.brain_transparency'),f_alpha_aux));
set(h_surf_2{ab_ind},'FaceAlpha','interp');
set(h_surf_2{ab_ind},'AlphaDataMapping','none'); 
end

if ismember(i,aux_brain_ind) && cb_done == 0 && ismember(evalin('base','zef.visualization_type'),[3])
cb_done = 1;
h_colorbar = colorbar(h_axes_image,'EastOutside','Position',[0.94 0.675 0.01 0.29],'fontsize',1500);
h_axes_hist = axes('position',[0.03 0.035 0.2 0.1],'visible','off');
h_bar = bar(h_axes_hist,b_hist+(max_rec-min_rec)/(2*50),a_hist,'hist');
set(h_bar,'facecolor',[0.5 0.5 0.5]);
set(h_axes_hist,'ytick',[]);
hist_ylim =  [0 log10(length_reconstruction)];
set(h_axes_hist,'ylim',hist_ylim);
set(h_axes_hist,'fontsize',1500);
set(h_axes_hist,'xlim',[min_rec max_rec]);
set(h_axes_hist,'linewidth',200);
set(h_axes_hist,'ticklength',[0 0]);
axes(h_axes_image); set(h_fig_aux,'visible','on');

h_axes_text = axes('position',[0.03 0.94 0.01 0.05],'visible','off');
h_text = text(0, 0.5, ['Time: ' num2str(evalin('base','zef.inv_time_1') + evalin('base','zef.inv_time_2')/2 + frame_step*(f_ind - 1)*evalin('base','zef.inv_time_3'),'%0.6f') ' s, Frame: ' num2str(f_ind) ' / ' num2str(length_reconstruction_cell) '.']);
  set(h_text,'visible','on','fontsize',1500);
axes(h_axes_image); set(h_fig_aux,'visible','on');
end
%set(h_colorbar,'layer','bottom');
lighting phong;



else
    
h_surf = trimesh(reuna_t{i},reuna_p{i}(:,1),reuna_p{i}(:,2),reuna_p{i}(:,3),'edgecolor','none','facecolor',color_str);
set(h_surf,'specularstrength',0.1);
set(h_surf,'diffusestrength',0.5);
set(h_surf,'ambientstrength',0.85);
set(h_surf,'facealpha',evalin('base','zef.layer_transparency'));
%if not(evalin('base','zef.visualization_type')==3);
lighting flat;
%end

end
end
end
end

view(evalin('base','zef.azimuth'),evalin('base','zef.elevation'));
axis('image');
camva(evalin('base','zef.cam_va'));
if evalin('base','zef.axes_visible')
set(gca,'visible','on');
set(gca,'xGrid','on');
set(gca,'yGrid','on');
set(gca,'zGrid','on');
else
set(gca,'visible','off');
set(gca,'xGrid','off');
set(gca,'yGrid','off');
set(gca,'zGrid','off');
end  

if iscell(evalin('base','zef.reconstruction')) &  evalin('base','zef.visualization_type') == 3

  if is_video  
  if file_index == 1; 
  print(h_fig_aux,'-djpeg95','-r1',[file_path  file_name(1:end-4) '_' int2str(frame_start) file_name(end-3:end)]); 
  elseif file_index ==2; 
  print(h_fig_aux,'-dtiff','-r1',[file_path  file_name(1:end-4) '_' int2str(frame_start) file_name(end-3:end)]); 
  elseif file_index ==3;
  print(h_fig_aux,'-dpng','-r1',[file_path  file_name(1:end-4) '_' int2str(frame_start) file_name(end-3:end)]); 
  elseif file_index ==4;
  bmp_file_temp = [file_path  file_name(1:end-4) '_temp.bmp'];
  print(h_fig_aux,'-dbmp','-r1',bmp_file_temp);
  [movie_frame] = imread(bmp_file_temp,'bmp');
  h_frame = im2frame(movie_frame);
  writeVideo(h_aviobj,h_frame);
   delete(bmp_file_temp);
  end;
  else
  if file_index == 1; 
  print(h_fig_aux,'-djpeg95','-r1',[file_path  file_name]); 
  elseif file_index ==2; 
  print(h_fig_aux,'-dtiff','-r1',[file_path  file_name]); 
  elseif file_index==3; 
  print(h_fig_aux,'-dpng','-r1',[file_path  file_name]); 
  end;
  end;
 
else
camva(c_va);
campos(c_pos);
camtarget(c_ta);
camproj(c_p); 
camup(c_u);
drawnow;
if file_index == 1; 
print(h_fig_aux,[file_path file_name],'-djpeg95','-r1'); 
elseif file_index ==2; 
print(h_fig_aux,[file_path file_name],'-dtiff','-r1'); 
elseif file_index ==3; 
print(h_fig_aux,[file_path file_name],'-dpng','-r1'); 
end;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ismember(evalin('base','zef.visualization_type'),[3])

tic;
f_ind_aux = 1;
for f_ind = frame_start + frame_step : frame_step : frame_stop
f_ind_aux = f_ind_aux + 1;
if f_ind_aux > 2
time_val = toc;
waitbar(f_ind_aux/number_of_frames,h_waitbar,['Frame ' int2str(f_ind_aux) ' of ' int2str(number_of_frames) '. Ready approx: ' datestr(datevec(now+((number_of_frames-1)/(f_ind_aux-2) - 1)*time_val/86400)) '.']); 
else
waitbar(f_ind_aux/number_of_frames,h_waitbar,['Frame ' int2str(f_ind_aux) ' of ' int2str(number_of_frames) '.'])
end
set(h_waitbar,'handlevisibility','off');
%******************************************************
if iscell(evalin('base','zef.reconstruction')) 
reconstruction = evalin('base',['zef.reconstruction{' int2str(f_ind) '}']);
else
reconstruction = evalin('base','zef.reconstruction');
end
reconstruction = reconstruction(:);  
reconstruction = reshape(reconstruction,3,length(reconstruction)/3);
reconstruction = sqrt(sum(reconstruction.^2))';
reconstruction = sum(reconstruction(s_i_ind_2),2)/4;
if evalin('base','zef.inv_scale') == 1
reconstruction = 10*log10(max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
elseif evalin('base','zef.inv_scale') == 2
reconstruction = (max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
elseif evalin('base','zef.inv_scale') == 3
reconstruction = sqrt(max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
end
[a_hist, b_hist] = hist(reconstruction,min_rec:(max_rec-min_rec)/50:max_rec - (max_rec - min_rec)/50);
a_hist = max(0,real(log10(a_hist)));
length_reconstruction = length(reconstruction);
%******************************************************

for i = intersect(aux_brain_ind, aux_brain_visible_ind)
ab_ind = find(aux_brain_ind == i);
reconstruction = single(evalin('base',['zef.reconstruction{' int2str(f_ind) '}']));
reconstruction = reconstruction(:);  
reconstruction = reshape(reconstruction,3,length(reconstruction)/3);

if ismember(evalin('base','zef.reconstruction_type'),[1 7])
reconstruction = sqrt(sum(reconstruction.^2))';
elseif evalin('base','zef.reconstruction_type') == 6
reconstruction = (1/sqrt(3))*sum(reconstruction)';
end
if ismember(evalin('base','zef.reconstruction_type'), [1 6 7])
reconstruction = sum(reconstruction(s_i_ind{ab_ind}),2)/3;
end

if ismember(evalin('base','zef.reconstruction_type'), [2 3 4 5])
rec_x = reconstruction(1,:)';
rec_y = reconstruction(2,:)';
rec_z = reconstruction(3,:)';
rec_x = sum(rec_x(s_i_ind{ab_ind}),2)/3;
rec_y = sum(rec_y(s_i_ind{ab_ind}),2)/3;
rec_z = sum(rec_z(s_i_ind{ab_ind}),2)/3;
n_vec_aux = cross(reuna_p{i}(reuna_t{i}(:,2),:)' - reuna_p{i}(reuna_t{i}(:,1),:)',... 
 reuna_p{i}(reuna_t{i}(:,3),:)' - reuna_p{i}(reuna_t{i}(:,1),:)')';
n_vec_aux = n_vec_aux./repmat(sqrt(sum(n_vec_aux.^2,2)),1,3);
end

if ismember(evalin('base','zef.reconstruction_type'), [2 3 4 5])
reconstruction = sqrt((rec_x.*n_vec_aux(:,1)).^2 + (rec_y.*n_vec_aux(:,2)).^2 + (rec_z.*n_vec_aux(:,3)).^2);
end

if evalin('base','zef.reconstruction_type') == 3
reconstruction = sqrt((rec_x - rec_x.*abs(n_vec_aux(:,1))).^2 + (rec_y - rec_y.*abs(n_vec_aux(:,2))).^2 + (rec_z - rec_z.*abs(n_vec_aux(:,3))).^2);
end

if evalin('base','zef.reconstruction_type') == 4
aux_rec = rec_x.*n_vec_aux(:,1) + rec_y.*n_vec_aux(:,2) + rec_z.*n_vec_aux(:,3);
I_aux_rec = find(aux_rec > 0);
reconstruction(I_aux_rec) = 0;
%reconstruction = reconstruction./max(abs(reconstruction(:)));
end

if evalin('base','zef.reconstruction_type') == 5 
aux_rec = rec_x.*n_vec_aux(:,1) + rec_y.*n_vec_aux(:,2) + rec_z.*n_vec_aux(:,3);
I_aux_rec = find(aux_rec <= 0);
reconstruction(I_aux_rec) = 0;
%reconstruction = reconstruction./max(abs(reconstruction(:)));
end

if ismember(evalin('base','zef.reconstruction_type'), [2 3 4 5 7])
reconstruction = smooth_field(reuna_t{i}, reconstruction, size(reuna_p{i}(:,1),1),3);
end

if not(ismember(evalin('base','zef.reconstruction_type'), [6]))
if evalin('base','zef.inv_scale') == 1
reconstruction = 10*log10(max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));
elseif evalin('base','zef.inv_scale') == 2
reconstruction = (max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
elseif evalin('base','zef.inv_scale') == 3
reconstruction = sqrt(max(reconstruction/max_abs_reconstruction,1/evalin('base','zef.inv_dynamic_range')));    
end
end

if evalin('base','zef.use_parcellation')
reconstruction_aux = zeros(size(reconstruction));
for p_ind = selected_list
reconstruction_aux(p_i_ind{p_ind}{2}{ab_ind}) = reconstruction(p_i_ind{p_ind}{2}{ab_ind});
end
reconstruction = reconstruction_aux;
end

delete(h_surf_2{ab_ind});

axes(h_axes_image); set(h_fig_aux,'visible','on');
h_surf_2{ab_ind} = trisurf(reuna_t{i},reuna_p{i}(:,1),reuna_p{i}(:,2),reuna_p{i}(:,3),reconstruction,'edgecolor','none');
set(h_surf_2{ab_ind},'edgecolor','none','facecolor','flat','facelighting','flat','CDataMapping','scaled');
set(gca,'CLim',[min_rec max_rec]); 
set(h_surf_2{ab_ind},'specularstrength',0.2);
set(h_surf_2{ab_ind},'specularexponent',0.8);
set(h_surf_2{ab_ind},'SpecularColorReflectance',0.8);
set(h_surf_2{ab_ind},'diffusestrength',1);
set(h_surf_2{ab_ind},'ambientstrength',1);
if evalin('base','zef.brain_transparency') < 1
f_alpha_aux = zeros(size(reuna_p{i},1),1);
if evalin('base','zef.inv_scale') == 1
r_alpha_aux = (reconstruction-min(reconstruction))/(max(reconstruction)-min(reconstruction));
else
r_alpha_aux = abs(reconstruction)/max(abs(reconstruction));
end
r_alpha_aux= max(0,r_alpha_aux-min(r_alpha_aux));
r_alpha_aux = r_alpha_aux/max(r_alpha_aux);
f_alpha_aux(reuna_t{i}(:,1)) = f_alpha_aux(reuna_t{i}(:,1)) + r_alpha_aux/3;
f_alpha_aux(reuna_t{i}(:,2)) = f_alpha_aux(reuna_t{i}(:,2)) + r_alpha_aux/3;
f_alpha_aux(reuna_t{i}(:,3)) = f_alpha_aux(reuna_t{i}(:,3)) + r_alpha_aux/3; 
set(h_surf_2{ab_ind},'FaceVertexAlpha',max(evalin('base','zef.brain_transparency'),f_alpha_aux));
set(h_surf_2{ab_ind},'FaceAlpha','interp');
set(h_surf_2{ab_ind},'AlphaDataMapping','none'); 
end
camorbit(frame_step*evalin('base','zef.orbit_1')/15,frame_step*evalin('base','zef.orbit_2')/15);
lighting phong;

end


delete(h_text);
axes(h_axes_text);set(h_fig_aux,'visible','on');
set(h_axes_text,'tag','image_details');
h_text = text(0, 0.5, ['Time: ' num2str(evalin('base','zef.inv_time_1') + evalin('base','zef.inv_time_2')/2 + frame_step*(f_ind - 1)*evalin('base','zef.inv_time_3'),'%0.6f') ' s, Frame: ' num2str(f_ind) ' / ' num2str(length_reconstruction_cell) '.']);
set(h_text,'visible','on','fontsize',1500);
set(h_text,'visible','on');
set(h_axes_text,'layer','bottom');
%drawnow;

delete(h_bar);
h_bar = bar(h_axes_hist,b_hist+(max_rec-min_rec)/(2*50),a_hist,'hist');
set(h_bar,'facecolor',[0.5 0.5 0.5]);
set(h_axes_hist,'ytick',[]);
hist_ylim =  [0 log10(length_reconstruction)];
set(h_axes_hist,'ylim',hist_ylim);
set(h_axes_hist,'fontsize',1500);
set(h_axes_hist,'xlim',[min_rec max_rec]);
set(h_axes_hist,'linewidth',200);
set(h_axes_hist,'ticklength',[0 0]);


  if file_index == 1; 
  print(h_fig_aux,'-djpeg95','-r1',[file_path  file_name(1:end-4) '_' int2str(f_ind) file_name(end-3:end)]); 
  elseif file_index ==2;  
  print(h_fig_aux,'-dtiff','-r1',[file_path  file_name(1:end-4) '_' int2str(f_ind) file_name(end-3:end)]); 
  elseif file_index ==3; 
  print(h_fig_aux,'-dpng','-r1',[file_path  file_name(1:end-4) '_' int2str(f_ind) file_name(end-3:end)]); 
  elseif file_index == 4;
  print(h_fig_aux,'-dbmp','-r1',bmp_file_temp);
  [movie_frame] = imread(bmp_file_temp,'bmp');
  h_frame = im2frame(movie_frame);
  writeVideo(h_aviobj,h_frame);
   delete(bmp_file_temp);
  end;

end

end

    
else
    
    

i = 0;

for k = 1 : 18
switch k
    case 1
        on_val = evalin('base','zef.d1_on');
        visible_val = evalin('base','zef.d1_visible');
        color_str =  evalin('base','zef.d1_color');
    case 2
        on_val = evalin('base','zef.d2_on');
        visible_val = evalin('base','zef.d2_visible');
        color_str = evalin('base','zef.d2_color');
    case 3
        on_val = evalin('base','zef.d3_on');
        visible_val = evalin('base','zef.d3_visible');
        color_str = evalin('base','zef.d3_color');
    case 4
        on_val = evalin('base','zef.d4_on');
        visible_val = evalin('base','zef.d4_visible');
        color_str = evalin('base','zef.d4_color');
            case 5
        on_val = evalin('base','zef.d5_on');
        visible_val = evalin('base','zef.d5_visible');
        color_str =  evalin('base','zef.d5_color');
    case 6
        on_val = evalin('base','zef.d6_on');
        visible_val = evalin('base','zef.d6_visible');
        color_str = evalin('base','zef.d6_color');
    case 7
        on_val = evalin('base','zef.d7_on');
        visible_val = evalin('base','zef.d7_visible');
        color_str = evalin('base','zef.d7_color');
    case 8
        on_val = evalin('base','zef.d8_on');
        visible_val = evalin('base','zef.d8_visible');
        color_str = evalin('base','zef.d8_color');
            case 9
        on_val = evalin('base','zef.d9_on');
        visible_val = evalin('base','zef.d9_visible');
        color_str =  evalin('base','zef.d9_color');
    case 10
        on_val = evalin('base','zef.d10_on');
        visible_val = evalin('base','zef.d10_visible');
        color_str = evalin('base','zef.d10_color');
    case 11
        on_val = evalin('base','zef.d11_on');
        visible_val = evalin('base','zef.d11_visible');
        color_str = evalin('base','zef.d11_color');
    case 12
        on_val = evalin('base','zef.d12_on');
        visible_val = evalin('base','zef.d12_visible');
        color_str = evalin('base','zef.d12_color');
 case 13
        on_val = evalin('base','zef.d13_on');
        visible_val = evalin('base','zef.d13_visible');
        color_str = evalin('base','zef.d13_color');
   case 14
        on_val = evalin('base','zef.w_on');
        visible_val = evalin('base','zef.w_visible');
        color_str = evalin('base','zef.w_color');
    case 15
        on_val = evalin('base','zef.g_on');
        visible_val = evalin('base','zef.g_visible');
        color_str = evalin('base','zef.g_color');
    case 16
        on_val = evalin('base','zef.c_on');
        visible_val = evalin('base','zef.c_visible');
        color_str = evalin('base','zef.c_color');
    case 17
        on_val = evalin('base','zef.sk_on');
        visible_val = evalin('base','zef.sk_visible');
        color_str = evalin('base','zef.sk_color');
    case 18
        on_val = evalin('base','zef.sc_on');
        visible_val = evalin('base','zef.sc_visible');
        color_str = evalin('base','zef.sc_color');
    end
if on_val 
i = i + 1;    
if visible_val
h_surf = trimesh(reuna_t{i},reuna_p{i}(:,1),reuna_p{i}(:,2),reuna_p{i}(:,3),'edgecolor','none','facecolor',color_str);
set(h_surf,'specularstrength',0.1);
set(h_surf,'diffusestrength',0.5);
set(h_surf,'ambientstrength',0.85);
set(h_surf,'facealpha',evalin('base','zef.layer_transparency'));
lighting flat;
end
end
end

view(evalin('base','zef.azimuth'),evalin('base','zef.elevation'));
axis('image');
camva(evalin('base','zef.cam_va'));
if evalin('base','zef.axes_visible')
set(gca,'visible','on');
set(gca,'xGrid','on');
set(gca,'yGrid','on');
set(gca,'zGrid','on');
else
set(gca,'visible','off');
set(gca,'xGrid','off');
set(gca,'yGrid','off');
set(gca,'zGrid','off');
end
drawnow;

camva(c_va);
campos(c_pos);
camtarget(c_ta);
camproj(c_p); 
camup(c_u);

drawnow;

if file_index == 1; 
print(h_fig_aux,[file_path file_name],'-djpeg95','-r1'); 
elseif file_index ==2; 
print(h_fig_aux,[file_path file_name],'-dtiff','-r1'); 
elseif file_index ==3; 
print(h_fig_aux,[file_path file_name],'-dpng','-r1'); 
end;


end


if iscell(evalin('base','zef.reconstruction')) &  evalin('base','zef.visualization_type') == 3 & is_video & file_index == 4
close(h_aviobj);
warning off; delete('avi_file'); warning on
movefile(avi_file_temp, avi_file, 'f');
end
close(h_fig_aux);
if  iscell(evalin('base','zef.reconstruction')) &  evalin('base','zef.visualization_type') == 3    
close(h_waitbar);     
end
%**************************************************************************

end

 


 