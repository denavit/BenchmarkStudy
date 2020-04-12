clear all; close all; clc;

study_path = fullfile(pathOf.BenchmarkStudyResults,'RC');
study = BenchmarkStudy('open',study_path);
load(study.path_of_data,'data')

fid = fopen('data.csv','w');

fprintf(fid,'index,fc,fy,shape,longitudinal_config,longitudinal_bar_size,transverse_config,L,frame_type\n');
for i = 1:length(data)
    fprintf(fid,'%i,',i);
    fprintf(fid,'%g,',data(i).fc);
    fprintf(fid,'%g,',data(i).fy);
    fprintf(fid,'%s,',data(i).shape);
    fprintf(fid,'%s,',data(i).longitudinal_config);
    fprintf(fid,'%s,',data(i).longitudinal_bar_size);
    fprintf(fid,'%s,',data(i).transverse_config);
    fprintf(fid,'%g,',data(i).L);
    fprintf(fid,'%s\n',data(i).frame_type);
end
fclose(fid);


