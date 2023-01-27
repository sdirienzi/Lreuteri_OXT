
%% Using loops 

channel_1_directory = dir("EcadherinChannel");
channel_3_directory = dir("DAPIChannel"); 

channel_3_list_nametemp = {channel_3_directory.name};
channel_3_list_name = channel_3_list_nametemp(1,4:length(channel_3_list_nametemp)); %need to look up the list to remove unused elements
channel_1_list_nametemp = {channel_1_directory.name}; 
channel_1_list_name = channel_1_list_nametemp(1,4:length(channel_1_list_nametemp)); %need to look up the list to remove unused elements

channel_3_list_foldertemp = {channel_3_directory.folder} ;
channel_3_list_folder = channel_3_list_foldertemp(1,4:length(channel_3_list_foldertemp));
channel_1_list_foldertemp = {channel_1_directory.folder} ;
channel_1_list_folder = channel_1_list_foldertemp(1,4:length(channel_1_list_foldertemp));


backslash ={'/'}; %\

channel_1_list_list = strcat(channel_1_list_folder,backslash, channel_1_list_name);

channel_3_list_list = strcat(channel_3_list_folder,backslash, channel_3_list_name);


nucleus_count_by_area_list = []; 
%imagename_list = [];
for i = 1:size(channel_1_list_list, 2)
   % nucleus_count_by_area = countCell_usingArea (channel_3_list_list{i}, channel_1_list_list{i} , 900, 850, 50, 1); %mouse stomach
   % nucleus_count_by_area = countCell_usingArea (channel_3_list_list{i}, channel_1_list_list{i} , 900, 1050, 300, 1); %mouse not stomach
     %  nucleus_count_by_area = countCell_usingArea (channel_3_list_list{i}, channel_1_list_list{i} , 900, 650, 300, 0); %mouse M6 Mid LI
     nucleus_count_by_area = countCell_usingArea (channel_3_list_list{i}, channel_1_list_list{i} , 900, 1050, 350, 0);  %human
    nucleus_count_by_area_list = [nucleus_count_by_area_list;nucleus_count_by_area];
    i 
end 
