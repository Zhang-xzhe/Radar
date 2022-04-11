function video_path = choose_video(base_path)
%CHOOSE_VIDEO
%   Allows the user to choose a video (sub-folder in the given path).
%   Returns the full path to the selected sub-folder.
%
%   Jo鉶 F. Henriques, 2012
%   http://www.isr.uc.pt/~henriques/

	%process path to make sure it's uniform
	if ispc(), base_path = strrep(base_path, '\', '/'); end
	if base_path(end) ~= '/', base_path(end+1) = '/'; end
	
	%list all sub-folders
	contents = dir(base_path);%显示base_path目录下的文件和文件夹
	names = {};
	for k = 1:numel(contents),
		name = contents(k).name;
		if isdir([base_path name]) && ~strcmp(name, '.') && ~strcmp(name, '..'),%每个文件夹下都默认含有“.”,“..”两个隐藏的系统文件夹，前者指向该文件夹，后者指向该文件夹的父文件夹
			names{end+1} = name;  %#ok
		end
	end
	
	%no sub-folders found
	if isempty(names), video_path = []; return; end
	
	%choice GUI
	choice = listdlg('ListString',names, 'Name','Choose video', 'SelectionMode','single');%列表选择对话框
	
	if isempty(choice),  %user cancelled
		video_path = [];
	else
		video_path = [base_path names{choice} '/'];
	end
	
end

