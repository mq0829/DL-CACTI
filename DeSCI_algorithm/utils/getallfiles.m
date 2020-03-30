function [ filenames ] = getallfiles( filedir )
%GETALLFILES get all the file names of a directory
%   filenames=GETALLFILES(filedir) returns all the filenames of the filedir
%   in the form of cell. Each file can be access by filenames{$fileorder$}
%   then, where fileorder is the default order of the file in filedir.
alldir = dir(filedir);
alldir = alldir(~ismember({alldir.name},{'.','..','Thumbs.db'}));
filenames = {alldir.name}';
filenames = strcat({[filedir, '\']}, filenames);
end

