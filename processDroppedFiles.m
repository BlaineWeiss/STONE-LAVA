function FileData = processDroppedFiles(app, file)
if endsWith(file, '.tif', 'IgnoreCase', true)
    FileData = LoadFile(app, file);
else
    uialert(app.LAVA, 'Only .tif files are supported.', 'File Type Error');
end
end