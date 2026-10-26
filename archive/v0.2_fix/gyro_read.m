function d = gyro_read(filename, cfg)
% Read ASCII numeric TXT/CSV, preserve original rows, NO imputation.
if ~isfile(filename)
    error('File not found: %s. Run E00 first or edit cfg.file.', filename);
end
try
    A = readmatrix(filename, 'FileType', 'text');
catch ME
    error('Unable to read file (%s): %s', filename, ME.message);
end
if isempty(A) || size(A,2) < max(cfg.cols)
    error('Data must have >= %d numeric columns. Current: %d.', max(cfg.cols),size(A,2));
end
needed = cfg.cols;
valid = all(isfinite(A(:,needed)),2);
nDropped = sum(~valid);
first=find(valid,1,'first');last=find(valid,1,'last');
if isempty(first),error('No finite a/s data rows found');end
if any(~valid(first:last))
    error('NaN/Inf inside the measurement sequence: possible lost/corrupt samples. Do not silently compress timestamps.');
end
A = A(first:last,:);
if size(A,1) < 600
    error('Too few valid rows (%d); require >=600.',size(A,1));
end
if nDropped > 0
    warning('%d nonnumeric leading/trailing rows removed (not internal samples). Check file headers.',nDropped);
end
y = double(A(:,cfg.cols));
if strcmpi(cfg.biasMode,'manual')
    y = y - reshape(cfg.bias,1,2);
elseif ~strcmpi(cfg.biasMode,'none')
    error('cfg.biasMode must be none or manual');
end
d.A = A;
d.y = y;
d.t = (0:size(A,1)-1)'/cfg.fs;
d.N = size(A,1);
d.source = filename;
d.removedRows = nDropped;
d.fs = cfg.fs;
end
