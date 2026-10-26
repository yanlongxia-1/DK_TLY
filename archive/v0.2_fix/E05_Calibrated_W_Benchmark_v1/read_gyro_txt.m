function record = read_gyro_txt(filePath)
%READ_GYRO_TXT 读取 20200824 数据并识别单整数文件头。
% 不把首行整数当作采样点；保留异常行计数，不静默填值。

fid = fopen(filePath, 'r');
if fid < 0
    error('read_gyro_txt:OpenFailed', '无法打开数据文件：%s', filePath);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
firstLine = fgetl(fid);
if ~ischar(firstLine)
    error('read_gyro_txt:EmptyFile', '空数据文件：%s', filePath);
end
record.headerPresent = ~isempty(regexp(firstLine, '^\s*[+-]?\d+\s*$', 'once'));
if record.headerPresent
    record.headerValue = str2double(strtrim(firstLine));
    raw = readmatrix(filePath, 'FileType', 'text', 'Delimiter', ',', 'NumHeaderLines', 1);
else
    record.headerValue = NaN;
    raw = readmatrix(filePath, 'FileType', 'text', 'Delimiter', ',');
end
if isempty(raw)
    error('read_gyro_txt:EmptyFile', '空数据文件：%s', filePath);
end
if isvector(raw)
    raw = raw(:).';
end

body = raw;

record.bodyRows = size(body, 1);
record.detectedColumnCount = max([0; sum(isfinite(body), 2)]);
if size(body, 2) < 7
    error('read_gyro_txt:TooFewColumns', '文件不足 7 列：%s', filePath);
end

nonblank = any(isfinite(body), 2);
valid = all(isfinite(body(:, 1:7)), 2);
record.invalidRows = sum(nonblank & ~valid);
record.validRows = sum(valid);
record.data = body(valid, 1:7);
record.extraFiniteValues = 0;
if size(body, 2) > 7
    record.extraFiniteValues = sum(sum(isfinite(body(:, 8:end))));
end

record.headerBlockRuleMatches = record.headerPresent && ...
    record.headerValue == ceil(record.validRows / 100);
end

