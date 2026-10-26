clear;clc;
cfg=gyro_config();
[~,stem]=fileparts(cfg.file);
src=fullfile(cfg.outputRoot,stem);
if ~isfolder(src),error('No results for %s; run one experiment first.',stem);end
zipFile=fullfile(cfg.root,['RESULTS_',stem,'.zip']);
if isfile(zipFile),delete(zipFile);end
zip(zipFile,{stem},cfg.outputRoot);
fprintf('Exported all completed stages for this data file:\n%s\n',zipFile);
