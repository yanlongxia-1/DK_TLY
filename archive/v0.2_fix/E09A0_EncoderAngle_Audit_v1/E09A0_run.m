function E09A0_run()
%E09A0_RUN 同一测试区间的外框角度积分评估，复用E08-K已保存角速度。
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'lib'));
cfg=E09A0_config();
assert(isfolder(cfg.dataRoot),'数据目录不存在: %s',cfg.dataRoot);
assert(isfile(cfg.rateCsv),'E08-K角速度序列CSV不存在: %s',cfg.rateCsv);
assert(cfg.referenceAngleCol==1,'未验证其他转台轴角度映射，不允许静默修改');
series=readtable(cfg.rateCsv);
allFiles=[dir(fullfile(cfg.dataRoot,'**','*.TXT'));dir(fullfile(cfg.dataRoot,'**','*.txt')); ...
          dir(fullfile(cfg.dataRoot,'*.TXT'));dir(fullfile(cfg.dataRoot,'*.txt'))];
if isempty(allFiles),error('没有找到原始TXT文件');end
match=arrayfun(@(s) strcmpi(s.name,cfg.encoderTxtName),allFiles);
items=allFiles(match);
if ~isempty(items)
    ps=fullfile({items.folder},{items.name});
    [~,u]=unique(lower(string(ps)),'stable');items=items(u);
end
assert(numel(items)==1,'应找到唯一的%s，实际匹配%d项。请清理重复拷贝或使用独立数据目录。',cfg.encoderTxtName,numel(items));
path=fullfile(items(1).folder,items(1).name);
r=read_gyro_txt(path);
assert(r.invalidRows==0 && r.headerBlockRuleMatches,'原始TXT异常/头行规则不匹配');
assert(size(r.data,2)>=7,'原始TXT列数错误');
[metrics,angleSeries]=E09A0_calculate(r.data(:,cfg.referenceAngleCol),series,cfg);
if ~isfolder(cfg.outputDir),mkdir(cfg.outputDir);end
writetable(metrics,fullfile(cfg.outputDir,'E09A0_encoder_angle_metrics.csv'));
writetable(angleSeries,fullfile(cfg.outputDir,'E09A0_encoder_angle_series.csv'));
save(fullfile(cfg.outputDir,'E09A0_frozen_configuration.mat'),'cfg');
if cfg.makeFigure
    fig=figure('Color','w','Visible','off','Position',[100 100 1200 750]);
    subplot(2,1,1);plot(angleSeries.time_s,angleSeries.encoder_outer_angle_deg,'k--','LineWidth',1.4);hold on;
    l={'outer encoder angle'};vars=angleSeries.Properties.VariableNames;
    for j=1:numel(vars)
        if strncmp(vars{j},'angle_',6)
            plot(angleSeries.time_s,angleSeries.(vars{j}),'LineWidth',.8);
            l{end+1}=vars{j}(7:end); %#ok<AGROW>
        end
    end
    legend(l,'Interpreter','none','Location','best');grid on;
    ylabel('Outer frame angle (deg)');title('E09A0: Integrated outer-frame angle vs encoder');
    subplot(2,1,2);hold on;
    for j=1:numel(vars)
        if strncmp(vars{j},'error_',6)
            plot(angleSeries.time_s,angleSeries.(vars{j}),'LineWidth',.8);
        end
    end
    grid on;xlabel('Time (s)');ylabel('Angular error vs encoder (deg)');
    saveas(fig,fullfile(cfg.outputDir,'E09A0_angle_comparison.png'));close(fig);
end
fid=fopen(fullfile(cfg.outputDir,'E09A0_scope.txt'),'w','n','UTF-8');
assert(fid>0,'无法创建说明文件');
fprintf(fid,'E09-A0：本结果为真实外框编码器角度与积分外框等效角度比较。\n');
fprintf(fid,'不是机体系俯仰/偏航角精度，不用于选择K/alpha/校准值。\n');
fprintf(fid,'数据源TXT: %s\n角速度CSV: %s\n',path,cfg.rateCsv);
fprintf(fid,'参考角度已unwrap; 首样本初始角共同对齐; 禁止事后去漂移或按参考角度反向校准。\n');
fprintf(fid,'下一阶段需转台外/中/内转轴定义、零位/安装矩阵、解调两自由度标定后，才能验证pitch/yaw。\n');
fclose(fid);
disp(metrics);
fprintf('E09A0完成：%s\n',cfg.outputDir);
end
