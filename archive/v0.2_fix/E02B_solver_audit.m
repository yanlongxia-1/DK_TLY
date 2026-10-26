%% E02B_safe_diagnostic.m - resource-bounded solver preflight (not final FVIM optimization)
% Put this file next to gyro_config.m, gyro_read.m, and gyro_mvmd.m in v0.2.
% It does not modify gyro_config.m, original data, or results/E02B.
% Uses ONLY the training interval, never the validation interval.
% Progress is saved to CSV after EACH case so an interruption loses at most one case.
% MATLAB R2019b+ recommended. Run in desktop MATLAB, no parallel pool required.
clearvars; clc;
cfg = gyro_config();

% ----- Conservative limits. The outputs are for numerical diagnosis only. -----
opts.maxSamples    = 1800;                  % 1.8 s; intentionally shorter than full record
opts.K              = [2 4 6];              % avoid expensive K=8/10 in first safe run
opts.alpha          = [200 786.4309 2000];
opts.maxIterations  = 250;                   % hard bound per MVMD call
opts.totalTimeSec   = 360;                   % soft stop between cases (not inside one call)
opts.matlabRAMStopGB = 4.0;                 % stop BEFORE next case if MATLAB too large
opts.minFreeDiskGB   = 8.0;                 % stop BEFORE next case if disk too full
opts.maxThreads      = 2;                   % limit computational pressure if supported

% Prevent accidental use of a parallel pool (this script itself is serial).
try
    prevThreads = maxNumCompThreads;
    maxNumCompThreads(opts.maxThreads);
    fprintf('MATLAB computational threads limited to %d (previous %d).\n', ...
        opts.maxThreads, prevThreads);
catch ME
    fprintf('Thread cap not supported (%s); continuing serial diagnostic.\n',ME.message);
end

% Use only the established, disjoint training interval.
d = gyro_read(cfg.file,cfg);
N=d.N;
i1=max(1,round(N*cfg.trainFraction(1)));
i2=min(N,round(N*cfg.trainFraction(2)));
if i2<i1 || (i2-i1+1)<600
    error('Not enough samples in training interval.');
end
selectedIdx = i1 : min(i2, i1+opts.maxSamples-1);
x = d.y(selectedIdx,:);
if size(x,2)~=2 || any(~isfinite(x(:)))
    error('Expected two finite sensor channels.');
end
[~,stem] = fileparts(cfg.file);
out=fullfile(cfg.outputRoot,stem,'E02B_SAFE');
if ~isfolder(out),mkdir(out);end
csvFile=fullfile(out,'E02B_safe_results.csv');
infoFile=fullfile(out,'E02B_safe_info.txt');

T=table('Size',[0 10], ...
    'VariableTypes',{'double','double','double','double','double', ...
                     'double','double','double','double','string'}, ...
    'VariableNames',{'K','alpha','converged','iterations','delta', ...
                     'centerDelta_Hz','AC_recon_relative','runtime_sec', ...
                     'MATLAB_used_GB','status'});
if isfile(csvFile)
    try
        existing=readtable(csvFile,'TextType','string');
        if isequal(T.Properties.VariableNames,existing.Properties.VariableNames)
            T=existing;
            fprintf('Resume: found %d finished cases in %s\n',height(T),csvFile);
        else
            error('CSV columns mismatch; please move prior file to a backup folder.');
        end
    catch ME
        error('Unable to resume existing CSV safely: %s', ME.message);
    end
end

fid=fopen(infoFile,'w');
if fid<0,error('Unable to write summary to %s',infoFile);end
fprintf(fid,'E02B_SAFE numerical diagnostic, NOT full FVIM or paper validation\n');
fprintf(fid,'Input: %s\n',cfg.file);
fprintf(fid,'Training indices used: %d .. %d (%d samples)\n', ...
    selectedIdx(1),selectedIdx(end),numel(selectedIdx));
fprintf(fid,'fs=%.8g, tau=%.8g, tol=%.8g, maxIterations=%d\n',...
    cfg.fs,cfg.mvmd.tau,cfg.mvmd.tolerance,opts.maxIterations);
fprintf(fid,'K=[%s], alpha=[%s]\n',num2str(opts.K),num2str(opts.alpha));
fprintf(fid,'Soft wall-time limit: %.0fs; MATLAB memory threshold: %.2f GB; free disk threshold: %.2f GB\n',...
    opts.totalTimeSec,opts.matlabRAMStopGB,opts.minFreeDiskGB);
fclose(fid);

totalStart=tic;
stopReason='finished';
for ik=1:numel(opts.K)
    for ia=1:numel(opts.alpha)
        k=opts.K(ik); alpha=opts.alpha(ia);
        if any(T.K==k & abs(T.alpha-alpha)<1e-6)
            fprintf('SKIP completed case K=%d, alpha=%.4f\n',k,alpha);
            continue
        end
        if toc(totalStart)>opts.totalTimeSec
            stopReason='Total time budget reached (safe stop between cases)'; break;
        end
        [usedGB,freeGB]=localResources(out);
        if isfinite(usedGB) && usedGB>opts.matlabRAMStopGB
            stopReason=sprintf('MATLAB RAM usage %.2f GB exceeds %.2f GB',usedGB,opts.matlabRAMStopGB);
            break;
        end
        if isfinite(freeGB) && freeGB<opts.minFreeDiskGB
            stopReason=sprintf('Free disk %.2f GB below %.2f GB',freeGB,opts.minFreeDiskGB);
            break;
        end
        fprintf('\nRUN K=%d alpha=%.4f | %d samples | used RAM=%.2fGB free disk=%.2fGB\n',...
            k,alpha,numel(selectedIdx),usedGB,freeGB);
        caseStart=tic;
        converged=0;iter=NaN;delta=NaN;centerDelta=NaN;acResidual=NaN;
        status='OK';
        try
            [m,cent,infom]=gyro_mvmd(x,cfg.fs,k,alpha,opts.maxIterations,...
                cfg.mvmd.tolerance,cfg.mvmd.tau); %#ok<ASGLU>
            converged=double(infom.converged);
            iter=infom.iterations;
            delta=infom.delta;
            centerDelta=infom.centerDeltaHz;
            acResidual=infom.AC_reconRelative;
            if ~infom.converged
                status='NONCONVERGED';
            end
            clear m cent infom
        catch ME
            status="ERROR: " + string(ME.identifier) + " " + string(ME.message);
            fprintf(2,'Candidate failure: %s\n',char(status));
            clear m cent infom
        end
        elapsed=toc(caseStart);
        [usedAfter,~]=localResources(out);
        row={k,alpha,converged,iter,delta,centerDelta,acResidual,...
             elapsed,usedAfter,string(status)};
        T=[T;cell2table(row,'VariableNames',T.Properties.VariableNames)]; %#ok<AGROW>
        writetable(T,csvFile); % checkpoint immediately after each candidate
        fprintf('SAVED %d/%d | conv=%d delta=%g AC=%.2f%% | %.1fs\n',...
            height(T),numel(opts.K)*numel(opts.alpha),converged,delta,100*acResidual,elapsed);
        if startsWith(string(status),'ERROR')
            stopReason='An MVMD candidate failed; stopped safely after checkpoint';break;
        end
    end
    if ~strcmp(stopReason,'finished'), break; end
end
fid=fopen(infoFile,'a');
if fid>=0
    fprintf(fid,'End status: %s\n',stopReason);
    fprintf(fid,'Number of completed cases: %d of %d\n',height(T),numel(opts.K)*numel(opts.alpha));
    fprintf(fid,'Elapsed total: %.2f seconds\n',toc(totalStart));
    fclose(fid);
end
fprintf('\nE02B SAFE finished: %s\nCSV: %s\n',stopReason,csvFile);
disp(T);
try
    if exist('prevThreads','var') && isnumeric(prevThreads) && isfinite(prevThreads)
        maxNumCompThreads(prevThreads);
    end
catch
end
fprintf('No FVIM multiseed run was attempted. Do NOT use these shortened tests to choose final parameters.\n');

function [matlabUsedGB,freeDiskGB]=localResources(path)
matlabUsedGB=NaN;freeDiskGB=NaN;
try
    u=memory;
    if isstruct(u) && isfield(u,'MemUsedMATLAB')
        matlabUsedGB=double(u.MemUsedMATLAB)/(1024^3);
    end
catch
end
try
    if usejava('jvm')
        checkedPaths={path};
        if ispc
            systemDrive=getenv('SystemDrive');
            if ~isempty(systemDrive)
                checkedPaths{end+1}=[systemDrive filesep]; %#ok<AGROW>
            end
        end
        for ii=1:numel(checkedPaths)
            f=java.io.File(checkedPaths{ii});
            v=double(f.getUsableSpace())/(1024^3);
            if isfinite(v) && v>0
                if isnan(freeDiskGB),freeDiskGB=v;else,freeDiskGB=min(freeDiskGB,v);end
            end
        end
    end
catch
end
end
