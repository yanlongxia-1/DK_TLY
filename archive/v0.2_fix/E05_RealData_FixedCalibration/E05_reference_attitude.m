function E05_reference_attitude(thetaGimbal, t, rateTable, cs, ref, out)
% OPTIONAL: attitude comparison, deliberately DISABLED by default.
% Enable only after verifying the gimbal axes/signs/sequence/installation
% and showing that (OmegaX,OmegaY) are body-frame transverse components.
if ~ref.enabled,return;end
if numel(ref.axes)~=3 || any(strcmp(ref.axes,'UNVERIFIED')) || ...
   any(~isfinite(ref.signs)) || ~strcmp(ref.omegaXYFrame,'body') || ...
   ~strcmp(ref.rollBodyAxis,'z')
    error(['E05 angle comparison BLOCKED: the sensor-body mapping is not ',...
           'confirmed. Do not claim physical angle RMSE yet.']);
end
for j=1:3
   if ~ismember(ref.axes{j},{'x','y','z'})
      error('Reference axes must be x/y/z, physically confirmed.');
   end
end
if norm(ref.sensorMount'*ref.sensorMount-eye(3),'fro')>1e-6 || ...
        abs(det(ref.sensorMount)-1)>1e-6
    error('sensorMount must be a proper rotation matrix.');
end
N=numel(t); if N<3,error('Reference comparison needs >= 3 samples.');end
% Gimbal col1/2/3 are ANGLES in degrees; no rate assumption!
% R = R_outer(theta1) * R_middle(theta2) * R_inner(theta3) * mount.
Rref=zeros(3,3,N);
for k=1:N
    R=eye(3);
    for j=1:3
        R=R*axRot(ref.axes{j},deg2rad(ref.signs(j)*thetaGimbal(k,j)));
    end
    Rref(:,:,k)=R*ref.sensorMount;
end
methods={'raw','KF_fixed','FVIM_MVMD_frozen'};
rows={};fig=figure('Visible','off');hold on;
for kk=1:numel(methods)
    nm=methods{kk};
    wx=[nm '_omegaX_deg_s'];wy=[nm '_omegaY_deg_s'];
    if ~ismember(wx,rateTable.Properties.VariableNames),continue;end
    ox=rateTable.(wx);oy=rateTable.(wy);
    wroll=cs.rollSign*360*cs.fRoll;
    Rh=Rref(:,:,1);err=zeros(N,1);
    for k=1:N
        Q=Rh'*Rref(:,:,k);
        err(k)=rad2deg(acos(max(-1,min(1,(trace(Q)-1)/2))));
        if k<N
            w=deg2rad([0.5*(ox(k)+ox(k+1));...
                         0.5*(oy(k)+oy(k+1));wroll]);
            Rh=Rh*rotVec(w*(t(k+1)-t(k)));
        end
    end
    % initial orientation is known from reference for this experiment;
    % this is a CONDITIONAL integration comparison, not full alignment estimation.
    rows(end+1,:)={nm,sqrt(mean(err.^2)),mean(abs(err)),max(err),N}; %#ok<AGROW>
    plot(t,err,'DisplayName',nm);
end
xlabel('Time (s)');ylabel('SO(3) orientation error (deg)');grid on;
legend('Location','best');
saveas(fig,fullfile(out,'E05_orientation_geodesic.png'));close(fig);
T=cell2table(rows,'VariableNames',{'method','RMSE_geodesic_deg',...
    'MAE_geodesic_deg','max_geodesic_deg','nSamples'});
writetable(T,fullfile(out,'E05_orientation_reference_metrics.csv'));
end

function R=axRot(axis,angle)
  switch axis
    case 'x', v=[1;0;0];
    case 'y', v=[0;1;0];
    case 'z', v=[0;0;1];
  end
  R=rotVec(v*angle);
end
function R=rotVec(u)
  a=norm(u);
  if a<1e-12
    W=skew(u); R=eye(3)+W+0.5*W*W;
  else
    W=skew(u/a);
    R=eye(3)+sin(a)*W+(1-cos(a))*W*W;
  end
end
function S=skew(v)
 S=[0,-v(3),v(2);v(3),0,-v(1);-v(2),v(1),0];
end
