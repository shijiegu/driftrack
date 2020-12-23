function [drift,loc]=align_block_1rank(st0,ops)

del_t=1; %isec
data=zeros(round(max(st0(:, 2)))+1,round(max(st0(:, 1)*del_t/ops.fs))+1);
for s=1:length(st0(:, 2))
    l=round(st0(s, 2))+1;
    t=round(st0(s, 1)*del_t/ops.fs)+1;
    data(l,t)=data(l,t)+st0(s, 3);
end

% save data to a location
data_location='/home/lorenlab/Documents/drift_data/run/test.npy';
writeNPY(data,data_location);

%[status, commandOut] = system('conda activate drift');
%commandStr = ['python /home/lorenlab/Documents/affine_warp_drift/affinewarp/drift_warp.py ',data_location];
%[status, commandOut] = system(commandStr);
%if status==0
%    fprintf('Drift Tracking is Done.');
%    drift=readNPY(commandOut);
%end
drift=readNPY('/home/lorenlab/Documents/drift_data/5_15min/shift.npy');
loc=1:round(max(st0(:, 2))+1);
