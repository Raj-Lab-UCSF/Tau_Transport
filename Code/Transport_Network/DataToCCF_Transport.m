function data_new = DataToCCF_Transport(data_old,regnames,vecormat,loadpath_)

% load([matdir_ filesep 'Mouse_Tauopathy_Data_HigherQ.mat'],'mousedata_struct');
% if isempty(data_old)
%     data_old = mousedata_struct.(studyname).data;
% end
load([loadpath_ filesep 'CCF_labels.mat'],'CCF_labels');
regnames_new = cell(length(regnames),2);
if strcmp(regnames{1}(end),'H')
    strfun_1 = @(x) x(1:end-3);
    regnames_new(:,1) = cellfun(strfun_1,regnames,'UniformOutput',false);
    strfun_2 = @(x) x(end-1:end);
    hems = cellfun(strfun_2,regnames,'UniformOutput',false);
    strfun_3 = @(x) strrep(x,'LH','Left Hemisphere');
    strfun_4 = @(x) strrep(x,'RH','Right Hemisphere');
    hems = cellfun(strfun_3,hems,'UniformOutput',false);
    hems = cellfun(strfun_4,hems,'UniformOutput',false);
    regnames_new(:,2) = hems;
end
rois = zeros(1,length(regnames));
for i = 1:size(regnames_new,1)
    regname = regnames_new{i,1}; hem = regnames_new{i,2};    
    regname_inds = ismember(CCF_labels(:,1),regname);
    hem_inds = ismember(CCF_labels(:,4),hem);
    rois(i) = find((regname_inds + hem_inds) == 2);
end
if vecormat
    data_new = NaN(size(CCF_labels,1),size(data_old,2)); % n ROI in CCF
    data_new(rois,:) = data_old;
else
    data_new = NaN(size(CCF_labels,1),size(CCF_labels,1),size(data_old,3)); % n ROI in CCF
    for i = 1:size(data_old,3)
        data_new(rois,rois,i) = data_old(:,:,i);
    end
end
end