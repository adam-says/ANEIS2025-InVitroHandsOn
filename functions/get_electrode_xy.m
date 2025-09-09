%% Helper function to determine the spatial location of the electrode in the MCS MEA system

function [x, y] = get_electrode_xy(elec_nr,metafile)

tmp = ['ID_' num2str(elec_nr-1)];
electrode_name = metafile.chan_names.(tmp);

code = 'ABCDEFGHJKLM';
Lia = ismember(code,electrode_name(1));
out = find(Lia);
if out == 0
    warning('Column not found!')
end

x = out;
y = str2double(extract(electrode_name, digitsPattern));

end