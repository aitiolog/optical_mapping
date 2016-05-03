function [dataMat] = APD_ROI_mat_conv(dataCell)
%APD_ROI_mat_conv Converts the AP data in cell arrays into matrix
% If the number of APs in the ROIs do not match, it returns the value -1

dataName = inputname(1);

% Convert the cell arrays to matrices for easier manipulation
% Check if the number of APs is the same in all ROIs (-1 if not)

dataCell_size = cellfun(@length, dataCell);
    if all(dataCell_size == dataCell_size(1));
        dataMat = vertcat(dataCell{:});
    else
        disp(['Number of APs in ', dataName,...
            ' ROIs do not match - no summary analysis']);
        %disp(dataName);
        dataMat = -1;
    end

end

