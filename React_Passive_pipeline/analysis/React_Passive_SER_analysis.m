function React_Passive_SER_analysis

Activation_maps
fake_PSTH
raw_signal_analysis
dataset_activation_video
raw_signal_figure

PSTH_figure

% Remove empty voxels
[t1, t2] = find(data == 0);
if t1 ~= 0
    emptyvox = unique(t1);
    fprintf('There are empty voxels. Removing them: %d \n', emptyvox);
    display('Having 0 voxels is not correct, please figure it out');
    data(emptyvox, :) = [];
end

%Check if the data have already been Z-scored at this point.
data = zscore(data,1,2); % zscore in time

% evoked response + decoding (if == 1) + projection on weights + projection on evoked response
decoding = 1;
[index{sess}, datacut{sess}, weights{sess}]  =...
    EvokedResponseAligned(data, 'trigs', decoding, sess, params, datapath, slot_code{n});


end