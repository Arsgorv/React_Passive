
function React_Passive_AG
%(slot_code, slice, comp, preprocessing)
%{
%
% Data created here can be used for figures in other script ALLFIGs.m
% I) Preprocessing all sessions of a slice (align and CCA)
% II) Run analyses and transform data for main figures (decod and PCA)
% III) Analyses by sound type ?
% IV) Control Analyses (outside mask, shuffles, decod-PCA, regressions) ?
% V) Additional Analyses (dPCA, Varproj inversed, Hippo, oscillations)
%
% INPUT
%
%
% OUTPUT
%
%
% EXAMPLE
%
%
% Coded by Arsenii Goriachenkov, ENS - PSL, LSP, Paris, 2022-2026
% github.com/arsgorv

At the moment (25/08/2023) I:

0. Recompose the acquisition data to raw data (.scan -> raw_data)
1. Realign data within the session with normcorr (raw_data -> data_cat)
2. NOT USING IT: I use realignment across day sessions (data_cat -> data_aligned). 
3. I CCA (data_cat -> data_CCA). If no realignment across sessions, I apply different masks for each session
4. I cut into trials (data_CCA -> data_cut_in_trials)
5. I calculate dCBV (data_cut_in_trials -> dCBV)
6. I split on categories (dCBV -> d1 & d2)


How to improve the pipeline: 
- remove all doubling and unnecessary variables
- rename variables so they make sense
- implement a local variable working_data, which you can easily access in
    the beginning of the script to select the data you work on
- put flags to preprocessing and every step of preprocessing
- put flags to all control figures
- make a cooking scripts which has links to all subscripts. You must be able to select a session(s) and to study it selectively

Notes:
- Be very careful with night sessions as after Repositioning, there is a region full of zeros, which must be taken into account.
%}

%% Select sessions
% Form the list of sessions
selection = 3;

Dir{1} = PathForExperimentsReactPassive('Chabichou', 'all', 'all', 'all');
Dir{2} = PathForExperimentsReactPassive('Edel', 'all', 'all', 'all');
Dir{3} = PathForExperimentsReactPassive('Kosichka', 'all', 'all', 'all');

sessions = Dir{selection}.path';

%% PreProcessing: Sleep Scoring
Master_SleepScoring_preproc(sessions)

%% PreProcessing: fUS
Master_fUS_preproc(sessions)

%% PreProcessing: align all datastreams -> extract trial information -> reconstruct missed triggers -> make Epochs
Master_data_sync_preproc(sessions, true)

%% PreProcessing: Behaviour (pupil ; baphy)
% React_Passive_epoch_preprocessing
% Master_behaviour_preproc(sessions)





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% Tonotopy

%% Sound-evoked response analysis
React_Passive_SER_analysis

%% Reactivation analysis
React_Passive_reactivations_analysis

%% Decoding analysis
% Decode sound identity from SER and React

%% Cross-correlation analysis
React_Passive_cross_corr_analysis

%% Arousal analysis
baseline_drift

%% Vasomotion analysis
React_Passive_vasomotion_analysis

%% Dataset figures

%% Paper figures



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% NOT USED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below needs to be sorted later as well as other Florian codes
%
% convert_data_to_A_Bergel_format

% Permanent parameters
% Switch between computers
% switch comp
%     case 'biggerguy'
%         prefix = '/mnt/working1/';
%         
%     case 'icelos'
%         prefix = '/home/';
% end
% animal = {'Edel'; 'Chabichou'};
% % Path to the unprocessed Edel data
% datapath = [prefix 'arsenii/data5/Arsenii/React_Passive_AG/fUS/Processed_Data/' animal{1}];
% dailypath = ['/' slot_code '/CCA/'];
% 
% addpath(genpath([datapath '/' slot_code]))
% cd([datapath '/' slot_code]);
%
% slices = {'05'};
% slot_code = {'O'};
% 
% for n = 1:size(slices,2)
%     dailypath = ['/' slot_code{n} '/CCA/'];
%     
%     load([datapath '/' slot_code{n} '/REA_AL_B_' slices{n} '.mat'], 'params')
%     
%     slice = slices{n};
%     filenames = [datapath dailypath 'REA_' slice '*.mat'];
%     
%     listings = dir(filenames);
%     
%     for sess = 1:size(listings,1)
%         load([datapath dailypath listings(sess).name], 'slice', 'pair', 'data')
%         disp(listings(sess).name)
%         
        
%     end
%     save(['REA_AL_B_' slice ], 'index', 'datacut', 'weights' ,'-append')
%     clear index datacut weights
% end
% 

% ---- Moved it to Figs_AG ---- III) Concatenate data and store it in 'compo' groupped by session.

%
% clear accuracy_all
% for expo = 1:2
% %     filenames = [datapath '/' slot_code '/CCA' '/REA_*_*'  params.session{expo} '.mat'];
%     filenames = [datapath '/DATA_Aligned_CCAed' '/REA_*_p*_*_' num2str(expo) '.mat'];
%     listings = dir(filenames);
%     clear mat_res_all; clear Var_proj_all; clear Ref;
%     for  i = 1:size(listings,1)
%         load([datapath '/DATA_Aligned_CCAed/' listings(i).name], 'VarProj', 'mat_res', 'accuracy') % how to compose these files?
%
%         %regroup all variance of projection in cell
%         Var_proj_all.pre{i} = VarProj(1,:,:);
%         Var_proj_all.post{i} = VarProj(2,:,:);
% %         accuracy_all{i} = nanmean(accuracy,2);
%         accuracy_all(expo, :, i) = nanmean(accuracy, 2);
%         %regroup all components in cell
%         mat_res_all{i} = mat_res;
%         %index
%         Ref{i} = listings(i).name;
%     end
%     compo{expo} = vertcat(mat_res_all{:});
%     VarAll{expo} = Var_proj_all;
% end

% a = 1;

% Perfusion analysis to control that reactivations are not due to an increase in the blood volume
% Load the preprocessed data
% Apply a hippocampus mask
% Average the CBV in each phase and store it in 'perfu'
% If you want to apply this analysis to the hippocampus, unmask the
% data by using NewResp instead of Din

% for n = 1:size(slice,2)
%     slice = slice{n};
%     load(['REA_AL_B_' slice '.mat'], 'NewResp', 'hippo', 'params');
%     
%     n_tbins = size(NewResp,4);
%     n_sess = size(NewResp,3);
%     Mask = hippo;
%     param.ManualMask = Mask;
%     Din = Mat2Pixs(NewResp.*repmat(Mask, [1 1 n_sess n_tbins]), param);
%     
%     perfu(n,1:n_sess,:,1) = [snm(Din(:,:,1:3250), [1 3])' snm(Din(:,:,3251:6500), [1 3])' snm(Din(:,:,6501:9750), [1 3])'];
% end
% save('perfuH.mat', 'perfu')

% ADDITIONAL exploratory analysis
% % Do analyses on sessions on the same slice (dPCA and RS with next)
% for n = 1:size(slice,2)
%     slice=slice{n};
%     load(['REA_AL_B_' slice '.mat'], 'NewResp');
%     NewData = NewRespC;
%
%     R  = ReactivationAnaAligned7(NewData, slice);
%     % dpca
%    [option, NewDataDem, dPCs]  = dPCAna(slice);
%    save(['REA_AL_B_' slice ], 'option', 'dPCs', '-append')
% end

% No idea what is this #1. Looks like part 2 + something else

% % Evoked Response Analysis
% [index, datacut, weights] = EvokedResponse(data, 'trigs', slice, pair, repet);
% [index, datacut, weights] = EvokedResponseNorma(data, 'trigs', slice, pair, repet);
% %% Reactivation Strength Analysis
% [Rcat, phi, lMax, mat_res, PCs]  = ReactivationAna(data, slice, pair, repet);
%
% % close all
% %% RS timelock analysis
% [filt_dat, RS_lock] = fUSFreq(slice, pair, repet);
% %% Reactivation hippo
% slice    = '08'
% pair     = 'p14'
% repet    = 1
% load([datapath '\REA_' slice '_' pair '_' num2str(repet) '.mat'], 'D_c', 'Ct');
% n_tbins = size(D_c, 3)
% param.ManualMask = (Ct == 'hippo');
% Mask=param.ManualMask;
% data = Mat2Pixs(D_c.*repmat(Mask,[1 1 n_tbins]),param);
%
% [index, datacut, weights] = EvokedResponseH(data, 'trigs', slice, pair, repet);
% [Rcat, phi, lMax, mat_res, PCs]  = ReactivationAnaH(data, slice, pair, repet);
% close all
% %% Redo 1  analysis
% filenames = [datapath '\REA_*.mat'];
% listings = dir(filenames) ;
%
% for  i = 1:size(listings,1)
%     load([datapath '\' listings(i).name], 'slice', 'pair', 'repet', 'data')
%     [index, datacut, weights] = EvokedResponseNorma(data, 'trigs', slice, pair, repet);
% end
% %% Correlation hippo analysis
% filenames = [datapath '\REA_*p12*.mat'];
% listings = dir(filenames) ;
%
% for  i = 1:size(listings,1)
%     load([datapath '\' listings(i).name], 'slice', 'pair', 'repet')
%     mean_corr  = CorrHippo(slice, pair, repet)
% end
% %% fig corr hippo all
% filenames = [datapath '\REA_*_1.mat'];
% listings = dir(filenames) ;
% clear mean_corr_all
% for  i = 1:size(listings,1)
%     load([datapath '\' listings(i).name], 'mean_corr')
%     mean_corr_all(i, :,:) = mean_corr;
% end
%
% %
% f1 = figure
%
% plot(-10:1:10,snm(mean_corr_all(:,1,:),1), 'Linewidth', 1.5);
% hold on;
% plot(-10:1:10,snm(mean_corr_all(:,2,:),1), 'Linewidth', 1.5);
% plot(-10:1:10,snm(mean_corr_all(:,3,:),1), 'Linewidth', 1.5);
% legend ('pre',  'post','snd')
% xlabel('Time offset with hippo ROI (s)')
% ylabel ('Correlation Coefficient')
% exportgraphics(f1, ['H_CorrHippo_mean_all_1.jpg'], 'Resolution', '600')
%
% %
% f1 = figure
% plot(-10:1:10,snm(abs(mean_corr_all(:,1,:)),1), 'Linewidth', 1.5);
% hold on;
% plot(-10:1:10,snm(abs(mean_corr_all(:,2,:)),1), 'Linewidth', 1.5);
% plot(-10:1:10,snm(abs(mean_corr_all(:,3,:)),1), 'Linewidth', 1.5);
% legend ('pre',  'post','snd')
% xlabel('Time offset with hippo ROI (s)')
% ylabel ('Correlation Coefficient')
% exportgraphics(f1, ['H_CorrHippo_abs_mean_all_1.jpg'], 'Resolution', '600')
%
% %% look at variance and decoder weights look
% filenames = [datapath '\REA_*_1.mat'];
% listings = dir(filenames) ;
% clear data_all
% clear weights_all
% clear Acx_all
% for  i = 1:size(listings,1)
%     load([datapath '\' listings(i).name], 'data', 'weights', 'Acx')
%     data_all{i}=data;
%     weights_all{i} = weights.coeffs;
%     Acx_all{i}= Acx;
% end
% %%
% for  i = 1:size(listings,1)
%     mean_all(i,1) = mean(var(data_all{i}(:,1:3000)'))
%     mean_all(i,2) = mean(var(data_all{i}(:,6001:9000)'))
%     mean_all(i,3) = mean(var(data_all{i}(:,3001:6000)'))
% end
% plot(mean_all)
% %%
% i=2
% Mask = Acx_all{i};
% param.ManualMask = Mask;
% vox_var1 = Pixs2Mat(var(data_all{i}(:,1:3000)')',param);
% vox_var2 = Pixs2Mat(var(data_all{i}(:,6001:9000)')',param);
% weights = Pixs2Mat(weights_all{i}',param);
% %%
% tp = 15
% %background
% imAlpha=ones(size(vox_var1(5:60,30:120)));
% imAlpha(isnan(vox_var1(5:60,30:120)))=0;
%
% figure
% subplot(3,1,1)
% imagesc(vox_var1(5:60,30:120),'AlphaData',imAlpha); axis equal tight;
% colorbar; caxis([0 30]);% colormap('turbo')
% title('variance pre')
%
% subplot(3,1,2)
% imagesc(vox_var2(5:60,30:120),'AlphaData',imAlpha); axis equal tight;
% colorbar; caxis([0 30]); colormap('autumn')
% title('variance post')
%
% subplot(3,1,3)
% imagesc(weights(5:60,30:120,tp)/snm(weights(:,:,tp),[1 2]),'AlphaData',imAlpha);
% axis equal tight; colorbar; caxis([-1 1]*8)%*1e10);% colormap('hot')
% title(['decoder weights ' num2str(tp*.4) 's'])
%
% %% Realign across session %% WIP
% slice='04'
% data = Realign(slice)
%
% %% add data to big matrix\structure
% % 4 repet// 2*2sons
%
% Data_All(slice, repet).D_c(repet)
% Data_All(slice, repet).data_rea(repet)
% Data_All(slice, repet).datacut(repet)
%
% %% Figure with all data
%
% filenames = [datapath '\REA_*_1.mat'];
% listings = dir(filenames) ;
% clear mat_res_all
% clear Var_proj_all
% for  i = 1:size(listings,1)
%     load([datapath '\' listings(i).name], 'VarProj', 'mat_res', 'phi')
%     Var_proj_all.pre{i} = VarProj(1,:)
%     Var_proj_all.post{i} = VarProj(2,:)
%     mat_res(:,6)=phi
%     mat_res_all{i} = mat_res
%     Ref{i} = listings(i).name
% end
%
% %%
% % RSFig(vertcat(mat_res_all{:});, slice, pair, repet, [1.5 3 10])
% % RSFig(vertcat(mat_res_all{:}), 'all', 'all', 0, [1.5 3 10])
% RSFig(vertcat(mat_res_all{:}), '_1', '_', 0, [1.2 1.6 3 ])
%
% %%
% %% Project pre and post on decodeur max point
% maxx=28
% for i = 1:size(Var_proj_all.pre,2)
%     bis(1,i,:) = Var_proj_all.pre{i}(1:maxx)...
%         / snm( [Var_proj_all.pre{i}(1:maxx) Var_proj_all.post{i}(1:maxx)],2)
%     bis(2,i,:) = Var_proj_all.post{i}(1:maxx)...
%         / snm( [Var_proj_all.pre{i}(1:maxx) Var_proj_all.post{i}(1:maxx)],2)
% end
%
% f1 = figure
% patch([5 7 7 5], [0 0 100*[1 1]], [.91 .91 .91], 'EdgeColor','none', 'HandleVisibility', 'off')
% hold on
% patch([5 7 7 5], [0 0 -100*[1 1]], [.91 .91 .91], 'EdgeColor','none', 'HandleVisibility', 'off')
%
% xl = (1:maxx)/2.5;
% [p1, p2] = semshade(squeeze(bis(1,:,:)), .5,[0, 0.4470, 0.7410], xl);
% p2.HandleVisibility = 'off'
% hold on
%
% [p3, p4] = semshade(squeeze(bis(2,:,:)), .5,[0.8500, 0.3250, 0.0980], xl);
% p4.HandleVisibility = 'off'
% legend ({'Pre', 'Post'})
% title('Variance of period projected on decoder weights from each timepoint')
% lims = minmax([snm(bis(1,:,:),2) snm(bis(2,:,:),2)]);
% ylim([min(lims(:,1))-.25*min(lims(:,1)), max(lims(:,2))+ .25*max(lims(:,2))])
% xlabel('Time (s)')
% ylabel('Variance (Normalized by session)')
% xlim(minmax(xl))
%
% exportgraphics(f1, ['VarianceOfProjection_1' '' '_' '' '_' num2str(0) '.jpg'], 'Resolution', '600')
%
% No idea what is this #2. Maybe part 4?

%     filename = [datapath '\' listings(i).name]
%     AllPat  = EvokedResponseSpatialv4(filename)
%     HippoAna(slice, pair, repet)
%     ControlAna(slice, pair, repet)
%     type = 'ds'
%     [Rcat, phi, lMax, mat_res, PCs]  = ReactivationAnaSpe(data, datacut, slice, pair, repet, type);
% %     mean_corr  = CorrHippo(slice, pair, repet)
% %     repets(i)=repet;
%     perfu(i,:) = [ snm(data(:,1:3000), [1 2]) snm(data(:,3001:6000), [1 2]) snm(data(:,6001:9000), [1 2])];
%     VarProji = ProjRest(slice, pair, repet)
%     AreaAna(slice, pair, repet)
%     VarControl = EvokedResponseControlVarProj(filename)
%     VarControlAll(i,:) = VarControl

end