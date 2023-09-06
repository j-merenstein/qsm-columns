%% Load all data %%
cd('/Volumes/Data/Madden/QSM.02/Analysis/columns/');

load('QSM_results/NoFilter_3CurvBins_21points_threshold_T2_36ROIs/QSM_all_0.5_6mm_T2_Lausanne_neg.mat');
QSM_ROIs_neg = QSM_ROIs_all;

load('QSM_results/NoFilter_3CurvBins_21points_threshold_T2_36ROIs/QSM_all_0.5_6mm_T2_Lausanne_pos.mat');
QSM_ROIs_pos = QSM_ROIs_all;

% Flip: so that GM/WM boundary is on right side (i.e., column 21)
QSM_neg = flip(QSM_ROIs_neg,2);
QSM_pos = flip(QSM_ROIs_pos,2);


%% Run t-tests between groups %%
% Positive QSM
for k1 = 1:size(QSM_pos,1) % Number of ROIs
    for k2 = 1:21 % Number of cortical depths
        for k3 = 1:4 % Gyral crown, sulcal bank, sulcal fundus, or combined
            for k4 = 1:3 % LH, RH, bilateral
                [h,p] = ttest2(squeeze(QSM_pos(k1,k2,k3,k4,1:22)),squeeze(QSM_pos(k1,k2,k3,k4,23:44))); % 1-22 = AD, 23-44 = controls
                ps_pos_qsm(k1,k2,k3,k4) = p; 
            end
        end
    end
end

clear k1 k2 k3 k4 n p

% Negative QSM
for k1 = 1:size(QSM_neg,1) % Number of ROIs
    for k2 = 1:21 % Number of depths
        for k3 = 1:4 % Crown, bank, fundus, or combined
            for k4 = 1:3 % LH, RH, bilateral
                [h,p] = ttest2(squeeze(QSM_neg(k1,k2,k3,k4,1:22)),squeeze(QSM_neg(k1,k2,k3,k4,23:44)));
                ps_neg_qsm(k1,k2,k3,k4) = p; 
            end
        end
    end
end

clear k1 k2 k3 k4 n p


%% FDR correction for 21 depths %%
% Positive QSM
for k = 1:34
    [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(squeeze(ps_pos_qsm(k,1:21,1,1)), 0.05, 'pdep', 'no');
    fdr_lh_crown_pos(k,:) = adj_p; 
end

for k = 1:34
    [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(squeeze(ps_pos_qsm(k,1:21,1,2)), 0.05, 'pdep', 'no');
    fdr_rh_crown_pos(k,:) = adj_p; 
end

for k = 1:34
    [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(squeeze(ps_pos_qsm(k,1:21,2,1)), 0.05, 'pdep', 'no');
    fdr_lh_bank_pos(k,:) = adj_p; 
end

for k = 1:34
    [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(squeeze(ps_pos_qsm(k,1:21,2,2)), 0.05, 'pdep', 'no');
    fdr_rh_bank_pos(k,:) = adj_p; 
end

for k = 1:34
    [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(squeeze(ps_pos_qsm(k,1:21,3,1)), 0.05, 'pdep', 'no');
    fdr_lh_fundus_pos(k,:) = adj_p; 
end

for k = 1:34
    [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(squeeze(ps_pos_qsm(k,1:21,3,2)), 0.05, 'pdep', 'no');
    fdr_rh_fundus_pos(k,:) = adj_p; 
end


% Negative QSM
for k = 1:34
    [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(squeeze(ps_neg_qsm(k,1:21,1,1)), 0.05, 'pdep', 'no');
    fdr_lh_crown_neg(k,:) = adj_p; 
end

for k = 1:34
    [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(squeeze(ps_neg_qsm(k,1:21,1,2)), 0.05, 'pdep', 'no');
    fdr_rh_crown_neg(k,:) = adj_p; 
end

for k = 1:34
    [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(squeeze(ps_neg_qsm(k,1:21,2,1)), 0.05, 'pdep', 'no');
    fdr_lh_bank_neg(k,:) = adj_p; 
end

for k = 1:34
    [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(squeeze(ps_neg_qsm(k,1:21,2,2)), 0.05, 'pdep', 'no');
    fdr_rh_bank_neg(k,:) = adj_p; 
end

for k = 1:34
    [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(squeeze(ps_neg_qsm(k,1:21,3,1)), 0.05, 'pdep', 'no');
    fdr_lh_fundus_neg(k,:) = adj_p; 
end

for k = 1:34
    [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(squeeze(ps_neg_qsm(k,1:21,3,2)), 0.05, 'pdep', 'no');
    fdr_rh_fundus_neg(k,:) = adj_p; 
end


%% QSM by depth, across curvature bins (Figure 2) %%
% Positive QSM
for k1 = 1:34 % ROIs
    for k2 = 1:44 % Participants
        for k3 = 1:21 % Depths
            for k4 = 1:2 % Hemisphere (1 = LH, 2 = RH)
                m = mean(QSM_ROIs_pos(k1,k3,:,k4,k2),3); % Takes the mean across the third dimension (curvature)
                mean_pos_column_qsm(k1,k3,k4,k2) = m; 
            end
        end
    end
end

clear k1 k2 k3 k4 m

lh_col1 = mean(squeeze(mean_pos_column_qsm(:,1,1,:)),1,'omitnan')'; % Takes the mean across first dimension (ROI)
lh_col2 = mean(squeeze(mean_pos_column_qsm(:,2,1,:)),1,'omitnan')';
lh_col3 = mean(squeeze(mean_pos_column_qsm(:,3,1,:)),1,'omitnan')';
lh_col4 = mean(squeeze(mean_pos_column_qsm(:,4,1,:)),1,'omitnan')';
lh_col5 = mean(squeeze(mean_pos_column_qsm(:,5,1,:)),1,'omitnan')';
lh_col6 = mean(squeeze(mean_pos_column_qsm(:,6,1,:)),1,'omitnan')';
lh_col7 = mean(squeeze(mean_pos_column_qsm(:,7,1,:)),1,'omitnan')';
lh_col8 = mean(squeeze(mean_pos_column_qsm(:,8,1,:)),1,'omitnan')';
lh_col9 = mean(squeeze(mean_pos_column_qsm(:,9,1,:)),1,'omitnan')';
lh_col10 = mean(squeeze(mean_pos_column_qsm(:,10,1,:)),1,'omitnan')';
lh_col11 = mean(squeeze(mean_pos_column_qsm(:,11,1,:)),1,'omitnan')';
lh_col12 = mean(squeeze(mean_pos_column_qsm(:,12,1,:)),1,'omitnan')';
lh_col13 = mean(squeeze(mean_pos_column_qsm(:,13,1,:)),1,'omitnan')';
lh_col14 = mean(squeeze(mean_pos_column_qsm(:,14,1,:)),1,'omitnan')';
lh_col15 = mean(squeeze(mean_pos_column_qsm(:,15,1,:)),1,'omitnan')';
lh_col16 = mean(squeeze(mean_pos_column_qsm(:,16,1,:)),1,'omitnan')';
lh_col17 = mean(squeeze(mean_pos_column_qsm(:,17,1,:)),1,'omitnan')';
lh_col18 = mean(squeeze(mean_pos_column_qsm(:,18,1,:)),1,'omitnan')';
lh_col19 = mean(squeeze(mean_pos_column_qsm(:,19,1,:)),1,'omitnan')';
lh_col20 = mean(squeeze(mean_pos_column_qsm(:,20,1,:)),1,'omitnan')';
lh_col21 = mean(squeeze(mean_pos_column_qsm(:,21,1,:)),1,'omitnan')';

col_lh_pos = table(lh_col1, lh_col2, lh_col3, lh_col4, lh_col5, lh_col6, lh_col7, lh_col8, lh_col9, lh_col10, lh_col11, lh_col12, lh_col13, lh_col14, lh_col15, lh_col6, lh_col7, lh_col18, lh_col19, lh_col20, lh_col21);
clear lh_*

rh_col1 = mean(squeeze(mean_pos_column_qsm(:,1,2,:)),1,'omitnan')'; % Takes the mean across first dimension (ROI)
rh_col2 = mean(squeeze(mean_pos_column_qsm(:,2,2,:)),1,'omitnan')';
rh_col3 = mean(squeeze(mean_pos_column_qsm(:,3,2,:)),1,'omitnan')';
rh_col4 = mean(squeeze(mean_pos_column_qsm(:,4,2,:)),1,'omitnan')';
rh_col5 = mean(squeeze(mean_pos_column_qsm(:,5,2,:)),1,'omitnan')';
rh_col6 = mean(squeeze(mean_pos_column_qsm(:,6,2,:)),1,'omitnan')';
rh_col7 = mean(squeeze(mean_pos_column_qsm(:,7,2,:)),1,'omitnan')';
rh_col8 = mean(squeeze(mean_pos_column_qsm(:,8,2,:)),1,'omitnan')';
rh_col9 = mean(squeeze(mean_pos_column_qsm(:,9,2,:)),1,'omitnan')';
rh_col10 = mean(squeeze(mean_pos_column_qsm(:,10,2,:)),1,'omitnan')';
rh_col11 = mean(squeeze(mean_pos_column_qsm(:,11,2,:)),1,'omitnan')';
rh_col12 = mean(squeeze(mean_pos_column_qsm(:,12,2,:)),1,'omitnan')';
rh_col13 = mean(squeeze(mean_pos_column_qsm(:,13,2,:)),1,'omitnan')';
rh_col14 = mean(squeeze(mean_pos_column_qsm(:,14,2,:)),1,'omitnan')';
rh_col15 = mean(squeeze(mean_pos_column_qsm(:,15,2,:)),1,'omitnan')';
rh_col16 = mean(squeeze(mean_pos_column_qsm(:,16,2,:)),1,'omitnan')';
rh_col17 = mean(squeeze(mean_pos_column_qsm(:,17,2,:)),1,'omitnan')';
rh_col18 = mean(squeeze(mean_pos_column_qsm(:,18,2,:)),1,'omitnan')';
rh_col19 = mean(squeeze(mean_pos_column_qsm(:,19,2,:)),1,'omitnan')';
rh_col20 = mean(squeeze(mean_pos_column_qsm(:,20,2,:)),1,'omitnan')';
rh_col21 = mean(squeeze(mean_pos_column_qsm(:,21,2,:)),1,'omitnan')';

col_rh_pos = table(rh_col1, rh_col2, rh_col3, rh_col4, rh_col5, rh_col6, rh_col7, rh_col8, rh_col9, rh_col10, rh_col11, rh_col12, rh_col13, rh_col14, rh_col15, rh_col6, rh_col7, rh_col18, rh_col19, rh_col20, rh_col21);
clear rh_*


% Negative QSM
for k1 = 1:34 % ROIs
    for k2 = 1:44 % Participants
        for k3 = 1:21 % Depths
            for k4 = 1:2 % Hemisphere (1 = LH, 2 = RH)
                m = mean(QSM_ROIs_neg(k1,k3,:,k4,k2),3); % Takes the mean across the third dimension (curvature)
                mean_neg_column_qsm(k1,k3,k4,k2) = m; 
            end
        end
    end
end

clear k1 k2 k3 k4 m

lh_col1 = mean(squeeze(mean_neg_column_qsm(:,1,1,:)),1,'omitnan')'; % Takes the mean across first dimension (ROI)
lh_col2 = mean(squeeze(mean_neg_column_qsm(:,2,1,:)),1,'omitnan')';
lh_col3 = mean(squeeze(mean_neg_column_qsm(:,3,1,:)),1,'omitnan')';
lh_col4 = mean(squeeze(mean_neg_column_qsm(:,4,1,:)),1,'omitnan')';
lh_col5 = mean(squeeze(mean_neg_column_qsm(:,5,1,:)),1,'omitnan')';
lh_col6 = mean(squeeze(mean_neg_column_qsm(:,6,1,:)),1,'omitnan')';
lh_col7 = mean(squeeze(mean_neg_column_qsm(:,7,1,:)),1,'omitnan')';
lh_col8 = mean(squeeze(mean_neg_column_qsm(:,8,1,:)),1,'omitnan')';
lh_col9 = mean(squeeze(mean_neg_column_qsm(:,9,1,:)),1,'omitnan')';
lh_col10 = mean(squeeze(mean_neg_column_qsm(:,10,1,:)),1,'omitnan')';
lh_col11 = mean(squeeze(mean_neg_column_qsm(:,11,1,:)),1,'omitnan')';
lh_col12 = mean(squeeze(mean_neg_column_qsm(:,12,1,:)),1,'omitnan')';
lh_col13 = mean(squeeze(mean_neg_column_qsm(:,13,1,:)),1,'omitnan')';
lh_col14 = mean(squeeze(mean_neg_column_qsm(:,14,1,:)),1,'omitnan')';
lh_col15 = mean(squeeze(mean_neg_column_qsm(:,15,1,:)),1,'omitnan')';
lh_col16 = mean(squeeze(mean_neg_column_qsm(:,16,1,:)),1,'omitnan')';
lh_col17 = mean(squeeze(mean_neg_column_qsm(:,17,1,:)),1,'omitnan')';
lh_col18 = mean(squeeze(mean_neg_column_qsm(:,18,1,:)),1,'omitnan')';
lh_col19 = mean(squeeze(mean_neg_column_qsm(:,19,1,:)),1,'omitnan')';
lh_col20 = mean(squeeze(mean_neg_column_qsm(:,20,1,:)),1,'omitnan')';
lh_col21 = mean(squeeze(mean_neg_column_qsm(:,21,1,:)),1,'omitnan')';

col_lh_neg = table(lh_col1, lh_col2, lh_col3, lh_col4, lh_col5, lh_col6, lh_col7, lh_col8, lh_col9, lh_col10, lh_col11, lh_col12, lh_col13, lh_col14, lh_col15, lh_col6, lh_col7, lh_col18, lh_col19, lh_col20, lh_col21);
clear lh_*

rh_col1 = mean(squeeze(mean_neg_column_qsm(:,1,2,:)),1,'omitnan')'; % Takes the mean across first dimension (ROI)
rh_col2 = mean(squeeze(mean_neg_column_qsm(:,2,2,:)),1,'omitnan')';
rh_col3 = mean(squeeze(mean_neg_column_qsm(:,3,2,:)),1,'omitnan')';
rh_col4 = mean(squeeze(mean_neg_column_qsm(:,4,2,:)),1,'omitnan')';
rh_col5 = mean(squeeze(mean_neg_column_qsm(:,5,2,:)),1,'omitnan')';
rh_col6 = mean(squeeze(mean_neg_column_qsm(:,6,2,:)),1,'omitnan')';
rh_col7 = mean(squeeze(mean_neg_column_qsm(:,7,2,:)),1,'omitnan')';
rh_col8 = mean(squeeze(mean_neg_column_qsm(:,8,2,:)),1,'omitnan')';
rh_col9 = mean(squeeze(mean_neg_column_qsm(:,9,2,:)),1,'omitnan')';
rh_col10 = mean(squeeze(mean_neg_column_qsm(:,10,2,:)),1,'omitnan')';
rh_col11 = mean(squeeze(mean_neg_column_qsm(:,11,2,:)),1,'omitnan')';
rh_col12 = mean(squeeze(mean_neg_column_qsm(:,12,2,:)),1,'omitnan')';
rh_col13 = mean(squeeze(mean_neg_column_qsm(:,13,2,:)),1,'omitnan')';
rh_col14 = mean(squeeze(mean_neg_column_qsm(:,14,2,:)),1,'omitnan')';
rh_col15 = mean(squeeze(mean_neg_column_qsm(:,15,2,:)),1,'omitnan')';
rh_col16 = mean(squeeze(mean_neg_column_qsm(:,16,2,:)),1,'omitnan')';
rh_col17 = mean(squeeze(mean_neg_column_qsm(:,17,2,:)),1,'omitnan')';
rh_col18 = mean(squeeze(mean_neg_column_qsm(:,18,2,:)),1,'omitnan')';
rh_col19 = mean(squeeze(mean_neg_column_qsm(:,19,2,:)),1,'omitnan')';
rh_col20 = mean(squeeze(mean_neg_column_qsm(:,20,2,:)),1,'omitnan')';
rh_col21 = mean(squeeze(mean_neg_column_qsm(:,21,2,:)),1,'omitnan')';

col_rh_neg = table(rh_col1, rh_col2, rh_col3, rh_col4, rh_col5, rh_col6, rh_col7, rh_col8, rh_col9, rh_col10, rh_col11, rh_col12, rh_col13, rh_col14, rh_col15, rh_col6, rh_col7, rh_col18, rh_col19, rh_col20, rh_col21);
clear rh_*


%% QSM across depths, but by curvature (Figure 3) %%
% Positive QSM
for k1 = 1:34 % ROIs
    for k2 = 1:44 % Participants
        for k3 = 1:3 % Curvature
            for k4 = 1:2 % Hemisphere
                m = mean(QSM_ROIs_pos(k1,:,k3,k4,k2),2); % Takes the mean across the second dimension (depth)
                mean_pos_qsm(k1,k3,k4,k2) = m; 
            end
        end
    end
end

clear k1 k2 k3 k4 m

pos_crown_RH = mean(squeeze(mean_pos_qsm(:,1,2,:)),1,'omitnan')'; % Takes the mean across first dimension (ROI)
pos_crown_LH = mean(squeeze(mean_pos_qsm(:,1,1,:)),1,'omitnan')';

pos_bank_RH = mean(squeeze(mean_pos_qsm(:,2,2,:)),1,'omitnan')';
pos_bank_LH = mean(squeeze(mean_pos_qsm(:,2,1,:)),1,'omitnan')';

pos_fundus_RH = mean(squeeze(mean_pos_qsm(:,3,2,:)),1,'omitnan')';
pos_fundus_LH = mean(squeeze(mean_pos_qsm(:,3,1,:)),1,'omitnan')';

% Negative QSM
for k1 = 1:34
    for k2 = 1:44
        for k3 = 1:3
            for k4 = 1:2
                m = mean(QSM_ROIs_neg(k1,:,k3,k4,k2),2); % Takes the mean across the second dimension (depth)
                mean_neg_qsm(k1,k3,k4,k2) = m; 
            end
        end
    end
end

clear k1 k2 k3 k4 n p

neg_crown_RH = mean(squeeze(mean_neg_qsm(:,1,2,:)),1,'omitnan')'; % Takes the mean across first dimension (ROI)
neg_crown_LH = mean(squeeze(mean_neg_qsm(:,1,1,:)),1,'omitnan')';

neg_bank_RH = mean(squeeze(mean_neg_qsm(:,2,2,:)),1,'omitnan')';
neg_bank_LH = mean(squeeze(mean_neg_qsm(:,2,1,:)),1,'omitnan')';

neg_fundus_RH = mean(squeeze(mean_neg_qsm(:,3,2,:)),1,'omitnan')';
neg_fundus_LH = mean(squeeze(mean_neg_qsm(:,3,1,:)),1,'omitnan')';


%% Susceptibility by ROI %% 
% LH positive
ROI_1_f = squeeze(mean_pos_qsm(1,3,1,:));
ROI_5_b = squeeze(mean_pos_qsm(5,2,1,:));
ROI_5_f = squeeze(mean_pos_qsm(5,3,1,:));
ROI_7_f = squeeze(mean_pos_qsm(7,3,1,:));
ROI_8_f = squeeze(mean_pos_qsm(8,3,1,:));
ROI_12_f = squeeze(mean_pos_qsm(12,3,1,:));
ROI_13_f = squeeze(mean_pos_qsm(13,3,1,:));
ROI_16_f = squeeze(mean_pos_qsm(16,3,1,:));
ROI_16_b = squeeze(mean_pos_qsm(16,2,1,:));
ROI_21_c = squeeze(mean_pos_qsm(21,1,1,:));
ROI_21_b = squeeze(mean_pos_qsm(21,2,1,:));
ROI_21_f= squeeze(mean_pos_qsm(21,3,1,:));
ROI_25_c = squeeze(mean_pos_qsm(25,1,1,:));
ROI_27_c = squeeze(mean_pos_qsm(27,1,1,:));
ROI_27_b = squeeze(mean_pos_qsm(27,2,1,:));
ROI_29_c = squeeze(mean_pos_qsm(29,1,1,:));
ROI_33_c = squeeze(mean_pos_qsm(33,1,1,:));
ROI_33_b = squeeze(mean_pos_qsm(33,2,1,:));
ROI_33_f = squeeze(mean_pos_qsm(33,3,1,:));

lh_pos_ROI = table(ROI_1_f, ROI_5_b, ROI_5_f, ROI_7_f, ROI_8_f, ROI_12_f, ROI_13_f, ROI_16_f, ROI_16_b, ROI_21_c, ROI_21_b, ROI_21_f, ROI_25_c, ROI_27_c, ROI_27_b, ROI_29_c, ROI_33_c, ROI_33_b, ROI_33_f);
writetable(lh_pos_ROI);
clear ROI_*

% RH positive
ROI_4_b = squeeze(mean_pos_qsm(4,2,2,:));
ROI_12_f = squeeze(mean_pos_qsm(12,3,2,:));
ROI_13_b = squeeze(mean_pos_qsm(13,2,2,:));
ROI_25_c = squeeze(mean_pos_qsm(25,1,2,:));
ROI_30_b = squeeze(mean_pos_qsm(30,2,2,:));
ROI_30_f = squeeze(mean_pos_qsm(30,3,2,:));
ROI_33_b = squeeze(mean_pos_qsm(33,2,2,:));
ROI_33_c = squeeze(mean_pos_qsm(33,1,1,:));

rh_pos_ROI = table(ROI_4_b,ROI_12_f,ROI_13_b,ROI_25_c,ROI_30_b,ROI_30_f,ROI_33_b,ROI_33_c);
writetable(rh_pos_ROI);
clear ROI_*

% LH negative
ROI_10_b = squeeze(mean_neg_qsm(10,2,1,:));
ROI_11_f = squeeze(mean_neg_qsm(11,3,1,:));
ROI_12_b = squeeze(mean_neg_qsm(12,2,1,:));
ROI_12_f = squeeze(mean_neg_qsm(12,3,1,:));
ROI_13_b = squeeze(mean_neg_qsm(13,2,1,:));
ROI_13_f = squeeze(mean_neg_qsm(13,3,1,:));
ROI_15_b = squeeze(mean_neg_qsm(15,2,1,:));
ROI_16_f = squeeze(mean_neg_qsm(16,3,1,:));
ROI_18_f = squeeze(mean_neg_qsm(18,3,1,:));
ROI_23_c = squeeze(mean_neg_qsm(23,1,1,:));
ROI_25_c = squeeze(mean_neg_qsm(25,1,1,:));
ROI_25_b = squeeze(mean_neg_qsm(25,2,1,:));
ROI_25_f = squeeze(mean_neg_qsm(25,3,1,:));
ROI_27_f = squeeze(mean_neg_qsm(27,3,1,:));
ROI_30_b = squeeze(mean_neg_qsm(30,2,1,:));

lh_neg_ROI = table(ROI_10_b,ROI_11_f,ROI_12_b,ROI_12_f,ROI_13_b,ROI_13_f,ROI_15_b,ROI_16_f,ROI_18_f,ROI_23_c,ROI_25_c,ROI_25_b,ROI_25_f,ROI_27_f,ROI_30_b);
writetable(lh_neg_ROI);
clear ROI_*

% RH negative
ROI_12_f = squeeze(mean_neg_qsm(12,3,2,:));
ROI_13_f = squeeze(mean_neg_qsm(13,3,2,:));
ROI_16_b = squeeze(mean_neg_qsm(16,2,2,:));
ROI_24_b = squeeze(mean_neg_qsm(24,2,2,:));
ROI_25_b = squeeze(mean_neg_qsm(25,2,2,:));
ROI_25_f = squeeze(mean_neg_qsm(25,3,2,:));
ROI_31_c = squeeze(mean_neg_qsm(31,1,2,:));
ROI_33_f = squeeze(mean_neg_qsm(33,3,2,:));
ROI_34_b = squeeze(mean_neg_qsm(34,2,1,:));

rh_neg_ROI = table(ROI_12_f,ROI_13_f,ROI_16_b,ROI_24_b,ROI_25_b,ROI_25_f,ROI_31_c,ROI_33_f,ROI_34_b);
writetable(rh_neg_ROI);
clear ROI_*

%% Susceptibility by ROI: Figures 4 and 5 %%
load('/Volumes/Data/Madden/QSM.02/Analysis/columns/CustomMap.mat');

% LH positive
fh1 = figure(3)
tiledlayout(1,5,'Padding', 'compact', 'TileSpacing', 'compact');
ax1=nexttile;
hold on
colormap(ax1,CustomMap);
imagesc(cat(2,fdr_lh_crown_pos,fdr_lh_bank_pos,fdr_lh_fundus_pos),[0 0.05]);
set(gca,'xtick',[])
set(gca,'ytick',[])
yline(11.5)
yline(16.5)
yline(25.5)
yline(29.5)
xline(21.5)
xline(42.5)
xline(63.5)
ylim([0.5 size(fdr_lh_crown_pos,1)+0.5])
xlim([0.5 size(fdr_lh_crown_pos,2).*3+0.5])
yticks([1:size(fdr_lh_crown_pos,1)])
ax1.XAxis.FontSize = 20;
xlabel('Crown  Bank  Fundus')
ylabel('Cortical Gray Matter ROIs')
set(gca,'xaxisLocation','top')
set(gca,'TickLength',[0 0])
set(gca, 'YDir','reverse')
colorbar('eastoutside','Color','k')
set(gcf, 'InvertHardCopy', 'off'); 
set(gcf,'Color',[1 1 1]); 
set(gca, 'Color','k', 'XColor','k', 'YColor','k')
hold off
fh1.WindowState = 'maximized';

% RH positive
fh1 = figure(3)
tiledlayout(1,5,'Padding', 'compact', 'TileSpacing', 'compact'); 
ax1=nexttile;
hold on
colormap(ax1,CustomMap);
imagesc(cat(2,fdr_rh_crown_pos,fdr_rh_bank_pos,fdr_rh_fundus_pos),[0 0.05]);
set(gca,'ytick',[])
yline([11.5, 16.5, 25.5,29.5])
xline([21.5, 42.5, 63.5])
ylim([0.5 size(fdr_rh_crown_pos,1)+0.5])
xlim([0.5 size(fdr_rh_crown_pos,2).*3+0.5])
yticks([1:size(fdr_rh_crown_pos,1)])
ax1.XAxis.FontSize = 12;
ax1.YAxis.FontSize = 12;
xlabel('Crown  Bank  Fundus')
ylabel('Cortical Gray Matter ROIs')
set(gca, 'YDir','reverse')
colorbar('eastoutside','Color','k','FontSize',12)
set(gcf, 'InvertHardCopy', 'off'); 
set(gcf,'Color',[1 1 1]); 
set(gca, 'Color','k', 'XColor','k', 'YColor','k')
xticks([5.5 10.5 15.5 26.5 31.5 36.5 47.5 52.5 57.5])
xticklabels({'25%','50%','75%','25%','50%','75%','25%','50%','75%'})
xtickangle(90)
hold off
fh1.WindowState = 'maximized';

% LH negative
fh1 = figure(3)
tiledlayout(1,5,'Padding', 'compact', 'TileSpacing', 'compact');
ax1=nexttile;
hold on
colormap(ax1,CustomMap);
imagesc(cat(2,fdr_lh_crown_neg,fdr_lh_bank_neg,fdr_lh_fundus_neg),[0 0.05]);
set(gca,'xtick',[])
set(gca,'ytick',[])
yline(11.5)
yline(16.5)
yline(25.5)
yline(29.5)
xline(21.5)
xline(42.5)
xline(63.5)
ylim([0.5 size(fdr_lh_crown_neg,1)+0.5])
xlim([0.5 size(fdr_lh_crown_neg,2).*3+0.5])
yticks([1:size(fdr_lh_crown_neg,1)])
ax1.XAxis.FontSize = 20;
xlabel('Crown  Bank  Fundus')
ylabel('Cortical Gray Matter ROIs')
set(gca,'xaxisLocation','top')
set(gca,'TickLength',[0 0])
set(gca, 'YDir','reverse')
colorbar('eastoutside','Color','k')
set(gcf, 'InvertHardCopy', 'off'); 
set(gcf,'Color',[1 1 1]); 
set(gca, 'Color','k', 'XColor','k', 'YColor','k')
hold off
fh1.WindowState = 'maximized';

% RH negative
fh1 = figure(3)
tiledlayout(1,5,'Padding', 'compact', 'TileSpacing', 'compact');
ax1=nexttile;
hold on
colormap(ax1,CustomMap);
imagesc(cat(2,fdr_rh_crown_neg,fdr_rh_bank_neg,fdr_rh_fundus_neg),[0 0.05]);
set(gca,'xtick',[])
set(gca,'ytick',[])
yline(11.5)
yline(16.5)
yline(25.5)
yline(29.5)
xline(21.5)
xline(42.5)
xline(63.5)
ylim([0.5 size(fdr_rh_crown_neg,1)+0.5])
xlim([0.5 size(fdr_rh_crown_neg,2).*3+0.5])
yticks([1:size(fdr_rh_crown_neg,1)])
ax1.XAxis.FontSize = 20;
xlabel('Crown  Bank  Fundus')
ylabel('Cortical Gray Matter ROIs')
set(gca,'xaxisLocation','top')
set(gca,'TickLength',[0 0])
set(gca, 'YDir','reverse')
colorbar('eastoutside','Color','k')
set(gcf, 'InvertHardCopy', 'off'); 
set(gcf,'Color',[1 1 1]); 
set(gca, 'Color','k', 'XColor','k', 'YColor','k')
hold off
fh1.WindowState = 'maximized';

