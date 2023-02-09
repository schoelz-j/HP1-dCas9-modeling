#!/usr/bin/env python
# coding: utf-8

# # 'HP1-Only' Expression Models

import pandas as pd
import numpy as np
import sklearn
import warnings

from sklearn.model_selection import train_test_split
from sklearn.linear_model import ElasticNet
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import cross_val_predict

def main():
    # Model inputs
    # TSS HP1a, HP1B and HP1C from S2-DRSC cells measured by ChIP-seq
    # (quantified by deepTools). The target variable in this dataset is gene
    # expression in S2-DRSC cells measured by RNA-seq (quantified by Salmon).
    # Input datasets were prepared using HP1-only_Models_DataPrep.R
    # Data is not scaled

    s2 = pd.read_csv('data/model_inputs/preprocessed_data_for_modeling/HP1-only-S2-Inputs.csv')
    bg3 = pd.read_csv('data/model_inputs/preprocessed_data_for_modeling/HP1-only-BG3-Inputs.csv')
    cme = pd.read_csv('data/model_inputs/preprocessed_data_for_modeling/HP1-only-CME-Inputs.csv')

    # First split the data into training and test sets. For both sets, work with
    #  a copy where the FBgn ID column has been dropped. Leaving an original that
    # still has the FBgn column will be useful when writing model output.
    # Train Test Split

    s2_train, s2_test = train_test_split(s2, test_size = 0.2,
                                         random_state = 42)

    y_train_s2 = s2_train['log_TPM']
    y_test_s2 = s2_test['log_TPM']

    for set_ in (s2_train, s2_test):
        set_.drop("log_TPM", axis=1, inplace=True)

    # a version of each of these that retains the FBgn identifier
    s2_train2 = s2_train.drop(['FBgn'], axis = 1)
    s2_test2 = s2_test.drop(['FBgn'], axis = 1)

    # Data transformations and grid search
    # The pipeline below will scale the training data and fine-tune model
    # hyperparameters to values that generalize best.
    # Grid search to optimize model parameters:
    # alpha - penalty strength
    # l1_ratio - weights for lasso vs ridge regularization mean

    parameters = {'model__alpha': [0.1, 0.3, 0.5, 0.7, 0.9],
                  'model__l1_ratio': [0.1, 0.25, 0.5, 0.75, 0.9]}

    steps = list()
    steps.append(('scaler', StandardScaler()))
    steps.append(('model', ElasticNet()))
    pipeline = Pipeline(steps = steps)
    s2_grid_search = GridSearchCV(pipeline, parameters, n_jobs=-1, verbose=1)
    s2_grid_search.fit(s2_train2, y_train_s2)
    best_parameters = s2_grid_search.best_estimator_.get_params()

    # Generate Predictions
    # use fine-tuned model to generate predictions.
    s2_y_train_pred = cross_val_predict(s2_grid_search.best_estimator_, s2_train2, y_train_s2,
                                        cv=10)

    # Test set predictions
    s2_y_test_pred = s2_grid_search.best_estimator_.predict(s2_test2)

    # Extract coefficients
    feature_names = ['HP1a', 'HP1B', 'HP1C', 'A_B', 'A_C', 'B_C', 'A_B_C']

    s2_train_out = pd.DataFrame({'FBgn': s2_train['FBgn'],
                                 'Y_true': y_train_s2,
                                 'Y_pred': s2_y_train_pred})
    s2_train_out.to_csv('data/model_outputs/hp1_only/predictions/HP1-only-S2-training-output.csv', index = False)
    s2_test_out = pd.DataFrame({'FBgn': s2_test['FBgn'],
                                'Y_true': y_test_s2,
                                'Y_pred': s2_y_test_pred})
    s2_test_out.to_csv('data/model_outputs/hp1_only/predictions/HP1-only-S2-testing-output.csv', index = False)
    s2_coefs = pd.DataFrame({'Feature': feature_names,
                             'Coefficient': s2_grid_search.best_estimator_.named_steps['model'].coef_})
    s2_coefs.to_csv('data/model_outputs/hp1_only/predictions/HP1-only-S2-coefficients.csv', index = False)

    # Permutations
    s2_permut = pd.read_csv('data/model_inputs/preprocessed_data_for_modeling/HP1-only-S2-permut-Inputs.csv'))
    s2_ptrain, s2_ptest = train_test_split(s2_permut, test_size = 0.2,
                                           random_state = 42)
    yp_train_s2 = s2_ptrain['log_TPM']
    yp_test_s2 = s2_ptest['log_TPM']

    for set_ in (s2_ptrain, s2_ptest):
        set_.drop("log_TPM", axis=1, inplace=True)

    s2p_grid_search = GridSearchCV(pipeline, parameters, n_jobs=1, verbose=1)
    s2p_grid_search.fit(s2_ptrain, yp_train_s2)
    best_parameters = s2p_grid_search.best_estimator_.get_params()

    # Training set predictions
    s2p_y_train_pred = cross_val_predict(s2p_grid_search.best_estimator_, s2_ptrain, yp_train_s2,
                                         cv=10)

    # Test set predictions
    s2p_y_test_pred = s2p_grid_search.best_estimator_.predict(s2_ptest)

    # Write Output
    s2p_train_out = pd.DataFrame({'Y_true': yp_train_s2,
                                  'Y_pred': s2p_y_train_pred})
    s2p_train_out.to_csv('data/model_outputs/hp1_only/predictions/HP1-only-S2-permut-training-output.csv', index = False)
    s2p_test_out = pd.DataFrame({'Y_true': yp_test_s2,
                                 'Y_pred': s2p_y_test_pred})
    s2p_test_out.to_csv('data/model_outputs/hp1_only/predictions/HP1-only-S2-permut-testing-output.csv', index = False)
    s2p_coefs = pd.DataFrame({'Feature': feature_names,
                              'Coefficient': s2p_grid_search.best_estimator_.named_steps['model'].coef_})
    s2p_coefs.to_csv('data/model_outputs/hp1_only/predictions/HP1-only-S2-permut-coefficients.csv', index = False)

    # ### BG3 and CME-W1 cell data
    #
    # Re-run this analysis with the original datasets before performing permutations.
    bg3_train, bg3_test = train_test_split(bg3, test_size = 0.2,
                                           random_state = 42)


    y_train_bg3 = bg3_train['log_TPM']
    y_test_bg3 = bg3_test['log_TPM']

    for set_ in (bg3_train, bg3_test):
        set_.drop("log_TPM", axis=1, inplace=True)

        # Just keeping a version of each of these that retains the FBgn identifier
    bg3_train2 = bg3_train.drop(['FBgn'], axis = 1)
    bg3_test2 = bg3_test.drop(['FBgn'], axis = 1)
    bg3_grid_search = GridSearchCV(pipeline, parameters, n_jobs=1, verbose=1)
    bg3_grid_search.fit(bg3_train2, y_train_bg3)
    best_parameters = bg3_grid_search.best_estimator_.get_params()

    bg3_y_train_pred = cross_val_predict(bg3_grid_search.best_estimator_, bg3_train2, y_train_bg3,
                                   cv=10)
    bg3_y_test_pred = bg3_grid_search.best_estimator_.predict(bg3_test2)

    bg3_train_out = pd.DataFrame({'FBgn': bg3_train['FBgn'],
                                  'Y_true': y_train_bg3,
                                  'Y_pred': bg3_y_train_pred})
    bg3_train_out.to_csv('data/model_outputs/hp1_only/predictions/HP1-only-BG3-training-output.csv', index = False)
    bg3_test_out = pd.DataFrame({'FBgn': bg3_test['FBgn'],
                                 'Y_true': y_test_bg3,
                                 'Y_pred': bg3_y_test_pred})
    bg3_test_out.to_csv('data/model_outputs/hp1_only/predictions/HP1-only-BG3-testing-output.csv', index = False)
    bg3_coefs = pd.DataFrame({'Feature': feature_names,
                              'Coefficient': bg3_grid_search.best_estimator_.named_steps['model'].coef_})
    bg3_coefs.to_csv('data/model_outputs/hp1_only/predictions/HP1-only-BG3-coefficients.csv', index = False)
    cme_train, cme_test = train_test_split(cme, test_size = 0.2,
                                           random_state = 42)
    y_train_cme = cme_train['log_TPM']
    y_test_cme = cme_test['log_TPM']

    for set_ in (cme_train, cme_test):
        set_.drop("log_TPM", axis=1, inplace=True)

    # Just keeping a version of each of these that retains the FBgn identifier
    cme_train2 = cme_train.drop(['FBgn'], axis = 1)
    cme_test2 = cme_test.drop(['FBgn'], axis = 1)
    cme_grid_search = GridSearchCV(pipeline, parameters, n_jobs=1, verbose=1)
    cme_grid_search.fit(cme_train2, y_train_cme)
    best_parameters = cme_grid_search.best_estimator_.get_params()

    cme_y_train_pred = cross_val_predict(cme_grid_search.best_estimator_, cme_train2, y_train_cme,
                                         cv=10)
    cme_y_test_pred = cme_grid_search.best_estimator_.predict(cme_test2)

    cme_train_out = pd.DataFrame({'FBgn': cme_train['FBgn'],
                                  'Y_true': y_train_cme,
                                  'Y_pred': cme_y_train_pred})
    cme_train_out.to_csv('data/model_outputs/hp1_only/predictions/HP1-only-CME-training-output.csv', index = False)
    cme_test_out = pd.DataFrame({'FBgn': cme_test['FBgn'],
                                 'Y_true': y_test_cme,
                                 'Y_pred': cme_y_test_pred})
    cme_test_out.to_csv('data/model_outputs/hp1_only/predictions/HP1-only-CME-testing-output.csv', index = False)
    cme_coefs = pd.DataFrame({'Feature': feature_names,
                              'Coefficient': cme_grid_search.best_estimator_.named_steps['model'].coef_})
    cme_coefs.to_csv('data/model_outputs/hp1_only/predictions/HP1-only-CME-coefficients.csv', index = False)

    bg3_permut = pd.read_csv('../datasets/model_inputs/permutation_inputs/HP1-only-BG3-permut-Inputs.csv')

    bg3_ptrain, bg3_ptest = train_test_split(bg3_permut, test_size = 0.2,
                                             random_state = 42)
    yp_train_bg3 = bg3_ptrain['log_TPM']
    yp_test_bg3 = bg3_ptest['log_TPM']

    for set_ in (bg3_ptrain, bg3_ptest):
        set_.drop("log_TPM", axis=1, inplace=True)
    bg3p_grid_search = GridSearchCV(pipeline, parameters, n_jobs=1, verbose=1)
    bg3p_grid_search.fit(bg3_ptrain, yp_train_bg3)
    best_parameters = bg3p_grid_search.best_estimator_.get_params()
    bg3p_y_train_pred = cross_val_predict(bg3p_grid_search.best_estimator_, bg3_ptrain, yp_train_bg3,
                                          cv=10)
    bg3p_y_test_pred = bg3p_grid_search.best_estimator_.predict(bg3_ptest)
    bg3p_train_out = pd.DataFrame({'Y_true': yp_train_bg3,
                                   'Y_pred': bg3p_y_train_pred})
    bg3p_train_out.to_csv('data/model_outputs/hp1_only/predictions/HP1-only-BG3-permut-training-output.csv', index = False)
    bg3p_test_out = pd.DataFrame({'Y_true': yp_test_bg3,
                                  'Y_pred': bg3p_y_test_pred})
    bg3p_test_out.to_csv('data/model_outputs/hp1_only/predictions/HP1-only-BG3-permut-testing-output.csv', index = False)
    bg3p_coefs = pd.DataFrame({'Feature': feature_names,
                               'Coefficient': bg3p_grid_search.best_estimator_.named_steps['model'].coef_})
    bg3p_coefs.to_csv('data/model_outputs/hp1_only/predictions/HP1-only-BG3-permut-coefficients.csv', index = False)

    cme_permut = pd.read_csv('../datasets/model_inputs/permutation_inputs/HP1-only-CME-permut-Inputs.csv')

    cme_ptrain, cme_ptest = train_test_split(cme_permut, test_size = 0.2,
                                             random_state = 42)
    yp_train_cme = cme_ptrain['log_TPM']
    yp_test_cme = cme_ptest['log_TPM']

    for set_ in (cme_ptrain, cme_ptest):
        set_.drop("log_TPM", axis=1, inplace=True)
    cmep_grid_search = GridSearchCV(pipeline, parameters, n_jobs=1, verbose=1)
    cmep_grid_search.fit(cme_ptrain, yp_train_cme)
    best_parameters = cmep_grid_search.best_estimator_.get_params()
    cmep_y_train_pred = cross_val_predict(cmep_grid_search.best_estimator_, cme_ptrain, yp_train_cme,
                                          cv=10)
    cmep_y_test_pred = cmep_grid_search.best_estimator_.predict(cme_ptest)

    cmep_train_out = pd.DataFrame({'Y_true': yp_train_cme,
                                   'Y_pred': cmep_y_train_pred})
    cmep_train_out.to_csv('data/model_outputs/hp1_only/predictions/HP1-only-CME-permut-training-output.csv', index = False)
    cmep_test_out = pd.DataFrame({'Y_true': yp_test_cme,
                                  'Y_pred': cmep_y_test_pred})
    cmep_test_out.to_csv('data/model_outputs/hp1_only/predictions/HP1-only-CME-permut-testing-output.csv', index = False)
    cmep_coefs = pd.DataFrame({'Feature': feature_names,
                               'Coefficient': cmep_grid_search.best_estimator_.named_steps['model'].coef_})
    cmep_coefs.to_csv('data/model_outputs/hp1_only/predictions/HP1-only-CME-permut-coefficients.csv', index = False)

     # ## dCas9 Simulations
    a_cpr11b = pd.read_csv('../datasets/dcas9_simulations/cpr11b_hp1a_hp1_only_inputs.csv')
    b_cpr11b = pd.read_csv('../datasets/dcas9_simulations/cpr11b_hp1b_hp1_only_inputs.csv')
    c_cpr11b = pd.read_csv('../datasets/dcas9_simulations/cpr11b_hp1c_hp1_only_inputs.csv')

    a_dgt3 = pd.read_csv('../datasets/dcas9_simulations/dgt3_hp1a_hp1_only_inputs.csv')
    b_dgt3 = pd.read_csv('../datasets/dcas9_simulations/dgt3_hp1b_hp1_only_inputs.csv')
    c_dgt3 = pd.read_csv('../datasets/dcas9_simulations/dgt3_hp1c_hp1_only_inputs.csv')

    a_ceca1 = pd.read_csv('../datasets/dcas9_simulations/ceca1_hp1a_hp1_only_inputs.csv')
    b_ceca1 = pd.read_csv('../datasets/dcas9_simulations/ceca1_hp1b_hp1_only_inputs.csv')
    c_ceca1 = pd.read_csv('../datasets/dcas9_simulations/ceca1_hp1c_hp1_only_inputs.csv')

    a_mtk = pd.read_csv('../datasets/dcas9_simulations/mtk_hp1a_hp1_only_inputs.csv')
    b_mtk = pd.read_csv('../datasets/dcas9_simulations/mtk_hp1b_hp1_only_inputs.csv')
    c_mtk = pd.read_csv('../datasets/dcas9_simulations/mtk_hp1c_hp1_only_inputs.csv')

    a_egr = pd.read_csv('../datasets/dcas9_simulations/egr_hp1a_hp1_only_inputs.csv')
    a_pyr = pd.read_csv('../datasets/dcas9_simulations/pyr_hp1a_hp1_only_inputs.csv')
    a_mats = pd.read_csv('../datasets/dcas9_simulations/mats_hp1a_hp1_only_inputs.csv')

    b_rab3 = pd.read_csv('../datasets/dcas9_simulations/rab3_hp1b_hp1_only_inputs.csv')
    b_shaw = pd.read_csv('../datasets/dcas9_simulations/shaw_hp1b_hp1_only_inputs.csv')
    b_cg76 = pd.read_csv('../datasets/dcas9_simulations/cg76_hp1b_hp1_only_inputs.csv')

    c_alpha = pd.read_csv('../datasets/dcas9_simulations/alpha_hp1c_hp1_only_inputs.csv')
    c_cg26 = pd.read_csv('../datasets/dcas9_simulations/cg26_hp1c_hp1_only_inputs.csv')

    a_cpr11b_preds = s2_grid_search.best_estimator_.predict(a_cpr11b)
    b_cpr11b_preds = s2_grid_search.best_estimator_.predict(b_cpr11b)
    c_cpr11b_preds = s2_grid_search.best_estimator_.predict(c_cpr11b)

    a_dgt3_preds = s2_grid_search.best_estimator_.predict(a_dgt3)
    b_dgt3_preds = s2_grid_search.best_estimator_.predict(b_dgt3)
    c_dgt3_preds = s2_grid_search.best_estimator_.predict(c_dgt3)

    a_ceca1_preds = s2_grid_search.best_estimator_.predict(a_ceca1)
    b_ceca1_preds = s2_grid_search.best_estimator_.predict(b_ceca1)
    c_ceca1_preds = s2_grid_search.best_estimator_.predict(c_ceca1)

    a_mtk_preds = s2_grid_search.best_estimator_.predict(a_mtk)
    b_mtk_preds = s2_grid_search.best_estimator_.predict(b_mtk)
    c_mtk_preds = s2_grid_search.best_estimator_.predict(c_mtk)

    a_egr_preds = s2_grid_search.best_estimator_.predict(a_egr)
    a_pyr_preds = s2_grid_search.best_estimator_.predict(a_pyr)
    a_mats_preds = s2_grid_search.best_estimator_.predict(a_mats)

    b_rab3_preds = s2_grid_search.best_estimator_.predict(b_rab3)
    b_shaw_preds = s2_grid_search.best_estimator_.predict(b_shaw)
    b_cg76_preds = s2_grid_search.best_estimator_.predict(b_cg76)

    c_alpha_preds = s2_grid_search.best_estimator_.predict(c_alpha)
    c_cg26_preds = s2_grid_search.best_estimator_.predict(c_cg26)

    cpr11b_out = pd.DataFrame({'HP1a': a_cpr11b_preds,
                               'HP1B': b_cpr11b_preds,
                               'HP1C': c_cpr11b_preds})
    dgt3_out = pd.DataFrame({'HP1a': a_dgt3_preds,
                             'HP1B': b_dgt3_preds,
                             'HP1C': c_dgt3_preds})
    ceca1_out = pd.DataFrame({'HP1a': a_ceca1_preds,
                              'HP1B': b_ceca1_preds,
                              'HP1C': c_ceca1_preds})
    mtk_out = pd.DataFrame({'HP1a': a_mtk_preds,
                            'HP1B': b_mtk_preds,
                            'HP1C': c_mtk_preds})
    hp1a_out = pd.DataFrame({'egr': a_egr_preds,
                             'pyr': a_pyr_preds,
                             'mats': a_mats_preds})
    hp1b_out = pd.DataFrame({'rab3': b_rab3_preds,
                             'shaw': b_shaw_preds,
                             'cg76': b_cg76_preds})
    hp1c_out = pd.DataFrame({'alpha': c_alpha_preds,
                             'cg26': c_cg26_preds})
    cpr11b_out.to_csv('data/model_outputs/hp1_only/predictionscpr11b_hp1_only_simulated.csv')
    dgt3_out.to_csv('data/model_outputs/hp1_only/predictionsdgt3_hp1_only_simulated.csv')
    ceca1_out.to_csv('data/model_outputs/hp1_only/predictionsceca1_hp1_only_simulated.csv')
    mtk_out.to_csv('data/model_outputs/hp1_only/predictionsmtk_hp1_only_simulated.csv')
    hp1a_out.to_csv('data/model_outputs/hp1_only/predictions/HP1a_hp1_only_simulated.csv')
    hp1b_out.to_csv('data/model_outputs/hp1_only/predictions/HP1b_hp1_only_simulated.csv')
    hp1c_out.to_csv('data/model_outputs/hp1_only/predictions/HP1c_hp1_only_simulated.csv')

    # Test set predictions
    # Regular model:
    elastic_net.fit(train_set2, y_train)
    ytest_pred = elastic_net.predict(test_set2)
    test_resid = ytest_pred - y_test
    hp1_test = {'Y_pred': ytest_pred,
                'Y_true': y_test,
                'HP1a': test_set2['HP1a'],
                'HP1B': test_set2['HP1B'],
                'HP1C': test_set2['HP1C'],
                'FBgn': test_set['FBgn']}
    hp1_test_df = pd.DataFrame(hp1_test)
    hp1_test_df.to_csv('HP1Only_test_predictions.csv', index=False)

if __name__ == '__main__':
    main()
