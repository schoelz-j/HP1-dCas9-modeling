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

# In[12]:


bg3_train, bg3_test = train_test_split(bg3, test_size = 0.2,
                                       random_state = 42)


y_train_bg3 = bg3_train['log_TPM']
y_test_bg3 = bg3_test['log_TPM']

for set_ in (bg3_train, bg3_test):
    set_.drop("log_TPM", axis=1, inplace=True)

    # Just keeping a version of each of these that retains the FBgn identifier
bg3_train2 = bg3_train.drop(['FBgn'], axis = 1)
bg3_test2 = bg3_test.drop(['FBgn'], axis = 1)

if __name__ == "__main__":
    bg3_grid_search = GridSearchCV(pipeline, parameters, n_jobs=1, verbose=1)

    #print("Performing grid search...")
    #print("pipeline:", [name for name, _ in pipeline.steps])
    #print("parameters:")
    pprint(parameters)

    t0 = time()
    bg3_grid_search.fit(bg3_train2, y_train_bg3)
    #print("done in %0.3fs" % (time() - t0))
    #print()

    #print("Best score: %0.3f" % bg3_grid_search.best_score_)
    #print("Best parameters set:")
    best_parameters = bg3_grid_search.best_estimator_.get_params()

    #for param_name in sorted(parameters.keys()):
    #    print("\t%s: %r" % (param_name, best_parameters[param_name]))

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

# In[13]:


cme_train, cme_test = train_test_split(cme, test_size = 0.2,
                                       random_state = 42)


y_train_cme = cme_train['log_TPM']
y_test_cme = cme_test['log_TPM']

for set_ in (cme_train, cme_test):
    set_.drop("log_TPM", axis=1, inplace=True)

    # Just keeping a version of each of these that retains the FBgn identifier
cme_train2 = cme_train.drop(['FBgn'], axis = 1)
cme_test2 = cme_test.drop(['FBgn'], axis = 1)

if __name__ == "__main__":
    cme_grid_search = GridSearchCV(pipeline, parameters, n_jobs=1, verbose=1)

    #print("Performing grid search...")
    #print("pipeline:", [name for name, _ in pipeline.steps])
    #print("parameters:")
    pprint(parameters)

    t0 = time()
    cme_grid_search.fit(cme_train2, y_train_cme)
    #print("done in %0.3fs" % (time() - t0))
    #print()

    #print("Best score: %0.3f" % cme_grid_search.best_score_)
    #print("Best parameters set:")
    best_parameters = cme_grid_search.best_estimator_.get_params()

    #for param_name in sorted(parameters.keys()):
    #    print("\t%s: %r" % (param_name, best_parameters[param_name]))

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

# ### BG3 and CME cells: dataset permutations
#
# Finally, run model fitting and predictions on permuted dataset to check the robustness of model interpretations.

# In[14]:


bg3_permut = pd.read_csv('../datasets/model_inputs/permutation_inputs/HP1-only-BG3-permut-Inputs.csv')

bg3_ptrain, bg3_ptest = train_test_split(bg3_permut, test_size = 0.2,
                                       random_state = 42)


yp_train_bg3 = bg3_ptrain['log_TPM']
yp_test_bg3 = bg3_ptest['log_TPM']

for set_ in (bg3_ptrain, bg3_ptest):
    set_.drop("log_TPM", axis=1, inplace=True)

if __name__ == "__main__":
    bg3p_grid_search = GridSearchCV(pipeline, parameters, n_jobs=1, verbose=1)

    #print("Performing grid search...")
    #print("pipeline:", [name for name, _ in pipeline.steps])
    #print("parameters:")
    pprint(parameters)

    t0 = time()
    bg3p_grid_search.fit(bg3_ptrain, yp_train_bg3)
    #print("done in %0.3fs" % (time() - t0))
    #print()

    #print("Best score: %0.3f" % bg3p_grid_search.best_score_)
    #print("Best parameters set:")
    best_parameters = bg3p_grid_search.best_estimator_.get_params()

    #for param_name in sorted(parameters.keys()):
    #    print("\t%s: %r" % (param_name, best_parameters[param_name]))

# Training set predictions
bg3p_y_train_pred = cross_val_predict(bg3p_grid_search.best_estimator_, bg3_ptrain, yp_train_bg3,
                                   cv=10)

# Test set predictions
bg3p_y_test_pred = bg3p_grid_search.best_estimator_.predict(bg3_ptest)

# Write Output
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

if __name__ == "__main__":
    cmep_grid_search = GridSearchCV(pipeline, parameters, n_jobs=1, verbose=1)

    #print("Performing grid search...")
    #print("pipeline:", [name for name, _ in pipeline.steps])
    #print("parameters:")
    #pprint(parameters)

    t0 = time()
    cmep_grid_search.fit(cme_ptrain, yp_train_cme)
    #print("done in %0.3fs" % (time() - t0))
    #print()

    #print("Best score: %0.3f" % cmep_grid_search.best_score_)
    #print("Best parameters set:")
    best_parameters = cmep_grid_search.best_estimator_.get_params()

    #for param_name in sorted(parameters.keys()):
    #    print("\t%s: %r" % (param_name, best_parameters[param_name]))

# Training set predictions
cmep_y_train_pred = cross_val_predict(cmep_grid_search.best_estimator_, cme_ptrain, yp_train_cme,
                                   cv=10)

# Test set predictions
cmep_y_test_pred = cmep_grid_search.best_estimator_.predict(cme_ptest)

# Write Output
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
#
# Generate predictions from simulated data

# In[15]:


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

# In[8]:


# Model with no HP1B data
no_b_set = train_set2.drop(['HP1B'], axis = 1)
no_b_set = no_b_set.drop(['A_B'], axis = 1)
no_b_set = no_b_set.drop(['A_B_C'], axis = 1)
no_b_set = no_b_set.drop(['B_C'], axis = 1)

no_b_y_pred = cross_val_predict(elastic_net, no_b_set, y_train, cv=10)
no_b_resid = no_b_y_pred - y_train

no_b_sse = sum(no_b_resid**2)
print(no_b_sse)

# In[9]:


no_c_set = train_set2.drop(['HP1C'], axis = 1)
no_c_set = no_c_set.drop(['A_C'], axis = 1)
no_c_set = no_c_set.drop(['A_B_C'], axis = 1)
no_c_set = no_c_set.drop(['B_C'], axis = 1)


no_c_y_pred = cross_val_predict(elastic_net, no_c_set, y_train, cv=10)
no_c_resid = no_c_y_pred - y_train

no_c_sse = sum(no_c_resid**2)
print(no_c_sse)

# In[10]:


plt.scatter(x = no_b_resid, y = no_c_resid, alpha = 0.3)
plt.show

# In[11]:


no_a_set = train_set2.drop(['HP1a'], axis = 1)
no_a_set = no_a_set.drop(['A_C'], axis = 1)
no_a_set = no_a_set.drop(['A_B_C'], axis = 1)
no_a_set = no_a_set.drop(['A_B'], axis = 1)

no_a_y_pred = cross_val_predict(elastic_net, no_a_set, y_train, cv=10)
no_a_resid = no_a_y_pred - y_train

no_a_sse = sum(no_a_resid**2)
print(no_a_sse)

# In[12]:


# Individual variables as predictors in linear models
from sklearn import linear_model
lin_model = linear_model.LinearRegression()
hp1b = train_set2.drop(['HP1C', 'HP1a', 'A_C', 'A_B',
                      'B_C', 'A_B_C'], axis = 1)
lin_model.fit(hp1b, y_train)
b_y_pred = lin_model.predict(hp1b)
b_resid = b_y_pred - y_train

hp1c = train_set2.drop(['HP1B', 'HP1a', 'A_C', 'A_B',
                      'B_C', 'A_B_C'], axis = 1)
lin_model.fit(hp1c, y_train)
c_y_pred = lin_model.predict(hp1c)
c_resid = c_y_pred - y_train

hp1a = train_set2.drop(['HP1C', 'HP1B', 'A_C', 'A_B',
                      'B_C', 'A_B_C'], axis = 1)
lin_model.fit(hp1a, y_train)
a_y_pred = lin_model.predict(hp1a)
a_resid = a_y_pred - y_train

print(sum(b_resid**2))
print(sum(c_resid**2))
print(sum(a_resid**2))

# In[13]:


plt.scatter(x = b_resid, y = c_resid, alpha = 0.3)
plt.show()

# In[14]:


plt.scatter(x = b_resid, y = a_resid, alpha = 0.3)
plt.show()

# In[15]:


plt.scatter(x = c_resid, y = a_resid, alpha = 0.3)
plt.show()

# In[16]:


# Let's do the same of HP1B vs HP1a
ba_cov_diff = train_set['HP1B'] - train_set['HP1a']
bc_cov_diff = train_set['HP1B'] - train_set['HP1C']
plt.scatter(x = ba_cov_diff, y = a_resid, alpha = 0.3)
plt.show()

# In[17]:


plt.scatter(x = ba_cov_diff, y = b_resid, alpha = 0.3)
plt.show()

# In[18]:


# So now I just want to add predictions, residuals, cov diffs
data = {'FBgn': train_set['FBgn'],
       'HP1a': train_set['HP1a'],
       'HP1B': train_set['HP1B'],
       'HP1C': train_set['HP1C'],
       'A_B': train_set['A_B'],
       'A_C': train_set['A_C'],
       'B_C': train_set['B_C'],
       'A_B_C': train_set['A_B_C'],
       'Y_true': y_train,
       'Y_Pred_Full': y_pred,
       'Y_Pred_HP1a': a_y_pred,
       'Y_Pred_HP1B': b_y_pred,
       'Y_Pred_HP1C': c_y_pred,
        'Full_Resid': resid,
       'HP1a_Resid': a_resid,
       'HP1B_Resid': b_resid,
       'HP1C_Resid': c_resid,
       'HP1B-HP1C': bc_cov_diff,
       'HP1B-HP1a': ba_cov_diff}
output = pd.DataFrame(data)
output.head()

# In[19]:


# Okay now I just need to output this thing
output.to_csv('Baseline_HP1_EN_Predictions.csv', index=False)

# In[20]:


# Train a model on permuted data to test robustness of model 1
random = pd.read_csv('HP1only_randomized.csv')
random.head()

# In[21]:


rhp1_train, rhp1_test = train_test_split(random, test_size = 0.2,
                                       random_state = 42)
rhp1_y_train = rhp1_train['log_S2']
rhp1_y_test = rhp1_test['log_S2']
for set_ in (rhp1_train, rhp1_test):
    set_.drop("log_S2", axis=1, inplace=True)

rhp1_train.head()

# In[22]:


random_enet = GridSearchCV(ElasticNet(), enet_grid, scoring = 'neg_mean_squared_error', cv=10,
                          refit=True)
random_enet.fit(rhp1_train, rhp1_y_train)
print(random_enet.best_params_)

# In[23]:


random_enet_final = ElasticNet(alpha = 0.9, l1_ratio = 0.0)
ry_pred = cross_val_predict(random_enet_final, rhp1_train, rhp1_y_train, cv=10)

rand_resid = ry_pred - rhp1_y_train
print(sum(rand_resid**2))

# In[24]:


random_output = {'Y_true': rhp1_y_train,
                'Y_pred': ry_pred,
                'Resid': rand_resid}
random_output_df = pd.DataFrame(random_output)
random_output_df.to_csv('HP1Only_Permutation_Predictions.csv', index=False)

# In[25]:


# Test set predictions
# Regular model:
elastic_net.fit(train_set2, y_train)
ytest_pred = elastic_net.predict(test_set2)
test_resid = ytest_pred - y_test
print(sum(test_resid**2))

# In[26]:


# Write out these predictions
hp1_test = {'Y_pred': ytest_pred,
           'Y_true': y_test,
           'HP1a': test_set2['HP1a'],
           'HP1B': test_set2['HP1B'],
           'HP1C': test_set2['HP1C'],
           'FBgn': test_set['FBgn']}
hp1_test_df = pd.DataFrame(hp1_test)
hp1_test_df.to_csv('HP1Only_test_predictions.csv', index=False)

# In[27]:


# Randomized model
random_enet_final.fit(rhp1_train, rhp1_y_train)
ry_test_pred = random_enet_final.predict(rhp1_test)
r_test_resid = ry_test_pred - rhp1_y_test
print(sum(r_test_resid**2))

# In[28]:


rhp1_test = {'Y_pred': ry_test_pred,
            'Y_true': rhp1_y_test,
            'Resid': r_test_resid}
rhp1_test_df = pd.DataFrame(rhp1_test)
rhp1_test_df.to_csv('HP1Only_Permuted_test_predictions.csv', index = False)

# In[29]:


# Obtain coefficients of each model
feature_names = ['HP1a', 'HP1B', 'HP1C', 'A_B', 'A_C', 'B_C', 'A_B_C']
print(elastic_net.coef_)
print(random_enet_final.coef_)

# In[30]:


coefs_data = {'Variable_Name': feature_names,
             'HP1_coef': elastic_net.coef_,
             'Random_coef': random_enet_final.coef_}
coefs = pd.DataFrame(coefs_data)
coefs.to_csv('HP1_Only_Model_Coefficients.csv', index = False)

# In[32]:


# Simulate predictions
cpr11b_a = pd.read_csv('Cpr11B_HP1a_HP1Only_sim.csv')
cpr11b_b = pd.read_csv('Cpr11B_HP1B_HP1Only_sim.csv')
cpr11b_c = pd.read_csv('Cpr11B_HP1C_HP1Only_sim.csv')
hp1a_predictions = elastic_net.predict(cpr11b_a)
hp1b_predictions = elastic_net.predict(cpr11b_b)
hp1c_predictions = elastic_net.predict(cpr11b_c)

cpr11b_predictions = pd.DataFrame({'HP1a': hp1a_predictions,
                                'HP1B': hp1b_predictions,
                                'HP1C': hp1c_predictions})
cpr11b_predictions.to_csv('Cpr11B_HP1Only_sim_predictions.csv', index = False)

# In[33]:


dgt3_a = pd.read_csv('dgt3_HP1a_HP1Only_sim.csv')
dgt3_b = pd.read_csv('dgt3_HP1B_HP1Only_sim.csv')
dgt3_c = pd.read_csv('dgt3_HP1C_HP1Only_sim.csv')

a_dgt3_preds = elastic_net.predict(dgt3_a)
b_dgt3_preds = elastic_net.predict(dgt3_b)
c_dgt3_preds = elastic_net.predict(dgt3_c)

dgt3_predictions = pd.DataFrame({'HP1a': a_dgt3_preds,
                                'HP1B': b_dgt3_preds,
                                'HP1C': c_dgt3_preds})
dgt3_predictions.to_csv('dgt3_HP1Only_sim_predictions.csv', index = False)

# In[34]:


cecA1_a = pd.read_csv('CecA1_HP1a_HP1Only_sim.csv')
cecA1_b = pd.read_csv('CecA1_HP1B_HP1Only_sim.csv')
cecA1_c = pd.read_csv('CecA1_HP1C_HP1Only_sim.csv')

a_cecA1_preds = elastic_net.predict(cecA1_a)
b_cecA1_preds = elastic_net.predict(cecA1_b)
c_cecA1_preds = elastic_net.predict(cecA1_c)

cecA1_predictions = pd.DataFrame({'HP1a': a_cecA1_preds,
                                 'HP1B': b_cecA1_preds,
                                 'HP1C': c_cecA1_preds})
cecA1_predictions.to_csv('CecA1_HP1Only_sim_predictions.csv', index = False)

# In[35]:


mtk_a = pd.read_csv('Mtk_HP1a_HP1Only_sim.csv')
mtk_b = pd.read_csv('Mtk_HP1B_HP1Only_sim.csv')
mtk_c = pd.read_csv('Mtk_HP1C_HP1Only_sim.csv')

a_mtk_preds = elastic_net.predict(mtk_a)
b_mtk_preds = elastic_net.predict(mtk_b)
c_mtk_preds = elastic_net.predict(mtk_c)

mtk_predictions = pd.DataFrame({'HP1a': a_mtk_preds,
                               'HP1B': b_mtk_preds,
                               'HP1C': c_mtk_preds})
mtk_predictions.to_csv('Mtk_HP1Only_sim_predictions.csv', index = False)

# In[ ]:
