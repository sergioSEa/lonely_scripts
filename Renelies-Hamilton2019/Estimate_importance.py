# Sergio Andre-Sanchez 12/02/2021
# Prediction of taxa presenting vertical transmision and interpretabiloity of the featurs that determine it
# Includes: two versions of the data with and without normalization
import numpy as np
import pandas as pd
import shap
import sklearn
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.model_selection import train_test_split
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from joblib import dump, load
import random
random.seed( 30 )


shap.initjs()

###Dataset with normalized abundances --> PD2

#If pickle does not exist, uncomment this part
#PD2  = pd.read_csv("DataFrame2.tsv", sep="\t")
#PD2.to_pickle("DataFrame2.pkl")
#If pickle does not exist, comment this part
PD2 = pd.read_pickle("./DataFrame2.pkl")

#Remove NA (mainly in the network columns) and standarize column name to use all the functions the same between the two DF
PD2 = PD2.dropna(axis=0)
PD2["Abundance_in_origin"] = PD2["Ab_normalized"]
PD2 = PD2.drop(["Ab_normalized"], axis=1)

#Comment in/out if pickle available

#PD = pd.read_csv("DataFrame.tsv", sep="\t") 
PD = pd.read_pickle("./DataFrame.pkl")
#PD.to_pickle("DataFrame.pkl")
PD = PD.dropna(axis=0)


def Preprocess(PD):
	'Make dummy variables of the taxonomy features'
	#Make dummy variables
	PD = PD.drop(["ASV", "Genus", "Presence_in_origin"], axis=1)

	Family_DF = pd.get_dummies(PD["Family"])
	Order_DF = pd.get_dummies(PD["Order"])
	Class_DF = pd.get_dummies(PD["Class"])
	Phylum_DF = pd.get_dummies(PD["Phylum"])
	Tax_DF = pd.get_dummies(PD["Full_taxonomy"])
	Lin_DF =  pd.get_dummies(PD["Lineage"])

	PD = PD.drop(["Family","Order","Class","Phylum","Full_taxonomy", "Lineage"], axis=1)

	for i in [Family_DF,Order_DF,Class_DF,Phylum_DF,Tax_DF, Lin_DF]:
		PD = pd.concat([PD, i], axis=1)
	return(PD)
PD = Preprocess(PD)
PD2 = Preprocess(PD2)

#1. Divide data in training and test
X = PD.drop(["VT","VT_random"], axis=1)
y = PD["VT"].astype(int)
y2 = PD["VT_random"].astype(int)


X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=7)


#PD2 = PD2.drop(["o_NA", "f_NA"], axis=1)
X_new = PD2.drop(["VT","VT_random"], axis=1)
y_new = PD2["VT"].astype(int)
y2_new = PD2["VT_random"].astype(int)


X_train_new, X_test_new, y_train_new, y_test_new = train_test_split(X_new, y_new, test_size=0.2, random_state=7)



#2. Apply models: first time I tested different models, after discovering that RF works best I just use RF
n_folds = 10

def Explain_model(Model, X_train):
	'Get shap in linear models'
	explainer = shap.Explainer(Model, X_train, feature_names=X_train.columns)
	shap_values = explainer(X_train)
	#Plot
	fig = shap.plots.waterfall(shap_values[1], show=False)
	plt.savefig('scratch.png')
	
	#Importance values
	vals= np.abs(shap_values[1].values).mean(0)
	feature_importance = pd.DataFrame(list(zip(X_train.columns,vals)),columns=['col_name','feature_importance_vals'])
	feature_importance.sort_values(by=['feature_importance_vals'],ascending=False,inplace=True)
	print(feature_importance)
def Explain_kernel(Model,  X_train):
	'Get shap in any model - really slow'
	#https://slundberg.github.io/shap/notebooks/Iris%20classification%20with%20scikit-learn.html
	explainer =shap.KernelExplainer(Model.predict_proba, X_train)
	shap_values = explainer.shap_values(X_test)
	fig = shap.plots.waterfall(shap_values[1], show=False)
	plt.savefig('scratch.png')
	
def Explain_Tree(Model, X_train, Prefix):
	'Get shap in tree models'
	#Prefix will tell: Complete, Step1, Step2  + Random/Regular
	Output_summary = "Plots/{P}_summaryPlot_Shap.pdf".format(P=Prefix)
	Output_bar = "Plots/{P}_barPlot_Shap.pdf".format(P=Prefix)
	Output_waterfall =  "Plots/{P}_waterfallPlot_Shap.pdf".format(P=Prefix)
	Output_table = "Tables/{P}_meanShap.pdf".format(P=Prefix)
	Completa_shap = "Tables/Complete_shap/{P}_shap.tsv".format(P=Prefix)
	Completa_data = "Tables/Complete_shap/{P}_data.tsv".format(P=Prefix)
	print("Calculating SHAP values")
	#https://shap.readthedocs.io/en/latest/example_notebooks/tabular_examples/model_agnostic/Diabetes%20regression.html?highlight=random%20forest#Random-forest
	explainer = shap.TreeExplainer(Model,check_additivity=False)
	shap_values = explainer.shap_values(X_test)	
	print("Drawing figures")
	plt.clf()
	fig = shap.summary_plot(shap_values[1], X_test, show=False, plot_size= (40,40) )
	plt.savefig(Output_summary)
	plt.clf()
	fig2 = shap.summary_plot(shap_values[1], X_test, plot_type="bar", show=False, plot_size= (40,40))
	plt.savefig(Output_bar)
	plt.clf()
	#try: 
	#	fig3 = shap.plots.waterfall(shap_values[1], show=False)
	#	plt.savefig(Output_waterfall)
	#	plt.clf()
	#except:
	#	print("Waterfall plot failed")

	print("Getting shap values")
	vals= np.abs(shap_values[1]).mean(0)
	np.savetxt(Completa_data, shap_values[1], delimiter="\t")
	X_test.to_csv(Completa_shap, sep= "\t")
	feature_importance = pd.DataFrame(list(zip(X_train.columns,vals)),columns=['col_name','feature_importance_vals'])
	feature_importance.sort_values(by=['feature_importance_vals'],ascending=False,inplace=True)
	feature_importance.to_csv(Output_table, sep = "\t")

def Do_lasso(X_train,y_train, n_folds):
	#Lostic regression penalty l1 means LASSO, solver is the algorithm  used for optimization, max_iter is the number of iterations it goes on
	# Set values of the grid search; C: float, default=1.0. Inverse of regularization strength, the smaller the C the stronger the regularization. 
	C_values = [0.001, 0.01, 0.1, 1, 10, 100, 1000]
	C_grid = {'C': C_values}
	Logistic_model = sklearn.linear_model.LogisticRegression(penalty="l1",solver='saga', tol=0.01, max_iter=150)
	grid_logReg = sklearn.model_selection.GridSearchCV(Logistic_model, C_grid, cv=n_folds, refit=True,n_jobs= -1) #scoring: method to score, n_jobs: jobs in paralel
	grid_logReg.fit(X_train,y_train)

	CV_results = pd.DataFrame( list(zip(C_values, grid_logReg.cv_results_['mean_test_score'])), columns=["C" , "CV"])
	CV_results = CV_results.sort_values(by="CV",ascending=False)

	best_logReg = grid_logReg.best_estimator_
	Accuracy =  best_logReg.score(X_train,y_train)
	return(Accuracy, best_logReg)
	#Logistic_model.predict() #Gives you the label to predict
	#Logistic_model.predict_proba() #Gives you the probability on the logit

	#Coeff = Logistic_model.coef_.transpose().flatten()
	#logReg_coeff = pd.DataFrame({'feature_name': full_col_names, 'model_coefficient': Logistic_model.coef_.transpose().flatten()})
	#logReg_coeff = logReg_coeff.sort_values('model_coefficient',ascending=False)
	#print(logReg_coeff)

def KNN(X_train,y_train, n_folds):
	Grid_params = {'n_neighbors' : [3, 5, 11, 19] }
	knn = sklearn.neighbors.KNeighborsClassifier(weights = "uniform", metric = 'euclidean' )
	gs = sklearn.model_selection.GridSearchCV(knn, Grid_params,  cv=n_folds, refit=True,n_jobs= -1)
	gs.fit(X_train,y_train)
		
	CV_results = pd.DataFrame( list(zip( [3, 5, 11, 19], gs.cv_results_['mean_test_score'])), columns=["NN" , "CV"])
	CV_results = CV_results.sort_values(by="CV",ascending=False)
	print(CV_results)
	
	best = gs.best_estimator_
	Accuracy =  best.score(X_train,y_train)
	return(Accuracy, best)

def RF(X_train,y_train, n_folds):
	n_features = X_train.shape[1]
	n_features = round(n_features**(1/2))
	n_F = [n_features-2,n_features-1,n_features, n_features+1,n_features+2]
	param_grid = { 'bootstrap': [True],
	'max_features': n_F,
	}
	rf = sklearn.ensemble.RandomForestClassifier(n_estimators = 500)
	grid_search = sklearn.model_selection.GridSearchCV(estimator = rf, param_grid = param_grid, cv = n_folds, n_jobs= -1) #n_jobs= -1
	grid_search.fit(X_train,y_train)
	
	CV_results = pd.DataFrame( list(zip( n_F, grid_search.cv_results_['mean_test_score'])), columns=["Features" , "CV"])
	CV_results = CV_results.sort_values(by="CV",ascending=False)
	print(CV_results)

	best = grid_search.best_estimator_
	Accuracy =  best.score(X_train,y_train)

	return(Accuracy, best)


def Choice(X_train,y_train, n_folds):
	'Function that compares KNN, RF and lasso performance'
	Results = {}
	print("Running RF")
	Acc, RaFo = RF(X_train,y_train, n_folds)
	Results["RF"]  = Acc
	dump(RaFo, 'Models/RF.joblib') 

	print("Running KNN")
	Acc, KNN = KNN(X_train,y_train, n_folds)
	Results["KNN"]  = Acc
	dump(KNN, 'Models/KNN.joblib')

	print("Running Logistic regression")
	Acc, Logistic_model = Do_lasso(X_train,y_train, n_folds)
	Results["Regression"]  = Acc
	dump(Logistic_model, 'Models/Logistic.joblib')
	print(Results)


#We performed Choice on the complete real dataset and decided to use random forest as the algorithm through the work


def Prediction(Model_choice, X_test, y_test):
	'Predict labels and compute 1. Accuracy, 2. Confusion matrix' 
	Yhat = Model_choice.predict(X_test)
	#Yhat= []
	#[ Yhat.append(1) if X>= 0.5 else Yhat.append(0) for X in yhat ]
	Accuracy = sklearn.metrics.accuracy_score(y_test,Yhat)
	conf_mat = sklearn.metrics.confusion_matrix(y_test, Yhat)
	return(Yhat, Accuracy, conf_mat)
def Gini(X_train, Model_choice, Prefix):
	'Get Gini values from random forest model'
	Output = "Tables/{Prefix}_Gini.tsv".format(Prefix=Prefix)
	Gini_importance = pd.DataFrame(list(zip(X_train.columns ,Model_choice.feature_importances_)), columns = ["Feature", "Gini"])
	print("Saving Gini importance of RF model")
	Gini_importance.to_csv(Output, sep="\t")
def Produce_ROC(X_test, y_test, Model_choice, Prefix):
	'Compute AUC, get values for ROC curve, make plot'
	AUC_figure = "Plots/{Prefix}_ROC.pdf".format(Prefix=Prefix)
	Data_save = "Tables/{Prefix}_ROC.tsv".format(Prefix=Prefix)

	yhat = Model_choice.predict_proba(X_test)[:,1]
	
	fpr, tpr, threshold =  sklearn.metrics.roc_curve(y_test, yhat)
	roc_auc = sklearn.metrics.auc(fpr, tpr)
	plt.title('Receiver Operating Characteristic')
	plt.plot(fpr, tpr, 'b', label = 'AUC = %0.2f' % roc_auc)
	plt.legend(loc = 'lower right')
	plt.plot([0, 1], [0, 1],'r--')
	plt.xlim([0, 1])
	plt.ylim([0, 1])
	plt.ylabel('True Positive Rate')
	plt.xlabel('False Positive Rate')
	plt.savefig('AUC.png')
	plt.clf()
		
	Data_ROC = pd.DataFrame(list(zip(fpr, tpr)), columns= ["fpr", "tpr"])
	Data_ROC["Model"] = [Prefix] * len(fpr)
	Data_ROC.to_csv(Data_save, sep= "\t")
	print("AUC: " + str(roc_auc))

def Complete_real(X_train, X_test, y_test):
	'This was a preliminary funciton just for complete real data, use Do_models instead'
	print("=========================COMPLETE REAL DATA====================")
	print("loading data")
	Model_choice =  load('Models/RF.joblib') 
	print("compute metrics")
	Yhat, Accuracy, conf_mat =  Prediction(Model_choice, X_test, y_test)
	# true negatives is C00, false negatives is C1,0 , true positives is C1,1  and false positives C0,1
	print("Confusion matrix of prediction")
	print(conf_mat)
	print("Accuracy of prediction: " + str(Accuracy))
	print("Extracting Gini values")
	Gini(X_train, Model_choice, "Complete")
	print("Estimating AUC")
	Produce_ROC(X_test, y_test, Model_choice, "Complete")
	print("Shap values")
	Explain_Tree(Model_choice, X_train, "Complete")

def Do_models(X, y, Prefix, n_folds=10, Do_model = True, Do_shap=True):
	'Main function to adjust models  and compute variable importance and accuracy metrics'
	print("========================={p}=====================".format(p=Prefix))
	X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=7)
	if Do_model == True:
		print("Running RF")
		Acc, RaFo = RF(X_train,y_train, n_folds)
		dump(RaFo, 'Models/RF_{p}.joblib'.format(p=Prefix))
		Model_choice =RaFo
	else:
		 Model_choice =  load('Models/RF_{p}.joblib'.format(p=Prefix))
	print("compute metrics")
	Yhat, Accuracy, conf_mat =  Prediction(Model_choice, X_test, y_test)
	np.savetxt('Tables/Confusion_matrix_{p}.csv'.format(p=Prefix), conf_mat.astype(int), delimiter=',',fmt='%i')
	return()
	print("Accuracy of prediction: " + str(Accuracy))
	print("Extracting Gini values")
	Gini(X_train, Model_choice, Prefix)
	Produce_ROC(X_test, y_test, Model_choice, Prefix)
	if Do_shap== True:
		print("Shap values")
		Explain_Tree(Model_choice, X_train, Prefix)

#Complete_real(X_train, X_test, y_test)

#Complete
Do_models(X_new, y_new, "Complete_new", n_folds, Do_model = False, Do_shap= True)
#Random complete
Do_models(X_new, y2_new, "Complete_random_new", n_folds, Do_model = False, Do_shap= True)
#Step1
PD_Step1 = PD2.loc[PD2['Step'] == 1]
X_Step1 = PD_Step1.drop(["VT","VT_random"], axis=1)
y_real = PD_Step1["VT"].astype(int)
y_random = PD_Step1["VT_random"].astype(int)
#Real
Do_models(X_Step1, y_real, "step1_new", n_folds, Do_model = False,Do_shap= True)
#Random
Do_models(X_Step1, y_random, "step1_random_new", n_folds, Do_model = False, Do_shap= True)
#Step2
PD_Step2 = PD2.loc[PD2['Step'] == 2]
X_Step2 = PD_Step2.drop(["VT","VT_random"], axis=1)
y_real = PD_Step2["VT"].astype(int)
y_random = PD_Step2["VT_random"].astype(int)
#Real
Do_models(X_Step2, y_real, "step2_new", n_folds, Do_model = False, Do_shap= True)
#Random
Do_models(X_Step2, y_random, "step2_random_new", n_folds, Do_model = False, Do_shap= True)

#Bash command to put all ROC stats togeter - then use an R script for plotting
#from subprocess import call
#call("cat Tables/*new*_ROC.tsv > ALL_ROC_new", shell=True)
#call("grep -v fpr ALL_ROC_new  > Tables/ALL_ROC_new.tsv", shell=True)



#Complete
Do_models(X, y, "Complete", n_folds, Do_model = False)

#Random complete
Do_models(X, y2, "Complete_random", n_folds, Do_model = False)


#Step1
PD_Step1 = PD.loc[PD['Step'] == 1]
X_Step1 = PD_Step1.drop(["VT","VT_random"], axis=1)

y_real = PD_Step1["VT"].astype(int)
y_random = PD_Step1["VT_random"].astype(int)

#Real
Do_models(X_Step1, y_real, "step1", n_folds, Do_model = False)
#Random
Do_models(X_Step1, y_random, "step1_random", n_folds, Do_model = False)


#Step2
PD_Step2 = PD.loc[PD['Step'] == 2]
X_Step2 = PD_Step2.drop(["VT","VT_random"], axis=1)

y_real = PD_Step2["VT"].astype(int)
y_random = PD_Step2["VT_random"].astype(int)

#Real
Do_models(X_Step2, y_real, "step2", n_folds, Do_model = False)
#Random
Do_models(X_Step2, y_random, "step2_random", n_folds, Do_model = False)





#from subprocess import call
#call("cat Tables/*_ROC.tsv > ALL_ROC", shell=True)
#call("grep -v fpr ALL_ROC  > Tables/ALL_ROC.tsv", shell=True)
