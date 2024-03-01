Instructions for running the software
 
The following steps are required:
	Fill in the file Input.xlsx (the file provided includes the data for the presented use case)
	Run the function RecOptim.m
	Open the file Output.xlsx to see the results.
Note that, every time the software is run, the Matlab file Inputs_MAT.mat is rewritten with the simulation results. In order to obtain the same results of the presented use case, it is required to include the original values of the presented use case in the file Inputs_MAT.mat. These original data are recorded in the file Inputs_MAT_0.mat which is included in the software folder. 
Note that, after the first interval of simulation has run, before running the following time interval, the file Input.xlsx must be updated as follows:
	The sheet RES must be updated with the new values of the forecasts and the operating period can be changed.
	The sheet Load must be updated with the new values of the forecasts and the operating period can be changed.
	The sheet DR_Load_Type_I must be updated according to the new load requirements and, if present, new aggregation of loads of this type.
	The sheet DR_Load_Type_II must be updated by deleting the load data related to the preceeding time interval and adding, if present, new aggregation of loads of this type.
	The sheet DR_Load_Type_III must be updated according to the new load requirements and, if present, new aggregation of loads of this type.
This procedure is based on the assumption that DR loads of type I and III can be rescheduled during the operation period whereas Dr loads of type II, once scheduled, cannot be modified for their own characteristics. 
