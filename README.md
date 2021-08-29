EMT (epithelial-mesenchymal transition) is a basic developmental process, which might promote cancer metastasis, has been studied from various perspectives. This is a matlab implementation of simulating EMT by MAP theory.

matrix.txt is the matrix representing the relationship between genes.

EMT_action.m is used to calculate the path action

EMTode_plus.m is used to solve the stable states of the ODEs.
	 
EMT_refine.m is used to refine the transition path. It does not affect the transition path, but only beautifies it.
	 
EMT_landscape_MAP.m is the main function which is uesd to establish the energy landscape of EMT and calculate the minimum action path (MAP) of EMT.
	 
Please run EMT_landscape_MAP.m and the program runs about 30 hours. Reducing the parameter time1/time in EMT_landscape_MAP.m (Line 7,8) will reduce the running time, but the quality of the transition path will be reduced.

