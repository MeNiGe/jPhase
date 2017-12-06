package jphase.fit;


import java.util.LinkedList;
import java.util.ListIterator;

import jphase.ContPhaseVar;
import jphase.CorrHyperErlangVar;
import jphase.DenseContPhaseVar;
import jphase.HyperErlangVar;
import jphase.fit.EMHyperErlangFit;
import jphase.MatrixUtils;
import jphase.Utils;


/**
 * This class implements the Maximum Likelihood method 
 * proposed by Qiu, Perez, Birke, Chen, Harrison", 2016. The method 
 * returns a Hyper-Erlang distribution, a subclass of 
 * Phase-type distributions.  
 * @author Juan F. Pérez
 * @version 1.0
 */
public class EMCHEFit {
    
	double[][] data;
	//double[] data2;
	double[] sumData; // sum of data1  and data2
	double[] avgData; // avg of data1  and data2
	double[] prodData; // product of data1  and data2
	
	int numSamples = 0; // number of samples
	int r = 0; // number of correlated variables
	
	ContPhaseVar var;
	
    
    /**
     * @see jphase.fit.MLContPhaseFitter#MLContPhaseFitter(double[])
     */
    public EMCHEFit(double[][] data){
    	this.data = data;
    	this.r = data.length;
    	this.numSamples = data[0].length;
		
		this.sumData = new double[this.numSamples];
		this.avgData = new double[this.numSamples];
		this.prodData = new double[this.numSamples];
		for ( int i = 0; i < this.numSamples; i++){
			double sumTemp = 0;
			double prodTemp = 1;
			for ( int j = 0; j < this.r ; j++){
				sumTemp += data[j][i];
				prodTemp *= data[j][i];
			}
			this.sumData[i] = sumTemp; 
			this.avgData[i] = sumTemp/this.r;
			this.prodData[i] = prodTemp;
		}
    }
    
    double logLH = 0;
    
    
    /**
     * Precision for other calculation, small nonzero number 
     */
    public static double precision = 1E-50;
    
    /**
     * Precision for the convergence criterion
     * in the algorithm
     */
    public static double convCrit = 1E-6;
    
    /**
     * Precision for the convergence criterion in
     * the coefficient of variance
     */
    public static double precisionCV = 10E-2;
    
    /**
     * Maximum number of iterations for the algorithm execution
     */
    public static int maxIter = 200;

    /**
     * Returns a HyperErlang variable with the best fit, 
     * in the form of a Dense Continuous Phase variable
     * @return HyperErlang variable with the best fit
     */
    //@Override
    public ContPhaseVar fit() {
    	int N = 10; 
        return fit(N); 
    }
    
    /**
     * Returns a CorrHyperErlang variable with the best fit, 
     * in the form of a Dense Continuous Phase variable
     * @param N number of phases in the distribution
     * @return HyperErlang variable with the best fit
     */
    public ContPhaseVar fit(int N) {
        if( this.data!=null ){
	        ContPhaseVar[] vars = new ContPhaseVar[N];
	        double[] LH = new double[N];
	        int best = 0;
	        for(int i = 1; i <= N; i++){
	            CorrHyperErlangVar temp = new CorrHyperErlangVar(N, i, new int[i], 
	                    new int[i], new double[i], new double[i], new double[i], true);
	            double tempLH = doFitNM(temp);
	            if(Math.abs(tempLH)>precision )LH[i-1] = tempLH ;
	            else LH[i-1] = Double.NEGATIVE_INFINITY;
	            vars[i-1] = temp.copy();
	            if(LH[i-1]>LH[best])best=i-1;
	        }
	        this.var = vars[best];
	        this.logLH = LH[best];
	        return this.var;
        }
        return null;
    }
    
    
        
    
    /**
     * This method returns a completely specified HyperErlang
     * variable, such that it has the best likelihood between 
     * all the possible combinations of N phases in M branches  
     * @param var HyperErlang variable with the parameters
     * N and M determined
     * @return Likelihood of the best variable found. The variable
     * is modified with the best parameters found.
     */
    public double doFitNM(CorrHyperErlangVar var){
    	//list of options combining N phases in M branches 
        LinkedList<int[]> listAllocBranch = new LinkedList<int[]>();
        int M = var.getM();
        int N = var.getN();
                
        double[] alphas = new double[M];
        //double[] lambdas = new double[M];
        double tasa = 0; 
        for (int i = 0; i < this.r; i++){
        	tasa += MatrixUtils.average(data[i]);
        }
        tasa /= this.r;
        tasa = 1/tasa; 
        //double tasa = this.r/MatrixUtils.average(sumData);
        double md = M; md=1/md; 
        CorrHyperErlangVar[] vars=null;
        double[] logLH=null;
        
        if (var.getAlphas()[0] == 0){
	        for(int i=0;i<M;i++){
	            alphas[i]=md;
	            //lambdas[i]=tasa;
	        }
        }else{
            alphas = var.getAlphas();
            //lambdas = var.getLambdas();
        }
        
               
        int best = 0;
        LinkedList<CorrHyperErlangVar> varList = new LinkedList<CorrHyperErlangVar>();
        if(M<=N){
            //Initialization 
            int[] a = new int [M];
            int[] b = new int [M];
            a[0]=N-M+1;
            for(int i=1; i<M;i++)a[i]=1;
                
            //Fill-up allocation list (phases to branches) 
            System.arraycopy(a,0,b,0,M);
            listAllocBranch.addLast(a);
            searchAdd(b,listAllocBranch,0);
            
            //Fill-up full allocation list (phases to paths)
            ListIterator<int[]> iter = listAllocBranch.listIterator();
            
            
            //for
            while(iter.hasNext()){ // traverse all phase-to-branch allocations
            	LinkedList<int[]>[] listAllocPerBranch = new LinkedList[M];
            	int[] thisAlloc = iter.next(); 
            	for(int j = 0; j < M; j ++){ // traverse each branch
            		listAllocPerBranch[j] = new LinkedList<int[]>(); 
            		for(int k = 1; k <= thisAlloc[j]; k++){ // consider all possible numbers of paths in this branch given the allocation of phases
            			// create allocation of thisAlloc[j] phases to k paths 
            			int[] aB = new int [k];
                        int[] bB = new int [k];
                        aB[0] = thisAlloc[j] - k + 1;
                        for(int i=1; i<k;i++)aB[i]=1;
                        System.arraycopy(aB,0,bB,0,k);
                        listAllocPerBranch[j].addLast(aB);
                        searchAdd(bB,listAllocPerBranch[j],0);
            		}
            	}
            	
            	// allocs for all branches have been obtained, now combine them into a list of variables
            	int[] numAllocs = new int[M]; 
            	for (int i = 0; i < M; i ++)
            		numAllocs[i] = listAllocPerBranch[i].size();
            	
            	int[] idxAlloc = new int[M]; // index to traverse all allocation
            	int currDim = M-1; //current dimension where analysis is being performed 
            	boolean ready = false; 
            	while (!ready){
            		//create variable for current index
            		// number of paths in each branch
            		int[] V = new int[M]; 
            		for(int j = 0; j < M; j++){
            			V[j] = listAllocPerBranch[j].get(idxAlloc[j]).length; 
            		}
            		           		
            		int sumV = 0;
            		for(int j = 0; j < M; j++)
            			sumV += V[j]; 
            		
            		// number of phases in each path
            		int[] r = new int[sumV]; 
            		int idxP = 0; 
            		for(int j = 0; j < M; j ++){
            			for(int i = 0; i < V[j]; i ++){ 
            				r[idxP] = listAllocPerBranch[j].get(idxAlloc[j])[i];
            				idxP++; 
            			}
                	}
            		
            		// path selection probabilities in each branch
            		double[] betas = new double[sumV];
            		idxP = 0; 
            		for(int j = 0; j < M; j ++){
            			for(int i = 0; i < V[j]; i ++){ 
            				betas[idxP] = 1/((double) V[j]); 
            				idxP++; 
            			}
                	}
            		
            		// rates for each path
            		double[] lambdas = new double[sumV];
            		idxP = 0; 
            		//int initExp = -(int)Math.floor(((double)sumV)/2);
            		int initExpB = -(int)Math.floor(((double)M)/2);
            		for(int j = 0; j < M; j ++){
            			int initExpP = -(int)Math.floor(((double)V[j])/2);
            			for(int i = 0; i < V[j]; i ++){ 
            				//lambdas[idxP] = tasa*Math.pow(10,initExp+idxP); 
            				lambdas[idxP] = tasa*Math.pow(10,initExpB+j)*Math.pow(V[j], initExpP+i);
            				idxP++; 
            			}
                	}
            		
            		// add variable to list 
            		varList.addLast(new CorrHyperErlangVar(N, M, V, r, alphas, betas, lambdas, true));
            		
            		
            		// advance index
            		if(idxAlloc[M-1] < numAllocs[M-1]-1 ){ // last dimension still in exploration
            			idxAlloc[M-1]++; 
            		}else{ // last dimension exhausted, check which idx to update
            			int j = M-2;
            			while(j > -1 && idxAlloc[j] == numAllocs[j]-1) j--; // first index with values to explore
            			if(j>=0){// first index non-negative
            				idxAlloc[j]++;
            				for(int i = j+1; i<M; i++) // reset all indexes in front
            					idxAlloc[i] = 0; 
            			}else // all indexes have been covered
            				ready = true; 
            		}
            			
            	}
            	
            }
            
            
                    
            //Lists of logLH and random variables 
            int numVar = varList.size(); 
            logLH = new double[numVar];
            //vars = new CorrHyperErlangVar[listAllocBranch.size()];
            
            
            ListIterator<CorrHyperErlangVar> iterVar = varList.listIterator();
            int k=0;
            while(iterVar.hasNext()){
                //int[] temp = iter.next();
                
                //Initialization
                /*for(int i=0;i<M;i++){
                	lambdas[i]=temp[i]*tasa+Math.pow(10,i-2);
                }*/
                //HyperErlangVar varTemp = new HyperErlangVar(var.getN(),var.getM(),temp,alphas,lambdas,true);
                //CorrHyperErlangVar varTemp = new CorrHyperErlangVar(var.getN(),var.getM(),temp,alphas,lambdas,true);
                
                //Fit
            	CorrHyperErlangVar thisVar = iterVar.next();
            	System.out.print("Fitting configuration "+k+" of "+numVar+" - N: "+N+" - M: "+M);
            	System.out.print("Vb: [");
            	for(int i = 0; i < M; i++){
            		System.out.print(" "+thisVar.getV()[i]);
            	}
            	System.out.print("] - hbj: [");
            	for(int i = 0; i < thisVar.getSumV(); i++){
            		System.out.print(" "+thisVar.getR()[i]);
            	}
            	System.out.print("]");
            	Long tini = System.currentTimeMillis(); 
                logLH[k]=doFitNMR(thisVar);
                System.out.println("\nExec time (s): "+((double)(System.currentTimeMillis()-tini))/1000);
                System.out.println("LogLH: "+ logLH[k]);
                System.out.println("Mean: "+thisVar.moment(1));
                //System.out.println("Mean: "+thisVar.toString());
                //logLH[k]=0;
                //vars[k]=varTemp;
                if(logLH[k]>=logLH[best])best=k;
                k++;
            }
        }
        var.setV(varList.get(best).getV());
        var.setSumV(varList.get(best).getSumV());
        var.setR(varList.get(best).getR());
        var.setAlphas(varList.get(best).getAlphas());
        var.setBetas(varList.get(best).getBetas());
        var.setLambdas(varList.get(best).getLambdas()); 
        
        System.out.print("\nBest configuration: N: "+N+" M: "+M);
    	System.out.print("Vb: [");
    	for(int i = 0; i < M; i++){
    		System.out.print(" "+var.getV()[i]);
    	}
    	System.out.print("]\nhbj: [");
    	for(int i = 0; i < var.getSumV(); i++){
    		System.out.print(" "+var.getR()[i]);
    	}
    	System.out.print("]\n");
        
        
        return logLH[best];
    }

    
    /**
     * This method searches all the possible combinations of integers
     * in vector a, such that they always sum to the same number, and
     * there is not redundant options. 
     * @param a array in which the modifications are made
     * @param lista list where the new combinations found are stored
     * @param k index of the position being examined
     */
    private void searchAdd(int[] a, LinkedList<int[]> lista, int k) {
        if(k<a.length-1){
            while(a[k]-1 >= a[k+1]+1){
                a[k]--;
                a[k+1]++;
                int[] b = new int[a.length];
                System.arraycopy(a,0,b,0,a.length);
                lista.addLast(b);
                searchAdd(a,lista,k+1); 
            }
        }
    }
    

    /**
     * This method returns a completely specified CorrHyperErlang
     * variable, such that it has the best likelihood after the
     * execution of the EM algorithm for the case where the 
     * variable has N phases in M branches, distributed as determined
     * by the vectors V and r   
     * experiments to be fitted
     * @param var HyperErlang variable with the parameters
     * N, M and r determined
     * @return Likelihood of the best variable found. The variable
     * is modified with the best parameters found.
     */
    public double doFitNMR(CorrHyperErlangVar var){
        double logLHprev =0;
        double logLHnow =0;
        int n = var.getN();
        int M=var.getM();
        int sumV=var.getSumV();
        int[] V = var.getV();
        //int k=sumData.length;
        int iter = 0;
        
        /*
        if(m == 1){
        	// 1 branch -> Hyper Erlang case 
            var.setAlphas(new double[]{1.0});
            EMHyperErlangFit fit2 = EMHyperErlangFit(this.data);
            HyperErlangVar var2 = HyperErlangFit(this.data); 
            
            double lambda = (n)/FitterUtils.powerMomentK(avgData, 1);
            var.setLambdas(new double[] {lambda});
            for(int i =0;i < this.numSamples; i++ )logLHnow += 2*n*Math.log(lambda)
                    + (n-1)*Math.log(prodData[i]) -lambda*sumData[i]
                    -2*Math.log( Utils.fact(n-1));
                                                       
        }else{*/
            while( (Math.abs(logLHnow-logLHprev)>convCrit || iter ==0) && iter <= maxIter ){
                iter++;
                //compute Pm
                double[][][] pm = new double[sumV][this.numSamples][this.r];		// p(xdl | b, j, hat Theta)
                
                //double[][] pm2 = new double[M][this.numSamples];					// p(xd  | b, hat Theta)
                double[][] pmB = new double[M][this.numSamples];					// p(xd  | b, hat Theta)
                for(int i = 0; i<sumV; i++){
                    for(int j = 0; j < this.numSamples; j++){
                    	for(int l = 0; l < this.r; l++){
	                        double temp = var.getR()[i]*Math.log(var.getLambdas()[i]) -   	// r_i ln(lambda_i) 
	                        Utils.lnFactorial(var.getR()[i]-1) - 							// - ln((n_i-1)!)
	                        var.getLambdas()[i]*data[l][j] +  								// - lambda_i( x_jl )
	                        (var.getR()[i]-1)*Math.log(data[l][j]);						  	// + (r_i-1)log( x_jl)
	                        
	                        //pm[i][j]=var.getLambdas()[i]*Math.exp(temp);
	                        pm[i][j][l]=Math.exp(temp);
                    	}
                    }
                }
                
                int idxP = 0;
                for(int i = 0; i<M; i++){
                	for(int k = 0; k<this.numSamples; k++){
                		pmB[i][k] = 1;
                		for(int l = 0; l < this.r; l++){
                			double temp = 0;
                			for(int j = 0; j < V[i]; j++){
                				temp += var.getBetas()[idxP+j] * pm [idxP+j][k][l]; 
                			}
                			pmB[i][k] *= temp; 
                		}
                		
                	}
                	idxP += V[i]; 
                }
        
                //E-STEP
                //double[][][] qm = new double[sumV][this.numSamples][this.r]; 	// delta(b, j | x_dl, hat Theta) 
                //double[][] sumaPm = new double[this.numSamples][this.r];	    // p(xdl  | hat Theta)
                
                double[][][] qm2 = new double[sumV][this.numSamples][this.r]; 	// delta(j | b, x_dl, hat Theta)
                double[][][] sumaPm2 = new double[this.numSamples][this.r][M];	// p(xdl  | hat Theta, b)
                double[][] qmB = new double[M][this.numSamples]; 				// delta(b | x_dl,    hat Theta)
                double[] sumaPmB = new double[this.numSamples];					// p(xd  | hat Theta)
                
                //init Version 
                /*
                double logLH = 0;
                for(int k = 0; k < this.numSamples; k++){
                	for(int l = 0; l < this.r; l++){
                		idxP = 0; //idxPath 
                		for(int i = 0; i < M; i++){
                			for(int j = 0; j < V[i]; j++){
                				sumaPm[k][l] += var.getAlphas()[i]*var.getBetas()[idxP]*pm[idxP][k][l];
                				idxP++; 
                			}
                			if(sumaPm[k][l]>precision)
                    			logLH+=Math.log(sumaPm[k][l]);
                		}
                		
                	}
                }
                
                
                
                
                idxP = 0; 
                for(int i = 0; i<M; i++){
                	for(int j = 0; j < V[i]; j++){
	                    for(int k = 0; k<this.numSamples; k++){
	                    	for(int l = 0; l < this.r; l++){
	                    		if(sumaPm[k][l] > precision)
	                    			qm[idxP][k][l]= var.getAlphas()[i]*var.getBetas()[idxP]*pm[idxP][k][l]/sumaPm[k][l];
	                    	}
	                    }
	                    idxP++; 
                	}
                }
                */
                
                
                // new Version 
                // sumaPm2
                //double logLH = 0;
                for(int k = 0; k < this.numSamples; k++){
                	for(int l = 0; l < this.r; l++){
                		idxP = 0; //idxPath 
                		for(int i = 0; i < M; i++){
                			for(int j = 0; j < V[i]; j++){
                				sumaPm2[k][l][i] += var.getBetas()[idxP]*pm[idxP][k][l];
                				idxP++; 
                			}
                			//if(sumaPm2[k][l][i]>precision)
                    		//	logLH+=Math.log(sumaPm2[k][l][i]);
                		}
                		
                	}
                }
                
                // sumaPmB
                double logLH = 0;
                for(int k = 0; k < this.numSamples; k++){
                	for(int i = 0; i < M; i++){
                		sumaPmB[k] += var.getAlphas()[i] * pmB[i][k]; 
                		if(sumaPmB[k]>precision)
                			logLH+=Math.log(sumaPmB[k]);
                	}
                }
                
                // qm2
                idxP = 0; 
                for(int i = 0; i<M; i++){
                	for(int j = 0; j < V[i]; j++){
	                    for(int k = 0; k < this.numSamples; k++){
	                    	for(int l = 0; l < this.r; l++){
	                    		if(sumaPm2[k][l][i] > precision)
	                    			qm2[idxP][k][l]= var.getBetas()[idxP]*pm[idxP][k][l]/sumaPm2[k][l][i];
	                    	}
	                    }
	                    idxP++; 
                	}
                }
                
                
                // qmB
                for(int i = 0; i<M; i++){
                	for(int k = 0; k<this.numSamples; k++){
                		if(sumaPmB[k] > precision)
                			qmB[i][k] = var.getAlphas()[i] * pmB[i][k] / sumaPmB[k];
                	}
            	}
                				
                
                
                
                
                
                    
                //M-STEP
                //mStep(qm, var.getR(), var.getAlphas(), var.getLambdas(), var);
                mStep(qm2, qmB, var);
                
                logLHprev=logLHnow;
                logLHnow=logLH;
                        
            }
        //}
    return logLHnow;
    }
    
    
    
    
    /**
     * This method calculates the new set of parameters 
     * through the maximization Step and assign them to the
     * parameter variable (HyperErlang) 
     * experiments to be fitted
     * @param Qm qm values (HyperErlang densities) from the 
     * Expectation step 
     * @param r distribution of the phases in the Erlang branches
     * @param alphas mass probabilities in the previous iteration
     * @param lambdas rates of the Erlang branches in the previous 
     * iteration  
     * @param var HyperErlang distribution to be characterized 
     * by the new set of parameters (result of the maximization Step) 
     */
    private void mStep(double[][][] Qm2, double[][] QmB, CorrHyperErlangVar var){
	//private void mStep(double[][][] Qm, double[][][] Qm2, double[][] QmB, CorrHyperErlangVar var){
    	double[] alphas = var.getAlphas();
    	int[] r = var.getR();
    	int[] V = var.getV();
    	double[] betas = var.getBetas();
    	double[] lambdas = var.getLambdas(); 
    	
    	
        int M = var.getM();
        int K = this.numSamples; //sumData.length;
        int sumV = var.getSumV(); 
        
        //qm = new double[sumV][this.numSamples][this.r];
        //qmB = new double[M][this.numSamples];
        //qm2 = new double[sumV][this.numSamples][this.r]; 	// delta(j | b, x_dl, hat Theta)
        
        //compute alpha
        //first compute phi[d,b,l] = sum_j Qm[b,j][d][l] Q is oh size sumV, phi is of size M
        /*
        double[][][] phi = new double[M][K][this.r];
        for(int j = 0; j < K; j++){
        	for(int l = 0; l < this.r; l++){
        		int idxB = 0; 
        		for(int i = 0; i < M; i++){
            		for(int k = 0; k < V[i]; k++){
            			phi[i][j][l] += Qm[idxB+k][j][l]; 
            		}
            		idxB += V[i]; 
            	}        
        	}
        }
        
        // obtain alpha as product of phi        
        for(int i = 0; i < M; i++){
        	alphas[i] = 0;
            for(int k = 0; k < K; k++){
            	double temp = 1;
            	for(int l = 0; l < this.r; l++){
            		temp *= phi[i][k][l];
            	}
                alphas[i]+=temp;
            }
            //lambdas[i]=alphas[i];
            alphas[i]= alphas[i]/K;
        }
        */
        
        //second way to compute alphas 
        for(int i = 0; i < M; i++){
        	alphas[i] = 0;
        	for(int k = 0; k < K; k++){
        		alphas[i] += QmB[i][k];
        	}
        	alphas[i] /= K; 
        }
        normalizeVector(alphas, 0, M); 
        
        
        
        //compute betas and lambdas - first sums 
        /*
        double[] sumaQm = new double[sumV];  // sum_d sum_l Qm[][d][l]
        double[] sumaQmx = new double[sumV]; // sum_d sum_l Qm[][d][l] * x[d][l]
        
        for(int i = 0; i < sumV; i++){
            for(int k = 0; k < K; k++){
            		for(int l = 0; l < this.r; l++){
            			sumaQm[i] += Qm[i][k][l];
            			sumaQmx[i] += Qm[i][k][l]*data[l][k];
            		}
            }
        }
        
        double[] sumaQmTot = new double[M]; // sum_d sum_l sum_j\inVb Qm[j][d][l]
        int idxB = 0; 
		for(int i = 0; i < M; i++){
    		for(int k = 0; k < V[i]; k++){
    			sumaQmTot[i] += sumaQm[idxB+k]; 
    		}
    		idxB += V[i]; 
    	}
        
		//compute beta
		idxB = 0; 
		for(int i = 0; i < M; i++){
			if(sumaQm[i]>precision){
	    		for(int k = 0; k < V[i]; k++){
	    			betas[idxB+k] = sumaQm[idxB+k]/sumaQm[i]; 
	    		}
			}
    		idxB += V[i]; 
    	}
        
		//compute lambda
        for(int i = 0; i < sumV; i++){
            if(sumaQm[i]>precision) lambdas[i]= r[i]*sumaQm[i]/sumaQmx[i];
        }
        var.setAlphas(alphas);
        var.setBetas(betas);
        var.setLambdas(lambdas);
        */
        
        
        // second way to compute betas and lambdas - first sums
        double[] sumaQm2 = new double[sumV];  // sum_d sum_l Qm[][d][l]
        double[] sumaQmx2 = new double[sumV]; // sum_d sum_l Qm[][d][l] * x[d][l]
        
        for(int i = 0; i < sumV; i++){
            for(int k = 0; k < K; k++){
            		for(int l = 0; l < this.r; l++){
            			sumaQm2[i] += Qm2[i][k][l];
            			sumaQmx2[i] += Qm2[i][k][l]*data[l][k];
            		}
            }
        }
        
        double[] sumaQmTot2 = new double[M]; // sum_d sum_l sum_j\inVb Qm[j][d][l]
        int idxB = 0; 
		for(int i = 0; i < M; i++){
    		for(int k = 0; k < V[i]; k++){
    			sumaQmTot2[i] += sumaQm2[idxB+k]; 
    		}
    		idxB += V[i]; 
    	}
        
		//compute beta
		idxB = 0; 
		for(int i = 0; i < M; i++){
			if(sumaQm2[i]>precision){
	    		for(int k = 0; k < V[i]; k++){
	    			betas[idxB+k] = sumaQm2[idxB+k]/sumaQmTot2[i]; 
	    		}
	    		normalizeVector(betas, idxB, V[i]); 
			}
    		idxB += V[i]; 
    	}
		
		
		//compute lambda
        for(int i = 0; i < sumV; i++){
            if(sumaQm2[i]>precision) lambdas[i]= r[i]*sumaQm2[i]/sumaQmx2[i];
        }
        var.setAlphas(alphas);
        var.setBetas(betas);
        var.setLambdas(lambdas);
        
        
    }
    
    void normalizeVector(double[] vector, int ini, int num){
    	double sum = 0;
    	for(int i = ini; i < ini+num; i++)
    		sum += vector[i];
    	if( sum > 0){
	    	for(int i=ini; i < ini+num; i++)
	    		vector[i] /= sum;
    	}
    }
}