package jphase.fit;


import java.util.LinkedList;
import java.util.ListIterator;

import jphase.DenseContPhaseVar;
import jphase.ContPhaseVar;
import jphase.CorrProbHyperErlangVar;
import jphase.MatrixUtils;
import jphase.Utils;


/**
 * This class implements the Maximum Likelihood method 
 * proposed by Thümmler, Buchholz and Telek in "A novel approach
 * for fitting probability distributions to real trace
 * data with the EM algorithm", 2005. The method 
 * returns a Hyper-Erlang distribution, a subclass of 
 * Phase-type distributions.  
 * @author Juan F. Pérez
 * @version 1.0
 */
public class EMCorrProbHyperErlangFit {
    
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
    public EMCorrProbHyperErlangFit(double[][] data){
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
     * Precision for the convergence criterion
     * in the algorithm
     */
    public static double precision = 10E-50;
    
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
        //return DenseContPhaseVar.CorrProbHyperErlang(doFitHyperErlang());
        return doFitHyperErlang();
    }
    
    /**
     * Returns a HyperErlang variable with the best fit, 
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
	            CorrProbHyperErlangVar temp = new CorrProbHyperErlangVar(N, i,new int[i], 
	                    new double[i], new double[i][i], new double[i], true);
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
     * Returns a HyperErlang variable with the best fit
     * experiments to be fitted
     * @return HyperErlang variable with the best fit
     */
    public CorrProbHyperErlangVar doFitHyperErlang() {
        //double CVC = MatrixUtils.CV(data);
        double CVC = MatrixUtils.CV(sumData);
        double CVT = 0;
        boolean ready = false;
        int N = 1;
        int M = 1;
        CorrProbHyperErlangVar res = null;
        double logLH =0;
        int k=0;
        while(!ready){
            k++;
            System.out.println("\n\nITERATION: "+k);
            System.out.println("N: "+N+"\t M: "+M);
            res = new CorrProbHyperErlangVar(N,M,new int[M], new double[M], new double [M][M], new double[M], true);
            //CorrProbHyperErlangVar(int N, int M, int[] r, double[] alphas, double[][] betas, double[] lambdas, boolean deep){
            logLH = doFitNM(res);
            CVT = res.CV();
            System.out.println("SCV (data): "+CVC+"\t SCV (variable): "+CVT);
            System.out.println("Mean: "+res.expectedValue());
            System.out.println("Variance: "+res.variance());
            if(Math.abs(CVT-CVC)<precisionCV)ready=true;
            else if(CVC>CVT){
                M=M+1;
                N=Math.max(M,N);
            }
            else{
                N=N+1;
            }
            if(k>100)ready = true;
        }
        System.out.println("Final LogLH: "+logLH);
        System.out.println(res.toString());
        return res;
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
    public double doFitNM(CorrProbHyperErlangVar var){
    	//list of options combining N phases in M branches 
        LinkedList<int[]> listaPosib = new LinkedList<int[]>();
        int M = var.getM();
        int N = var.getN();
                
        double[] alphas = new double[M];
        double[][] betas = new double[M][M];
        double[] lambdas = new double[M];
        //double tasa = 1.0/MatrixUtils.average(data);
        double tasa = this.r/MatrixUtils.average(sumData);
        double md = M; md=1/md; 
        CorrProbHyperErlangVar[] vars=null;
        double[] logLH=null;
        
        if (var.getLambdas()[0] == 0){
	        for(int i=0;i<M;i++){
	            alphas[i]=md;
	            lambdas[i]=tasa;
	            for(int j=0;j<M;j++){
	            	betas[i][j] = md;
	            }
	        }
        }else{
            alphas = var.getAlphas();
            betas = var.getBetas();
            lambdas = var.getLambdas();
        }
        
               
        int best = 0;
        
        if(M<=N){
            //Initialization 
            int[] a = new int [M];
            int[] b = new int [M];
            a[0]=N-M+1;
            for(int i=1; i<M;i++)a[i]=1;
                
            //Fill-up option list
            System.arraycopy(a,0,b,0,M);
            listaPosib.addLast(a);
            searchAdd(b,listaPosib,0);
                    
            //Lists of lohLH and random variables 
            logLH = new double[listaPosib.size()];
            vars = new CorrProbHyperErlangVar[listaPosib.size()];
            
            
            ListIterator<int[]> iter = listaPosib.listIterator();
            int k=0;
            while(iter.hasNext()){
                int[] temp = (int[])iter.next();
                
                //Initialization
                /*for(int i=0;i<M;i++){
                	lambdas[i]=temp[i]*tasa+Math.pow(10,i-2);
                }*/
                CorrProbHyperErlangVar varTemp = new CorrProbHyperErlangVar(var.getN(),var.getM(),temp,alphas,betas, lambdas,true);
                
                //Fit
                logLH[k]=doFitNMR(varTemp);
                vars[k]=varTemp;
                if(logLH[k]>=logLH[best])best=k;
                k++;
            }
        }
        var.setR(vars[best].getR());
        var.setAlphas(vars[best].getAlphas());
        var.setLambdas(vars[best].getLambdas());
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
     * This method returns a completely specified HyperErlang
     * variable, such that it has the best likelihood after the
     * execution of the EM algorithm for the case where the 
     * variable has N phases in M branches, distributed as determined
     * by the vector r   
     * experiments to be fitted
     * @param var HyperErlang variable with the parameters
     * N, M and r determined
     * @return Likelihood of the best variable found. The variable
     * is modified with the best parameters found.
     */
    public double doFitNMR(CorrProbHyperErlangVar var){
        double logLHprev =0;
        double logLHnow =0;
        int n = var.getN();
        int m=var.getM();
        //int k=sumData.length;
        int iter = 0;
        
        if(m == 1){
        	// 1 branch -> Erlang case ***TBD
            var.setAlphas(new double[]{1.0});
            double lambda = (n)/FitterUtils.powerMomentK(avgData, 1);
            var.setLambdas(new double[] {lambda});
            for(int i =0;i < this.numSamples; i++ )logLHnow += 2*n*Math.log(lambda)
                    + (n-1)*Math.log(prodData[i]) -lambda*sumData[i]
                    -2*Math.log( Utils.fact(n-1));
                                                       
        }else{
            while( (Math.abs(logLHnow-logLHprev)>precision || iter ==0) && iter <= maxIter ){
                iter++;
                //compute Pm
                double[][][] pm = new double[m][this.r][this.numSamples]; //pm[i][j][d]: prob of observing x_dj given that the j-th replica of the d-th obs takes branch i
                for(int i = 0; i<m; i++){
                	for(int v = 0; v<this.r; v++){
	                    for(int d = 0; d<this.numSamples; d++){
	                        double temp = var.getR()[i]*Math.log(var.getLambdas()[i]*data[v][d]) -  // r_i ln(lambda_i*x_dv) 
	                        Utils.lnFactorial(var.getR()[i]-1) - 									// - ln((r_i-1)!)
	                        var.getLambdas()[i]*data[v][d];   										// - lambda_i x_dv
	                        //Math.log(lambdaprodData[j]);											// + log( lambda_i)
	                        
	                        pm[i][v][d]=var.getLambdas()[i]*Math.exp(temp);
	                        //pm[i][j]=Math.exp(temp);
	                        
	                    }
                	}
                }
        
                //E-STEP
                double[][] qm = new double[m][this.numSamples]; 
                double[][][][] qm2 = new double[m][m][this.r-1][this.numSamples];
                double[] sumaPm = new double[this.numSamples]; 						// sum_i alpha_i * pm[i][1][d]
                double[][][] sumaPm2 = new double[m][this.r-1][this.numSamples]; 	// sum_j beta_ij * pm[i][j][d]
                
                double logLH = 0;
                for(int d = 0; d < this.numSamples; d++){
                	double tempi = 0;
                    for(int i = 0; i<m; i++){
                    	double temp = var.getAlphas()[i]*pm[i][0][d];  // first replica of obs d takes branch i
                    	sumaPm[d] += temp; 

	                    for(int v = 1; v<this.r; v++){
	                    	double tempSum = 0; 
	                    	for(int j = 0; j<m; j++){
	                    		tempSum += var.getBetas()[i][j]*pm[j][v][d]; // v-th replica of obs d takes branch j given that first rep took branch i    
	                    	}
	                    	sumaPm2[i][v-1][d] += tempSum;  // v-th replica of obs d takes branch j given that first rep took branch i
	                    	temp *= tempSum; 
	                    }
	                    tempi += temp; 
                    }
                    
                    
                    if(sumaPm[d]>precision)
                        logLH+=Math.log(tempi);
                }
                
                //qm
                for(int i = 0; i<m; i++){
                    for(int d = 0; d<this.numSamples; d++){
                        if(sumaPm[d]>precision)
                            qm[i][d]= var.getAlphas()[i]*pm[i][0][d]/sumaPm[d];
                    }
                }
                
                //qm2
                for(int i = 0; i<m; i++){
                	for(int j = 0; j<m; j++){
	                    for(int d = 0; d<this.numSamples; d++){
	                    	for(int v = 1; v<this.r; v++){
		                        if(sumaPm2[i][v-1][d]>precision)
		                            qm2[i][j][v-1][d]= var.getBetas()[i][j]*pm[j][v][d]/sumaPm2[i][v-1][d];
	                    	}
	                    }
                	}
                }
                    
                //M-STEP
                mStep(qm, qm2, var.getR(), var.getAlphas(), var.getBetas(), var.getLambdas(), var);
                //mStep(qm, qm2, var.getR(), System.arraycopy(var.getAlphas(),0,new double[m],0,m), var.getBetas(), var.getLambdas(), var);
                
                
                logLHprev=logLHnow;
                logLHnow=logLH;
                        
            }
        }
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
    private void mStep(double[][] Qm, double[][][][] Qm2, int[] r, double[] alphas, double[][] betas, double[] lambdas, CorrProbHyperErlangVar var){
        int m = r.length;
        int D = this.numSamples;
        double[] sumaQm = new double[m];
        double[] sumaQm2i = new double[m];
        //computing alpha
        for(int i = 0; i < m; i++){
        	alphas[i] = 0;
            for(int d = 0; d < D; d++){
                alphas[i]+=Qm[i][d];
            }
            lambdas[i]=alphas[i];
            alphas[i]= alphas[i]/D;
        }
        
        //computing betas
        for(int i = 0; i < m; i++){
        	for(int j = 0; j < m; j++){
        		betas[i][j] = 0;
        		for(int v = 1; v < this.r; v++){
        			for(int d = 0; d < D; d++){
        				betas[i][j] += Qm[i][d] * Qm2[i][j][v-1][d];
        				sumaQm2i[j] +=  Qm2[i][j][v-1][d]; //prob this obs takes branch j   
        			}
	            }
        		//lambdas[i]=alphas[i];
                betas[i][j] = betas[i][j]/( (this.r-1) * lambdas[i] );
        	}
        }
        
        //computing lambda 
        // denominator
        for(int i = 0; i < m; i++){
            for(int d = 0; d < D; d++){
                sumaQm[i] += Qm[i][d]*data[0][d];
                for(int j = 0; j < m; j++){
                	for(int v = 1; v < this.r; v++){
                		sumaQm[i] += Qm2[j][i][v-1][d]*data[v][d];
                	}
                }
            }
        }
        
        for(int i = 0; i<m; i++){
            if(sumaQm[i]>precision)lambdas[i]= r[i]*(lambdas[i] + sumaQm2i[i] )/sumaQm[i];
        }
        var.setAlphas(alphas);
        var.setBetas(betas);
        var.setLambdas(lambdas);
    }
}