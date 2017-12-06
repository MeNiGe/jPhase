package jphase.fit;


import java.util.LinkedList;
import java.util.ListIterator;

import jphase.DenseContPhaseVar;
import jphase.ContPhaseVar;
import jphase.HyperErlangVar;
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
public class EMHyperErlangPairedFit {
    
	double[] data1;
	double[] data2;
	double[] sumData; // sum of data1  and data2
	double[] avgData; // avg of data1  and data2
	double[] prodData; // product of data1  and data2
	
	ContPhaseVar var;
	
    
    /**
     * @see jphase.fit.MLContPhaseFitter#MLContPhaseFitter(double[])
     */
    public EMHyperErlangPairedFit(double[] data1, double[] data2){
    	this.data1 = data1;
		this.data2 = data2;
		
		this.sumData = new double[data1.length];
		this.avgData = new double[data1.length];
		this.prodData = new double[data1.length];
		for ( int i = 0; i < data1.length; i++){
			this.sumData[i] = data1[i] + data2[i];
			this.avgData[i] = this.sumData[i]/2;
			this.prodData[i] = data1[i]*data2[i];
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
        return DenseContPhaseVar.HyperErlang(doFitHyperErlang());
    }
    
    /**
     * Returns a HyperErlang variable with the best fit, 
     * in the form of a Dense Continuous Phase variable
     * @param N number of phases in the distribution
     * @return HyperErlang variable with the best fit
     */
    public ContPhaseVar fit(int N) {
        if( data1!=null && data2!=null){
	        ContPhaseVar[] vars = new ContPhaseVar[N];
	        double[] LH = new double[N];
	        int best = 0;
	        for(int i = 1; i <= N; i++){
	            HyperErlangVar temp = new HyperErlangVar(N, i,new int[i], 
	                    new double[i], new double[i], true);
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
    public HyperErlangVar doFitHyperErlang() {
        //double CVC = MatrixUtils.CV(data);
        double CVC = MatrixUtils.CV(sumData);
        double CVT = 0;
        boolean ready = false;
        int N = 1;
        int M = 1;
        HyperErlangVar res = null;
        double logLH =0;
        int k=0;
        while(!ready){
            k++;
            System.out.println("\n\nITERATION: "+k);
            System.out.println("N: "+N+"\t M: "+M);
            res = new HyperErlangVar(N,M,new int[M], new double[M], new double[M], true);
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
    public double doFitNM(HyperErlangVar var){
    	//list of options combining N phases in M branches 
        LinkedList<int[]> listaPosib = new LinkedList<int[]>();
        int M = var.getM();
        int N = var.getN();
                
        double[] alphas = new double[M];
        double[] lambdas = new double[M];
        //double tasa = 1.0/MatrixUtils.average(data);
        double tasa = 2.0/MatrixUtils.average(sumData);
        double md = M; md=1/md; 
        HyperErlangVar[] vars=null;
        double[] logLH=null;
        
        if (var.getLambdas()[0] == 0){
	        for(int i=0;i<M;i++){
	            alphas[i]=md;
	            lambdas[i]=tasa;
	        }
        }else{
            alphas = var.getAlphas();
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
            vars = new HyperErlangVar[listaPosib.size()];
            
            
            ListIterator<int[]> iter = listaPosib.listIterator();
            int k=0;
            while(iter.hasNext()){
                int[] temp = (int[])iter.next();
                
                //Initialization
                for(int i=0;i<M;i++){
                	lambdas[i]=temp[i]*tasa+Math.pow(10,i-2);
                }
                HyperErlangVar varTemp = new HyperErlangVar(var.getN(),var.getM(),temp,alphas,lambdas,true);
                
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
    public double doFitNMR(HyperErlangVar var){
        double logLHprev =0;
        double logLHnow =0;
        int n = var.getN();
        int m=var.getM();
        int k=sumData.length;
        int iter = 0;
        
        if(m == 1){
        	// 1 branch -> Erlang case 
            var.setAlphas(new double[]{1.0});
            double lambda = (n)/FitterUtils.powerMomentK(avgData, 1);
            var.setLambdas(new double[] {lambda});
            for(int i =0;i < k; i++ )logLHnow += 2*n*Math.log(lambda)
                    + (n-1)*Math.log(data1[i]*data2[i]) -lambda*sumData[i]
                    -2*Math.log( Utils.fact(n-1));
                                                       
        }else{
            while( (Math.abs(logLHnow-logLHprev)>precision || iter ==0) && iter <= maxIter ){
                iter++;
                //compute Pm
                double[][] pm = new double[m][k];
                for(int i = 0; i<m; i++){
                    for(int j = 0; j<k; j++){
                        double temp = 2*var.getR()[i]*Math.log(var.getLambdas()[i]) -   //2 r_i ln(lambda_i) 
                        2*Utils.lnFactorial(var.getR()[i]-1) - 							// - 2ln((r_i-1)!)
                        var.getLambdas()[i]*sumData[j] +  								// - lambda_i(x_n+y_n )
                        (var.getR()[i]-1)*Math.log(prodData[j]);						// + (r_i-1)log(x_ny_n)
                        
                        //pm[i][j]=var.getLambdas()[i]*Math.exp(temp);
                        pm[i][j]=Math.exp(temp);
                        
                    }
                }
        
                //E-STEP
                double[][] qm = new double[m][k];
                double[] sumaPm = new double[k];
                
                double logLH = 0;
                for(int j = 0; j < k; j++){
                    for(int i = 0; i<m; i++){
                        sumaPm[j] += var.getAlphas()[i]*pm[i][j];
                    }
                    if(sumaPm[j]>precision)
                        logLH+=Math.log(sumaPm[j]);
                }
                
                for(int i = 0; i<m; i++){
                    for(int j = 0; j<k; j++){
                        if(sumaPm[j]>precision)
                            qm[i][j]= var.getAlphas()[i]*pm[i][j]/sumaPm[j];
                    }
                }
                    
                //M-STEP
                mStep(qm, var.getR(), var.getAlphas(), var.getLambdas(), var);
                
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
    private void mStep(double[][] Qm, int[] r, double[] alphas, double[] lambdas, HyperErlangVar var){
        int m=r.length;
        int k=sumData.length;
        double[] sumaQm = new double[m];
        //computing alpha
        for(int i = 0; i < m; i++){
            for(int j = 0; j < k; j++){
                alphas[i]+=Qm[i][j];
            }
            lambdas[i]=alphas[i];
            alphas[i]= alphas[i]/k;
        }
        
        //computing lambda 
        for(int i = 0; i < m; i++){
            for(int j = 0; j < k; j++){
                sumaQm[i] += Qm[i][j]*sumData[j];
            }
        }
        
        for(int i = 0; i<m; i++){
            if(sumaQm[i]>precision)lambdas[i]= 2*r[i]*lambdas[i]/sumaQm[i];
        }
        var.setAlphas(alphas);
        var.setLambdas(lambdas);
    }
}