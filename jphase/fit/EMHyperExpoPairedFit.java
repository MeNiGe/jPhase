package jphase.fit;


import jphase.DenseContPhaseVar;


/**
 * This class implements the Maximum Likelihood method 
 * proposed by Khayari, Sadre and Haverkort in "Fitting
 * world-wide web request traces with the EM algorithm", 2003.
 * The method returns a Hyper-Exponential distribution, a 
 * subclass of Phase-type distributions.
 * @author Juan F. Pérez
 * @version 1.0 
 */
public class EMHyperExpoPairedFit{
	
	double[] data1;
	double[] data2;
	double[] sumData; // sum of data1  and data2
	
	DenseContPhaseVar var;
	
	/**
	 * @see jphase.fit.MLContPhaseFitter#MLContPhaseFitter(double[])
	 */
	public EMHyperExpoPairedFit(double[] data1, double[] data2){
		this.data1 = data1;
		this.data2 = data2;
	}
	
	double logLH = 0;
	
	
	/**
	 * Precision for the convergence criterion
	 * in the algorithm
	 */
	public static double precision = 10E-5;
	
	/**
	 * Precision for the convergence criterion in
	 * the coefficient of variance
	 */
	public static double precisionParam = 10E-3;
	
	/**
	 * Number of exponential phases in the hyper-exponential
	 * distribution in evaluation 
	 */
	private int n;
	
	
	/**
	 * Probability of choose each exponential phase in the 
	 * hyper-exponential distribution in evaluation 
	 */
	private double[] probs;
	
	/**
	 * Rate at each exponential phase in the 
	 * hyper-exponential distribution in evaluation 
	 */
	private double[] rates;
	
	
	
	//@Override
	/**
	 * Main fitting method. Uses data and number of phases 
	 * already defined. 
	 */
	public DenseContPhaseVar fit() {
		DenseContPhaseVar var = new DenseContPhaseVar();
		System.out.println("\n\nFitting a hyper-exponential distribution with "+this.n+" phases");
		double LH = fit_main();
		var = DenseContPhaseVar.HyperExpo(this.rates, this.probs);
		System.out.println("Distribution found\n: "+var.toString());
		System.out.println("Likelihood: "+LH);
		System.out.println("Mean: "+var.expectedValue());
		System.out.println("Variance: "+var.variance());
		System.out.println("CV: "+var.CV());
		this.var = var;
		return var;
	}
	
	
	/**
	 * 
	 * @return HyperExponential variable with the best fit
	 */
	public DenseContPhaseVar autofit() {
		int N = 15;
		int best = 0;
		double[] LH = new double[N]; 
		DenseContPhaseVar[] vars = new DenseContPhaseVar[N];
		for(int k=0; k < N; k++){
			System.out.println("\n\nITERATION: "+(k+1));
			System.out.println("N: "+(k+1));
			this.n=k+1;
			vars[k] = new DenseContPhaseVar(k+1);
			//LH[k] = doFitN(data);
			LH[k] = fit_main();
			vars[k] = DenseContPhaseVar.HyperExpo(rates, probs);
			System.out.println("Likelihood: "+LH[k]);
			System.out.println("Mean: "+vars[k].expectedValue());
			System.out.println("Variance: "+vars[k].variance());
			System.out.println("CV: "+vars[k].CV());
			if(LH[k]>LH[best])best=k;
		}
		System.out.println("The best variable found is: "
				+vars[best].toString());
		System.out.println("Likelihood: "+LH[best]);
		return vars[best];
	}
	
	/**
	 * Main fitting method. Uses data and number of phases 
	 * already defined. 
	 */
	public DenseContPhaseVar fit(int k) {
		DenseContPhaseVar var = new DenseContPhaseVar();
		if (k > 0){
			this.n = k;
			System.out.println("\n\nFitting a hyper-exponential distribution with "+this.n+" phases");
			double LH = fit_main();
			var = DenseContPhaseVar.HyperExpo(this.rates, this.probs);
			System.out.println("Distribution found:\n "+var.toString());
			System.out.println("Likelihood: "+LH);
			System.out.println("Mean: "+var.expectedValue());
			System.out.println("Variance: "+var.variance());
			System.out.println("CV: "+var.CV());
			this.var = var;
			this.logLH = LH;
		}else{
			System.out.println("The number of phases must be a positive integer");
		}
		return var;
	}
	
	
	
	/**
	 * This method implements the EM algorithm from the data, with 
	 * the specified number of exponential phases
	 * @param data non-negative data trace from independent experiments
	 * @return The loglikelihood of the solution found. The values of the
	 * parameters are stored in the attributes probs and rates
	 */
	private double fit_main(){
	//public double doFitN(double[] data){
		//Initialization
		probs = new double[n];
		rates = new double[n];
		for(int i = 0; i < this.n; i++){
			probs[i] = 1.0/this.n;
			//probs[i] = (0.9 * Math.pow(10, -i ))/(1-Math.pow(10, -this.n));
			//rates[i] = Math.pow(10, -i+1);
		}
		
		sumData = new double[data1.length];
		for (int i = 0; i < data1.length; i++){
			sumData[i] = data1[i] + data2[i]; 
		}
		
		double C = 0;
		//double mean = FitterUtils.powerMomentK(data1, 1);
		double mean = (FitterUtils.powerMomentK(data1, 1) + FitterUtils.powerMomentK(data2, 1))/2;
		if (this.n % 2 == 0){
			for(int i = 0; i < this.n/2; i++){
				rates[i] = Math.pow(10, -n/2+i)/mean;
				C += Math.pow(10, -n/2+i);
			}
			for(int i = this.n/2; i < this.n; i++){
				rates[i] = Math.pow(10, -n/2+i+1)/mean;
				C += Math.pow(10, -n/2+i+1);
			}
		}else{
			for(int i = 0; i < this.n; i++){
				rates[i] = Math.pow(10, (this.n-1)/2+i)/mean;
				C += Math.pow(10, (this.n-1)/2+i);
			}
		}
		C = ((double)this.n)/C;
		for(int i = 0; i < this.n; i++){
			rates[i] = rates[i]*C;
		}
		
		//EM algorithm
		int iter = 0;
		int maxIter = 300;
		boolean ready = false;
		double[][] p = new double[sumData.length][this.n];
		double LHold = 0;
		double LHnew = 0;
		while(ready == false){
			iter ++;
			eStep(sumData, p);
			LHold = LHnew;
			LHnew = mStep(sumData, p);
			double dif = Math.abs(LHold-LHnew);
			if(dif < precision || iter > maxIter )ready = true;
		}
		return LHnew;
	}

	
	/**
	 * Executes the Expectation step of the EM algorithm
	 * @param data non-negative data trace from independent experiments
	 * @param p likelihood estimate to be computed
	 */
	private void eStep(double[] sumData, double[][]p){
		for(int n = 0; n < sumData.length; n++){
			double denom =0.0;
			for(int i = 0; i < this.n; i++){
				 //double num = probs[i]*rates[i]*Math.exp(-rates[i]*sumData[n]);
				 double num = probs[i]*rates[i]*rates[i]*Math.exp(-rates[i]*sumData[n]); //p_i * lambda_i^2 * exp(lambda_i(x+y)
				 p[n][i]= num;
				 denom += num;
			}
			for(int i = 0; i < this.n; i++)p[n][i] /= denom;
		}

	}
	
	/**
	 * Executes the Maximization step of the EM algorithm
	 * @param data non-negative data trace from independent experiments
	 * @param p likelihood estimate calculated at the Expectation step
	 * @return log-likelihood of the resulting parameter set
	 */
	private double mStep(double[] sumData, double[][] p){
		int N = sumData.length;
		for(int i = 0; i < this.n; i++){
			double sum = 0;
			double sumX = 0;
			for(int n = 0; n < N; n++){
				sum += p[n][i];
				sumX += p[n][i]*sumData[n];
			}
			probs[i] = sum/N;
			rates[i] = 2*sum/sumX;
		}
		double LH=0; 
		for(int n = 0; n < N; n++){
			double sumLH = 0;
			for(int i = 0; i < this.n; i++){
				//sumLH += probs[i]*rates[i]*Math.exp(-rates[i]*data1[n]);
				sumLH += probs[i]*rates[i]*rates[i]*Math.exp(-rates[i]*sumData[n]); //p_i * lambda_i^2 * exp(lambda_i(x+y)
			}
			LH += Math.log(sumLH); 
		}
			return LH;
	}
	
	public double getLogLikelihood(){
		return this.logLH;
	}
}
