package jphase;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;


/**
 * This class allows the creation and manipulation 
 * of HyperErlang distributions. The associated matrix has  
 * dense representation 
 * @author Juan F. Perez
 * @version 1.0
 */
public class CorrHyperErlangVar  extends AbstractContPhaseVar implements PhaseVar {
    
    /**
     * Total Number of phases in the HyperErlang distribution 
     */
    private int N;
    
    /**
     * Number of branches 
     */
    private int M;
    
    /**
     * Number of Erlang paths in each branch
     */
    private int[] V;
    
    /**
     * Total Number of Erlang paths 
     */
    private int sumV;
        
    /**
     * Number of phases on each branch and Erlang path
     */
    private int[] r;
    
    /**
     * Probability mass function of choosing each one of the branches 
     */
    private double[] alphas;
    
    /**
     * Probability mass function of choosing each one of the path in each branch 
     */
    private double[] betas;
    
    /**
     * Rates of each Erlang path
     */
    private double[] lambdas;
    
    
    /**
     * Constructor of a Correlated Hyper Erlang variable in dense representation.
     * As default it has just one branch that is taken with probability
     * one. the unique branch has one phase with rate 1 per time unit. 
     */
    public CorrHyperErlangVar() {
        this.N=1;
        this.M=1;
        this.r=new int[this.M];
        this.r[0]=1;
        this.V = new int[this.M]; 
        this.V[0]=1;
        this.sumV = 1; 
        this.alphas=new double[this.M];
        this.alphas[0]=1;
        this.betas=new double[this.sumV];
        this.betas[0]=1;
        this.lambdas=new double[this.sumV];
        this.lambdas[0]=1;
    }
    
    /**
     * Constructor of a Hyper Erlang variable with n phases
     * in dense representation
     * @param n Total number of phases
     */
    public CorrHyperErlangVar(int n){
        this.N=n;
    }
    
    /**
     * Constructor of a Hyper Erlang variable in dense representation
     * @param N Total number of phases
     * @param M Number of branches
     * @param V number of paths in each branch
     * @param r Number of phases in each path
     * @param alphas Probability associated to each branch
     * @param betas Probability associated to each path 
     * @param lambdas Rate associated to each branch
     * @param deep True if this is a deep constructor, false if not
     */
    public CorrHyperErlangVar(int N, int M, int[] V, int[] r, double[] alphas, double[] betas, double[] lambdas, boolean deep){
        if(M==V.length && M==alphas.length && r.length==lambdas.length && r.length==betas.length){
            this.N=N;
            this.M=M;
            if(deep==true){
            	this.V = new int[M];
            	this.alphas=new double[M];
            	System.arraycopy(V,0,this.V,0,M);
            	System.arraycopy(alphas,0,this.alphas,0,M);
            	this.updateSumV();
            	
                this.r=new int[this.sumV];
                this.betas=new double[this.sumV];
                this.lambdas=new double[this.sumV];
                System.arraycopy(r,0,this.r,0,this.sumV);
                System.arraycopy(lambdas,0,this.lambdas,0,this.sumV);
                System.arraycopy(betas,0,this.betas,0,this.sumV);
            }else{
                this.r=r;
                this.V=V;
                this.updateSumV();
                this.alphas=alphas;
                this.betas=betas;
                this.lambdas=lambdas;
            }
        }else{
        	System.out.println("Parameters are not valid for Correlated Hyper Erlang distribution\nVariable not initialized.");
        }
        	
    }
    
    
    
    /**
     * Constructor of a Hyper Erlang variable in dense representation
     * @param V number of paths in each branch
     * @param r Number of phases in each branch
     * @param alphas Probability associated to each branch 
     * @param lambdas Rate associated to each branch
     * @param deep True if this is a deep constructor, false if not
     */
    public CorrHyperErlangVar(int[] V, int[] r, double[] alphas, double[] betas, double[] lambdas, boolean deep) {
        if(V.length==alphas.length && r.length==lambdas.length && r.length==betas.length){
            this.M=alphas.length;
            
            if(deep==true){
            	this.V = new int[M];
            	this.alphas=new double[M];
            	System.arraycopy(V,0,this.V,0,M);
            	System.arraycopy(alphas,0,this.alphas,0,M);
            	this.updateSumV();
            	
                this.r=new int[M];
                this.betas=new double[this.sumV];
                this.lambdas=new double[M];
                System.arraycopy(r,0,this.r,0,this.sumV);
                System.arraycopy(betas,0,this.betas,0,this.sumV);
                System.arraycopy(lambdas,0,this.lambdas,0,this.sumV);
            }else{
                this.r=r;
                this.V=V;
                this.updateSumV();
                this.alphas=alphas;
                this.betas=betas;
                this.lambdas=lambdas;
            }
            
            // total # of phases
            for(int i=0;i<this.sumV;i++){
                this.N+=this.r[i];
            }
        }else{
        	System.out.println("Parameters are not valid for Correlated Hyper Erlang distribution\nVariable not initialized.");
        }
    }
    
    
    /**
     * @return Total number of phases
     */
    public int getN(){
        return this.N;
    }
    
    /**
     * @param N Total number of phases to set
     */
    public void setN(int N){
        this.N=N;
    }
    
    /**
     * @return Number of branches
     */
    public int getM(){
        return this.M;
    }
    
    /**
     * @param M Number of branches to set
     */
    public void setM(int M){
        this.M=M;
    }
    
    /**
     * @return Number of Erlang paths in each branch
     */
    public int[] getV(){
        return this.V;
    }
    
    /**
     * @param V Number of Erlang paths in each branch
     */
    public void setV(int[] V){
        this.V=V;
    }
    
    
    /**
     * @return Total number of Erlang paths ( in all branches)
     */
    public int getSumV(){
        return this.sumV;
    }
    
    /**
     * @param sumV Total number of Erlang paths ( in all branches)
     */
    public void setSumV(int sumV){
        this.sumV=sumV;
    }
    
    
    /**
     * @return Number of phases in each branch
     */
    public int[] getR(){
        return this.r;
    }
    
    /**
     * @param r Number of phases in each branch to set
     */
    public void setR(int[] r){
        this.r=r;
    }
    

    /**
     * @return Probability associated to each branch
     */
    public double[] getAlphas(){
        return this.alphas;
    }
    
    /**
     * @param alphas Probability associated to each branch to set
     */
    public void setAlphas(double[] alphas){
        this.alphas=alphas;
    }
    
    /**
     * @return Probability of each path 
     */
    public double[] getBetas(){
        return this.betas;
    }
    
    /**
     * 
     * @param betas Probability of each path
     */
    public void setBetas(double[] betas){
        this.betas=betas;
    }
    
    /**
     * @return Rate associated to each branch
     */
    public double[] getLambdas(){
        return this.lambdas;
    }
    
    /**
     * @param lambdas Rates associated to each branch to set
     */
    public void setLambdas(double[] lambdas){
        this.lambdas=lambdas;
    }

	public Matrix getMatrix() {
        return new DenseMatrix(getDMatrix());        
    }
    
    /**
     * Returns the Double MAtrix that represents the variable
     * @return Double Matrix
     */
    public double[][] getDMatrix(){
        int N = this.getN();
        int[] V = this.getV(); 
        int[] h = this.getR();
        
        int idxPh = 0; //index phase
        int idxB = 0;  //index branch
        double[][] matrix = new double[N][N];
        for (int i = 0; i < this.getM(); i++) {
            for (int j = 0; j < V[i]; j++) {
            	for (int k = 0; k < h[idxB]-1; k++) {
	                matrix[idxPh][idxPh] = -this.getLambdas()[idxB];
	                matrix[idxPh][idxPh + 1] = this.getLambdas()[idxB];
	                idxPh ++;
            	}
                // final phase in the path
	            matrix[idxPh][idxPh] = -this.getLambdas()[idxB];
                idxPh++; 
                idxB++;
            }
        }

        return matrix;
    }

	public void setMatrix(Matrix A) {
        throw new UnsupportedOperationException();
    }

	public Vector getVector() {
        return new DenseVector(getDVector());
    }

	
    
    /**
     * Returns the Double MAtrix that represents the variable
     * @return Double Vector
     */
    public double[] getDVector(){
        int N = this.getN();
        int[] V = this.getV(); 
        int[] h = this.getR();
        
        double[] vector = new double[N];
//        System.arraycopy(this.getAlphas(),0,vector,0,this.getM());
        
        int idxPh = 0; //index phase
        int idxB = 0;  //index branch
        for (int i = 0; i < this.getM(); i++) {
            for (int j = 0; j < V[i]; j++) {
            	vector[idxPh] = this.getAlphas()[i] * this.getBetas()[idxB];	
                idxPh += h[idxB];
                idxB ++; 
            }
        }
        return vector;

//        double[] vector = new double[N];
//        System.arraycopy(this.getAlphas(),0,vector,0,this.getM());
//        return vector
    }

	public void setVector(Vector alpha) {
        throw new UnsupportedOperationException();
    }

	
	private void updateSumV(){
		int sumV = 0;
		for(int i = 0; i < this.M; i++){
			sumV += this.V[i];
		}
		this.sumV = sumV; 
	}
    
	public ContPhaseVar copy() {
        return new CorrHyperErlangVar(
                this.N, this.M, this.V, this.r,
                this.alphas, this.betas, this.lambdas, true
        );
    }
    
	public ContPhaseVar newVar(int n){
        return new CorrHyperErlangVar(n);
    }
    
    @Override
    public int getNumPhases(){
        return this.getN();
    }
   
        
    @Override
    public double expectedValue() {
        return moment(1); 
    }

    
    @Override
    public double moment(int k) {
        double sum = 0;
        int idxP = 0;
        for(int b = 0; b< this.M; b++ ){
        	for(int i = 0; i< this.V[b]; i++ ){
	            double temp1 = this.alphas[b]*this.betas[idxP+i]/(Math.pow(this.lambdas[idxP+i],k));
	            for(int j = 0; j <=k-1;j++)temp1*=(this.r[idxP+i]+j);
	            sum+=temp1;
        	}
        	idxP += V[b]; 
        }
        return sum;
    }
    
    /*
    @Override
    public String description(){
        String s = "__________________________________________________\n";
        s += "HyperErlang Distribution";
        s+= "\nNumber of Branches: "+this.M+"\n";
        s += "Number of Phases: "+this.N+"\n";
        s=s+"\tBranch\t Ri\t alphai\t\t lambdai\n";
        //for(int i = 0 ; i < this.M; i++)s+=""+(i+1)+"\t "+this.r[i]+"\t"+this.alphas[i]+"\t"+this.lambdas[i]+"\n";
        for(int i = 0 ; i < this.M; i++){
            s+="\t"+(i+1)+"\t "+this.r[i];
            s+=String.format("\t %5.4f",this.alphas[i]);
            s+=String.format("\t\t %6.4f",this.lambdas[i]);
            s+="\n";
        }
        s+="__________________________________________________\n";
        return s;
    }
    */
}
