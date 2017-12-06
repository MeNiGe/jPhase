package jphase;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;


/**
 * This class allows the creation and manipulation 
 * of probabilistic Correlated HyperErlang (pCHE) distributions. The associated matrix has  
 * dense representation 
 * @author Juan F. Perez
 * @version 1.0
 */
public class CorrProbHyperErlangVar  extends AbstractContPhaseVar implements PhaseVar {
    
    /**
     * Total Number of phases in the HyperErlang distribution 
     */
    private int N;
    
    /**
     * Number of branches 
     */
    private int M;
    
    /**
     * Number of phases on each Erlang branch
     */
    private int[] r;
    
    /**
     * Probability mass function of choosing each one of the branches (first replica) 
     */
    private double[] alphas;
    
    /**
     * Probability mass function of choosing each one of the branches (from second to last replica) 
     */
    private double[][] betas;
    
    /**
     * Rates of each Erlang path
     */
    private double[] lambdas;
    
    
    /**
     * Constructor of a Correlated Hyper Erlang variable in dense representation.
     * As default it has just one branch that is taken with probability
     * one. the unique branch has one phase with rate 1 per time unit. 
     */
    public CorrProbHyperErlangVar() {
        this.N=1;
        this.M=1;
        this.r=new int[this.M];
        this.r[0]=1;
        this.alphas=new double[this.M];
        this.alphas[0]=1;
        this.betas=new double[this.M][this.M];
        this.betas[0][0]=1;
        this.lambdas=new double[this.M];
        this.lambdas[0]=1;
    }
    
    /**
     * Constructor of a Hyper Erlang variable with n phases
     * in dense representation
     * @param n Total number of phases
     */
    public CorrProbHyperErlangVar(int n){
        this.N=n;
    }
    
    /**
     * Constructor of a Hyper Erlang variable in dense representation
     * @param N Total number of phases
     * @param M Number of branches
     * @param r Number of phases in each path
     * @param alphas Probability associated to each branch (first replica)
     * @param betas Probability associated to each path (replicas 2 to last)
     * @param lambdas Rate associated to each branch
     * @param deep True if this is a deep constructor, false if not
     */
    public CorrProbHyperErlangVar(int N, int M, int[] r, double[] alphas, double[][] betas, double[] lambdas, boolean deep){
        if(M==alphas.length && r.length==lambdas.length && r.length==betas.length){
            this.N=N;
            this.M=M;
            if(deep==true){
            	this.alphas=new double[M];
            	System.arraycopy(alphas,0,this.alphas,0,M);
            	
                this.r=new int[this.M];
                this.betas=new double[this.M][this.M];
                this.lambdas=new double[this.M];
                System.arraycopy(r,0,this.r,0,this.M);
                System.arraycopy(lambdas,0,this.lambdas,0,this.M);
                System.arraycopy(betas,0,this.betas,0,this.M);
            }else{
                this.r=r;
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
    public CorrProbHyperErlangVar(int[] r, double[] alphas, double[][] betas, double[] lambdas, boolean deep) {
        if(r.length==lambdas.length && r.length==betas.length){
            this.M=alphas.length;
            
            if(deep==true){
            	this.alphas=new double[M];
            	System.arraycopy(alphas,0,this.alphas,0,M);
            	
                this.r=new int[M];
                this.betas=new double[this.M][this.M];
                this.lambdas=new double[M];
                System.arraycopy(r,0,this.r,0,this.M);
                System.arraycopy(betas,0,this.betas,0,this.M);
                System.arraycopy(lambdas,0,this.lambdas,0,this.M);
            }else{
                this.r=r;
                this.alphas=alphas;
                this.betas=betas;
                this.lambdas=lambdas;
            }
            
            // total # of phases
            for(int i=0;i<this.M;i++){
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
    public double[][] getBetas(){
        return this.betas;
    }
    
    /**
     * 
     * @param betas Probability of each path
     */
    public void setBetas(double[][] betas){
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
        int l =0;
        double[][] matrix = new double[N][N];
        for (int i = 0; i < this.getM(); i++) {
            for (int j = 0; j < this.getR()[i]; j++) {
                matrix[l][l] = -this.getLambdas()[i];
                if (l < N - 1 && j < this.getR()[i]-1) {
                    matrix[l][l + 1] = this.getLambdas()[i];
                }
                l++;
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
	
	public double[][] getBetaArray() {
        return this.getBetas();
    }

	
    
    /**
     * Returns the Double MAtrix that represents the variable -> first replica
     * @return Double Vector
     */
    public double[] getDVector(){
        int N = this.getN();
        int[] h = this.getR();
        
        double[] vector = new double[N];
//        System.arraycopy(this.getAlphas(),0,vector,0,this.getM());
        int l =0;
        for (int i = 0; i < this.getM(); i++) {
            for (int j = 0; j < this.getR()[i]; j++) {
            	if(j == 0){
                    vector[l] = this.getAlphas()[i];	
            	}
            	else 
            		vector[l] = 0;
                l++;
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

	
	public ContPhaseVar copy() {
        return new CorrProbHyperErlangVar(
                this.N, this.M, this.r,
                this.alphas, this.betas, this.lambdas, true
        );
    }
    
	public ContPhaseVar newVar(int n){
        return new CorrProbHyperErlangVar(n);
    }
    
    @Override
    public int getNumPhases(){
        return this.getN();
    }
   
        
    @Override
    public double expectedValue() {
        return moment(1); 
    }

    /*
    @Override -- definition of moments not clear as it depends on the replica (first or not)
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
    }*/
    
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
