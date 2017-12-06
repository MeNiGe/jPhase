package jphase.GUI;

import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;

import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextArea;
import javax.swing.JTextField;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.*;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import jmarkov.jqbd.solvers.CLPAlgorithm;
import jmarkov.jqbd.solvers.ModBoundSolverFromMatrix;
import jmarkov.jqbd.solvers.MtjLogRedSolverFromMatrix;
import jmarkov.jqbd.solvers.QBDPhaseSolver;
import jphase.ContPhaseVar;
import jphase.DenseContPhaseVar;
import jphase.DiscPhaseVar;
import jphase.MatrixUtils;
import jphase.PhaseVar;


public class eqResidualTimec extends JDialog implements  ItemListener, ActionListener{
	public final String ACEPT = "Get Residual Time";	
	public final String CANCEL = "Cancel";
	public JButton aceptB;	
	public JButton cancelB;
	public ArrayList<PhaseVarInfo> variables;
	public JComboBox arrivals;
	public JComboBox arrivals2;
	public ContPhaseVar l; 
	//public ContPhaseVar n; 
	public ContPhaseVar s;
	String var;
	String var2;
	//public PhaseVarInfo var;
	//public PhaseVarInfo var2;
	ContPhaseVar res;
	private JTextArea log;
	private TreeManagerPanel treeContentPane;
	public PhaseVar var3 = null;
	public JTextField e;
	public String varName1;
	//public ButtonGroup algorith;
	//public JRadioButton Ualg;
	//public JRadioButton LR;
	//public JRadioButton ModBon;

	MainFrame window;
	public eqResidualTimec(MainFrame window){
		this.window = window;
		setTitle("JPhase - ReisudalTime Generator");

		/**
		 * algorith=new ButtonGroup();
		Ualg=new JRadioButton("Linear Progresion");
		algorith.add(Ualg);
		LR=new JRadioButton("Logarithmic Reduction");
		algorith.add(LR);
		ModBon=new JRadioButton("Modified Boundary");
		algorith.add(ModBon);
		LR.setSelected(true);
		
		 */

		aceptB = new JButton(ACEPT);
		aceptB.setActionCommand(ACEPT);
		aceptB.addActionListener( this );

		cancelB = new JButton(CANCEL);
		cancelB.setActionCommand(CANCEL);
		cancelB.addActionListener( this );
		setLayout(new BorderLayout());

		JPanel center = new JPanel();
		JPanel ncenter = new JPanel();
	
		/**
		 * JPanel algo = new JPanel();
		 */
		
		JPanel buttons = new JPanel();
		JPanel left = new JPanel();
		JPanel right = new JPanel();

		variables = window.getVariables();

		e= new JTextField("                     ");
		add(e);
		
		arrivals = new JComboBox( getNames(variables).toArray() );
		arrivals.setSelectedIndex( 0 );
		arrivals.addItemListener( this );
		add( arrivals );

		//arrivals2 = new JComboBox( getNames(variables).toArray() );
		//arrivals2.setSelectedIndex( 0 );
		//arrivals2.addItemListener( this );
		//add( arrivals2 );

		/**
		 * services = new JComboBox( getNames(variables).toArray() );
		services.setSelectedIndex( 0 );
		services.addItemListener( this );

		 */
		
		JLabel label1 = new JLabel("Please, select a Case");
		JLabel label2 = new JLabel("In order to get Min Values");
		//JLabel label3 = new JLabel("Min values");
		label1.setHorizontalAlignment(javax.swing.SwingConstants.CENTER); 
		label2.setHorizontalAlignment(javax.swing.SwingConstants.CENTER); 

		
		JLabel varname = new JLabel("New Var Name");
		JLabel arrival = new JLabel("Case");
		//JLabel service = new JLabel("Case2");      
		 

		//left.setLayout(new GridLayout(2,1));
		left.add( varname );
		left.add( e );
		right.add( arrival);
		right.add( arrivals );

		//right.setLayout(new GridLayout(2,1));
		
	
		 //right.add( arrival );
		 //right.add( arrivals2 );
		 
		   

		center.setLayout(new BorderLayout());
		center.add( left , BorderLayout.CENTER );
	
		
		center.add( right , BorderLayout.SOUTH);
		 

	/**
	 * 	algo.setLayout(new GridLayout(3,1));
		algo.add( LR );
		algo.add( ModBon );
		algo.add( Ualg );
	 ncenter.add( algo , BorderLayout.EAST );
	 */

		ncenter.setLayout(new BorderLayout());
		ncenter.add( center , BorderLayout.CENTER );
		

		JPanel labels = new JPanel();
		labels.setLayout(new GridLayout(3,1));
		labels.add(label1);
		labels.add(label2);
		

		add( labels , BorderLayout.NORTH );
		add( ncenter , BorderLayout.CENTER );

		buttons.add( cancelB ); 
		buttons.add( aceptB );

		add( buttons , BorderLayout.SOUTH );

		setSize(700, 600);
		setResizable(true);
		setLocationRelativeTo
		(null);
	}

	public void itemStateChanged(ItemEvent arg0) {		
	}

	public void actionPerformed(ActionEvent arg0) {
		String message = arg0.getActionCommand();
		
		 	//variables = new ArrayList<PhaseVarInfo>();
		 	//category1 = addObject(null, "Set 1");
		 	//var = new PhaseVarInfo("Expo 5", 
	        		//DenseContPhaseVar.expo(5));
	        //addObject(category1, var);
		 	if(message.equals(ACEPT)){
		 		this.varName1 = e.getText();
		 		var = (String) arrivals.getSelectedItem();
		 		var2 = (String) arrivals2.getSelectedItem();
		 		int arri = arrivals.getSelectedIndex();
				int arri2 = arrivals2.getSelectedIndex();
		 		l = (ContPhaseVar)variables.get(arri).var;
		 		//n = (ContPhaseVar)variables.get(arri2).var;
				this.res = l.eqResidualTime();
				
				window.getTreeConstentPane().addVar(varName1, res);
				window.pack();
				
			
			/**
			 * long ctime = System.currentTimeMillis();
			int arri = arrivals.getSelectedIndex();
			if(Ualg.isSelected()){
				test = new CLPAlgorithm(variables.get(arri).var, variables.get(serv).var);
			}
			else if(ModBon.isSelected()){
				test = new ModBoundSolverFromMatrix(variables.get(arri).var, variables.get(serv).var);
			}
			else{
				test = new MtjLogRedSolverFromMatrix(variables.get(arri).var, variables.get(serv).var);
			}
			if(test.unstableSystem()){
				JOptionPane.showMessageDialog( this, "Unstable system", "Performance Measures", JOptionPane.ERROR_MESSAGE );
			}
			else{
				double[][] R = test.getRmatrix();
				test.printMatrices();
				String resp = test.performanceMeasures(R);
				ctime = System.currentTimeMillis() - ctime;
				resp += "\nTime duration: " + ctime + " milisecond"; 
				resp += (ctime!=1)?"s":"";
				JOptionPane.showMessageDialog( this, resp, "Performance Measures", JOptionPane.INFORMATION_MESSAGE );


				InputFrame probFrame = new InputFrame(
						"Steady State",
				"Number of states");
				probFrame.setVisible(true);
				probFrame.setFocusable(true);

				if(probFrame.res){
					try{
						int n = Integer.parseInt(probFrame.getValue());						
						ArrayList<Double> pis = test.getLevelSteadyStateProbs();				
						ArrayList<DenseMatrix> pisP = test.getSteadyStateProbsPerLevel();

						String ans = "";

						for(int i = 0 ;i < n && i < pis.size(); i++ ){
							ans += "Pi(" + i + ")=" + String.format("%6.4f", pis.get(i)) + "\t{";
							DenseMatrix doub = pisP.get(i);
							double[][] matrix = Matrices.getArray(doub);
							for(int j = 0 ;j < doub.numColumns(); j++ ){
								ans += String.format("%6.4f", matrix[0][j]) + "|";
							}
							ans += "}\n";
						}
						JOptionPane.showMessageDialog( this, ans, "Steady State Probabilities", JOptionPane.INFORMATION_MESSAGE );
					}
					catch(NumberFormatException e){

					}
				}
			 */
			
			
		
	}else if(message.equals(CANCEL)){
			dispose();
		}
		 this.setVisible(false);
	}

	public ArrayList<String> getNames(ArrayList<PhaseVarInfo> vars){
		ArrayList<String> names = new ArrayList<String>();
		for(PhaseVarInfo var : vars)
			names.add(var.varName);
		return names;
	}

	public PhaseVarInfo getVar(String name){
		for(PhaseVarInfo var : variables){
			if(var.varName.equals(name))
				return var;			
		}
		return null;
	}
	
}
