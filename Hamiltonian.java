//Author: Reza Niazi 
//Last Updated: 8/3/15
//This is a program that that constructs the Hamiltonian for a 2 shell
//system for a variable amount of nucleons. 
//It then diagonalizes this Hamiltonian, and finds its
//eigenvalues and eigenvectors, representd with 
//as matrices. The matrix elements, eigenvalues,
//and eigenvectors are seperately printed in output
//files with their appropriate names. The files' 
//first two columns represent which term in the matri
//that the printed value corresponds to. The Hamiltonian is given by
//H = epsilon*N(the number operator)-Gsum(of a,b from 1 to 2)
//PCPD+NCND, where PC is the proton pair creation operator, PD
//is the proton pair destruction operator, and NC and ND are the 
//neutron analogs. 
import java.util.Arrays; 
import java.io.*;
import Jama.*;

public class Hamiltonian {

	public static void main(String[] args) throws IOException
	
	{
		int Z1 = 2; //Must always initialize however many protons and neutrons here
		int N1 = 4; 
		BufferedWriter HMatrixElements = null;
		BufferedWriter Eigenvalues = null;
		BufferedWriter Eigenvectors = null;
		HMatrixElements = new BufferedWriter(new FileWriter("Hamiltonian.txt"));
		Eigenvalues = new BufferedWriter(new FileWriter("Eigenvalues.txt"));
		Eigenvectors = new BufferedWriter(new FileWriter("Eigenvectors.txt"));
		final int omega1 = 2;
		final int omega2 = 1;
		int Z2 = 0; //will always start the states in the lower shell
		int N2 = 0; //and then through iteration add them to the upper shell to 
		final int iterations = (((Z1+2)*(N1+2))/(4));
		double[][] Hamiltonian = new double[iterations][iterations];
		final int neutronIterations = (N1/2) + 1;
		final int protonIterations = (Z1/2) + 1;
		double eigenValue1 = 0.0; //for nn11
		double eigenValue2 = 0.0; //for nn21 
		double eigenValue3 = 0.0; //for nn12
		double eigenValue4 = 0.0; //for nn 22
		double eigenValue5 = 0.0; //for pp11
		double eigenValue6 = 0.0; //for pp21
		double eigenValue7 = 0.0; //for pp12
		double eigenValue8 = 0.0; //for pp22
		double eigenValue9 = 0.0; //for number operator 
		double temp = 0.0; 
		int[] bra1 = ket(Z1,N1,omega1);
		int[] bra2 = ket(Z2,N2,omega2);
		int[] ket1_11_N = ket(Z1,N1,omega1); //must create a new ket state for each SdaggerS 
		int[] ket2_11_N = ket(Z2,N2,omega2); //since the arrays are objects and they will
		int[] ket1_12_N = ket(Z1,N1,omega1); //store what's previously done, which would throw
		int[] ket2_12_N = ket(Z2,N2,omega2); //off the answer 
		int[] ket1_21_N = ket(Z1,N1,omega1); //_(ij)_ is for which operator pair its acting on
		int[] ket2_21_N = ket(Z2,N2,omega2); //_(N or P) is for whether its for the proton or neutron operator
		int[] ket1_22_N = ket(Z1,N1,omega1);
		int[] ket2_22_N = ket(Z2,N2,omega2);
		int[] ket1_11_P = ket(Z1,N1,omega1);
		int[] ket2_11_P = ket(Z2,N2,omega2);
		int[] ket1_12_P = ket(Z1,N1,omega1);
		int[] ket2_12_P = ket(Z2,N2,omega2);
		int[] ket1_21_P = ket(Z1,N1,omega1);
		int[] ket2_21_P = ket(Z2,N2,omega2);
		int[] ket1_22_P = ket(Z1,N1,omega1);
		int[] ket2_22_P = ket(Z2,N2,omega2);
		int[] ket1_NO = ket(Z1,N1,omega1); //NO for number operator
		int[] ket2_NO = ket(Z2,N2,omega2);
		int nBraIncrement = 0; //used to increase the number of neutrons in second level
		int nBraIndex = 0; //used to increase the number of nucleons by a factor of 2
		int zBraIncrement = 0;
		int zBraIndex = 0;
		int nKetIncrement = 0;
		int nKetIndex = 0;
		int zKetIndex = 0 ;
		int zKetIncrement =0;
		double G = 1.0;
		double epsilon = 20.0;
		int row = -1;
		int columnCounter = -1; //need columnCounter to iterate and
		int column; //need column to get modulated and give the right index
		
		
		{
			for(int i=N1;i>=0;i=i-2) //2 for loops will create all possible bra states 
			{
				
				bra1[1] = i; //assigns the N1 as the N1 in the loop, so it goes through all possible values
				nBraIncrement = 2*nBraIndex; //this allows one to decrement N1 by 2 and inrement N2 by 2 
				bra2[1] = nBraIncrement; //so all configurations of N within the 2 shells will be achieved
				nBraIndex = nBraIndex + 1;
				zBraIndex =0; //this makes sure that zIndex returns to 0 at the end of the inner loop, so no over-counting
				
				for(int j=Z1;j>=0;j=j-2)
				{
					row++;
					bra1[0] = j;
					zBraIncrement = 2*zBraIndex;
					bra2[0] = zBraIncrement;
					zBraIndex = zBraIndex + 1;
					nKetIndex = 0;
					
					
					for(int k=N1;k>=0;k=k-2) //2 for loops will create all possible ket states 
					{
						ket1_11_N[1] = k; 
						ket1_12_N[1] = k; 
						ket1_21_N[1] = k; 
						ket1_22_N[1] = k;
						ket1_11_P[1] = k; 
						ket1_12_P[1] = k; 
						ket1_21_P[1] = k; 
						ket1_22_P[1] = k;
						ket1_NO[1] = k;
						nKetIncrement = 2*nKetIndex;  
						ket2_11_N[1] = nKetIncrement; 
						ket2_12_N[1] = nKetIncrement; 
						ket2_21_N[1] = nKetIncrement; 
						ket2_22_N[1] = nKetIncrement;
						ket2_11_P[1] = nKetIncrement; 
						ket2_12_P[1] = nKetIncrement; 
						ket2_21_P[1] = nKetIncrement; 
						ket2_22_P[1] = nKetIncrement;
						ket2_NO[1] = nKetIncrement;
						nKetIndex = nKetIndex + 1;
						zKetIndex =0; 
						temp = 0.0;
						
						for(int l=Z1;l>=0;l=l-2)
						{
							columnCounter++;
							column = columnCounter%iterations;	
							Hamiltonian[row][column] = 0.0;
							ket1_11_N[0] = l;
							ket1_12_N[0] = l;
							ket1_21_N[0] = l;
							ket1_22_N[0] = l;
							ket1_11_P[0] = l;
							ket1_12_P[0] = l;
							ket1_21_P[0] = l;
							ket1_22_P[0] = l;
							ket1_NO[0] = l;
							zKetIncrement = 2*zKetIndex;
							ket2_11_N[0] = zKetIncrement;
							ket2_12_N[0] = zKetIncrement;
							ket2_21_N[0] = zKetIncrement;
							ket2_22_N[0] = zKetIncrement;
							ket2_11_P[0] = zKetIncrement;
							ket2_12_P[0] = zKetIncrement;
							ket2_21_P[0] = zKetIncrement;
							ket2_22_P[0] = zKetIncrement;
							ket2_NO[0] = zKetIncrement;
							zKetIndex = zKetIndex + 1;
							//System.out.println("<" + " " + bra1[0] + " " + bra1[1] + " " + bra2[0] + " " + bra2[1] + "|" + " " + "|" + ket1_11_N[0] + " " + ket1_11_N[1] + " " + ket2_11_N[0] + " " + ket2_11_N[1] + ">");
							temp = temp + neutronHamiltonian(ket1_11_N,ket2_11_N,G,1,1);
							if((bra1[0] == ket1_11_N[0]) & (bra1[1] == ket1_11_N[1]) & (bra2[0] == ket2_11_N[0]) & (bra2[1] == ket2_11_N[1]))
							{
								eigenValue1 = eigenValue1 + temp;
							}
							Hamiltonian[row][column] = Hamiltonian[row][column] + eigenValue1;
							eigenValue1 = 0.0;
							temp = 0.0;
						
								
							temp = temp + neutronHamiltonian(ket1_12_N,ket2_12_N,G,1,2);
							if((bra1[0] == ket1_12_N[0]) & (bra1[1] == ket1_12_N[1]) & (bra2[0] == ket2_12_N[0]) & (bra2[1] == ket2_12_N[1]))
							{
								eigenValue2 = eigenValue2 + temp;
							}
							Hamiltonian[row][column] = Hamiltonian[row][column] + eigenValue2;
							eigenValue2 = 0.0;
							temp = 0.0;
							
							temp = temp + neutronHamiltonian(ket1_21_N,ket2_21_N,G,2,1);
								if((bra1[0] == ket1_21_N[0]) & (bra1[1] == ket1_21_N[1]) & (bra2[0] == ket2_21_N[0]) & (bra2[1] == ket2_21_N[1]))
								{
									eigenValue3 = eigenValue3 + temp;
								}
							
						    Hamiltonian[row][column] = Hamiltonian[row][column] + eigenValue3;
							eigenValue3 = 0.0;
							temp = 0.0;
							
							
							temp = temp + neutronHamiltonian(ket1_22_N,ket2_22_N,G,2,2);
							if((bra1[0] == ket1_22_N[0]) & (bra1[1] == ket1_22_N[1]) & (bra2[0] == ket2_22_N[0]) & (bra2[1] == ket2_22_N[1]))
							{
								eigenValue4 = eigenValue4 + temp;
							}
							Hamiltonian[row][column] = Hamiltonian[row][column] + eigenValue4;
							eigenValue4 = 0.0;
							temp = 0.0;
							
							
							temp = temp + protonHamiltonian(ket1_11_P,ket2_11_P,G,1,1);
							if((bra1[0] == ket1_11_P[0]) & (bra1[1] == ket1_11_P[1]) & (bra2[0] == ket2_11_P[0]) & (bra2[1] == ket2_11_P[1]))
							{
								eigenValue5 = eigenValue5 + temp;
							}
							Hamiltonian[row][column] = Hamiltonian[row][column] + eigenValue5;
							eigenValue5 = 0.0;
							temp = 0.0;
							
							temp = temp + protonHamiltonian(ket1_12_P,ket2_12_P,G,1,2);
							if((bra1[0] == ket1_12_P[0]) & (bra1[1] == ket1_12_P[1]) & (bra2[0] == ket2_12_P[0]) & (bra2[1] == ket2_12_P[1]))
							{
								eigenValue6 = eigenValue6 + temp;
							}
							Hamiltonian[row][column] = Hamiltonian[row][column] + eigenValue6;
							eigenValue6 = 0.0;
							temp = 0.0;
							
							
							temp = temp + protonHamiltonian(ket1_21_P,ket2_21_P,G,2,1);
							if((bra1[0] == ket1_21_P[0]) & (bra1[1] == ket1_21_P[1]) & (bra2[0] == ket2_21_P[0]) & (bra2[1] == ket2_21_P[1]))
							{
								eigenValue7 = eigenValue7 + temp;
							}
							Hamiltonian[row][column] = Hamiltonian[row][column] + eigenValue7;
							eigenValue7 = 0.0;
							temp = 0.0;
						
							temp = temp + protonHamiltonian(ket1_22_P,ket2_22_P,G,2,2);
							if((bra1[0] == ket1_22_P[0]) & (bra1[1] == ket1_22_P[1]) & (bra2[0] == ket2_22_P[0]) & (bra2[1] == ket2_22_P[1]))
							{
								eigenValue8 = eigenValue8 + temp;
							}
							Hamiltonian[row][column] = Hamiltonian[row][column] + eigenValue8;
							eigenValue8 = 0.0;
							temp = 0.0;
							
							temp = temp + energyGap(ket2_NO,epsilon);
							if((bra1[0] == ket1_NO[0]) & (bra1[1] == ket1_NO[1]) & (bra2[0] == ket2_NO[0]) & (bra2[1] == ket2_NO[1]))
							{
								eigenValue9 = eigenValue9 + temp;
							}
							Hamiltonian[row][column] = Hamiltonian[row][column] + eigenValue9;
							eigenValue9 = 0.0;
							temp = 0.0;
							HMatrixElements.write(row + " " + column + " " + Double.toString(Hamiltonian[row][column]));
							HMatrixElements.newLine();
						}
					
					
						}
					
					}
					
				}
			
		
		Matrix H = new Matrix(Hamiltonian);
		EigenvalueDecomposition H_Diag = new EigenvalueDecomposition(H);
		for(int m=0;m<=iterations-1;m++)
			{
			for(int n=0;n<=iterations-1;n++)
				{
					Eigenvalues.write(Double.toString(H_Diag.getD().get(m,n)));
					Eigenvalues.newLine();
					Eigenvectors.write(Double.toString(H_Diag.getV().get(m,n)));
					Eigenvectors.newLine();
				}
			
			}
		
		HMatrixElements.flush();  
		HMatrixElements.close();
		Eigenvalues.flush();
		Eigenvalues.close();
		Eigenvectors.flush();
		Eigenvalues.close();
		}
	}

	
	public static int[] ket(int Z1,int N1,int omega1)
	//Z is the number of protons number, N is the number of neutrons
	//omega is the omega value of that shell
	{
		int[] state1 = {Z1,N1,omega1};
		return state1;
	}
		
	public static int numberOperator(int[] ket2)
	//tells you the number of nucleons, but only for the second shell, so the ket2 is there as a reminder
	{
		int nucleonNumber = ket2[0] + ket2[1];
		return nucleonNumber;
	}
	
	
	public static double pairCreation(int[] ket1, int[] ket2, int shell,String nucleon)
	//if shell=1, then acting on ket1/shell 1
	//if shell=2, then acting on ket2/ shell 2
	//if nucleon = p, then acting on protons
	//if nucleon = n, then acting on neutrons
	
	{
		double eigenValue;
		
		if(shell == 1)
		{
			if(nucleon.equals("p"))
			{
				eigenValue = Math.sqrt((ket1[0]+2)*(2*ket1[2]-ket1[0]))/(2.0); //this is by definition
				ket1[0] = ket1[0] + 2; //the creation operator creates 2 more particles in a certain shell
			}
			
			else
			{
				eigenValue = Math.sqrt((ket1[1]+2)*(2*ket1[2]-ket1[1]))/(2.0);
				ket1[1] = ket1[1] + 2;
			}
		}
		
		else
		{
			if(nucleon.equals("p"))
			{
				eigenValue = Math.sqrt((ket2[0]+2)*(2*ket2[2]-ket2[0]))/(2.0);
				ket2[0] = ket2[0] + 2;
			}
			
			else
			{
				eigenValue = Math.sqrt((ket2[1]+2)*(2*ket2[2]-ket2[1]))/(2.0);
				ket2[1] = ket2[1] + 2;
			}	
		}
		
		return eigenValue;
		
	}
	
	public static double pairDestruction(int[] ket1, int[] ket2, int shell,String nucleon)
	
	{
		double eigenValue;
		
		if(shell == 1)
		{
			if(nucleon.equals("p"))
			{
				eigenValue = Math.sqrt(ket1[0]*(2*ket1[2]+2-ket1[0]))/(2.0); //these are by definition
				ket1[0] = ket1[0] - 2;	//the destruction operator destroys 2 particles in some shell
			}
			
			else
				
			{				
				eigenValue = Math.sqrt(ket1[1]*(2*ket1[2]+2-ket1[1]))/(2.0);
				ket1[1] = ket1[1] - 2;
				
			}
			
		}
	
		else
		{
			if(nucleon.equals("p"))
			{
				eigenValue = Math.sqrt(ket2[0]*(2*ket2[2]+2-ket2[0]))/(2.0);
				ket2[0] = ket2[0] - 2;
			}
			
			else
			{
				
				eigenValue = Math.sqrt(ket2[1]*(2*ket2[2]+2-ket2[1]))/(2.0);
				ket2[1] = ket2[1] - 2;
			}
			
		}	
		
		return eigenValue;
		
	}
	
	public static double energyGap(int[] ket, double epsilon)
	{
		//calculates the energy gap term in Hamiltonian
		double energyGap = epsilon*numberOperator(ket);
		return energyGap;
	}
	
	public static double protonHamiltonian(int[] ket1, int[] ket2, double G, int shell_1, int shell_2)
	{
		//calculates the proton contribution from the Hamiltonian based on the shell
		double pTerm = -G*((pairDestruction(ket1,ket2,shell_1,"p"))*(pairCreation(ket1,ket2,shell_2,"p")));
		return pTerm;
	}
	
	public static double neutronHamiltonian(int[] ket1, int[] ket2, double G, int shell_1,int shell_2)
	{
		//calculates the neutron contribution from the Hamiltonian based on the shell
		double nTerm = -G*((pairDestruction(ket1,ket2,shell_1,"n"))*(pairCreation(ket1,ket2,shell_2,"n")));
		return nTerm;
	}
	
	
	}
	
