package com.alessandrognoatto.matrixDecompositions;

import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.complex.*;

public class SchurDecomposition {
	FieldMatrix<Complex> U;
	FieldMatrix<Complex> T;
	
	/**
	 * This class computes the Schur decomposition of a complex matrix A by iterating the QR decomposition.
	 * For matrices whose eigenvalues are "very close" among each other this may not be reliable.
	 * 
	 * See http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter3.pdf
	 * 
	 * @param A
	 */
	public SchurDecomposition(FieldMatrix<Complex> A){
		int maxIter = 10;
		int n = A.getColumnDimension();
		
		FieldMatrix<Complex> Ak = A;
		FieldMatrix<Complex> Uk = MatrixUtils.createFieldIdentityMatrix(ComplexField.getInstance(), n);
		
		for(int i = 1; i <= maxIter; i++){
			HausHolderQR myHaus = new HausHolderQR(Ak);
			Ak = (myHaus.getR()).multiply(myHaus.getQ());
			Uk = Uk.multiply(myHaus.getQ());
					
		}
		
		this.T = Ak;
		this.U = Uk;
	}
	
	public FieldMatrix<Complex> getU(){
		return this.U;
	}
	
	public FieldMatrix<Complex> getT(){
		return this.T;
	}
	
	public static void main(String[] args){
		
		Complex A11 = new Complex(0.4854,0.9157);
		Complex A12 = new Complex(0.1419,0.9595);
		Complex A21 = new Complex(0.8003,0.7922);
		Complex A22 = new Complex(0.4218,0.6557);
		
		Complex[][] data2 = {{A11,A12},{A21,A22}};
		FieldMatrix<Complex> myMatrix2 = MatrixUtils.createFieldMatrix(data2);
		
		SchurDecomposition mySchur = new SchurDecomposition(myMatrix2);
		System.out.println("Schur decomposition");
		System.out.println("The matrix T is");
		System.out.println(mySchur.getT());
		System.out.println("Matlab result for T");
		System.out.println("0.925766799065930 + 1.726648232444361i  0.033664159886327 - 0.163938615885574i");
		System.out.println("0.000000000000000 + 0.000000000000000i -0.018566799065931 - 0.155248232444361i");
		
		System.out.println("The matrix U is");
		System.out.println(mySchur.getU());
		System.out.println("Matlab result for U");
		System.out.println("0.724277100667501 + 0.017630898571263i -0.652911034588147 - 0.220950251815401i");
		System.out.println("0.652911034588147 - 0.220950251815396i  0.724277100667501 - 0.017630898571256i");
	}
		
	


}
