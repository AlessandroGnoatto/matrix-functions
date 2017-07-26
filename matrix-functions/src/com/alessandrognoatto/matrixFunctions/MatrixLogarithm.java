package com.alessandrognoatto.matrixFunctions;

import java.util.Arrays;

import org.apache.commons.math3.complex.*;
import org.apache.commons.math3.linear.FieldLUDecomposition;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.linear.MatrixUtils;

import com.alessandrognoatto.basics.FindInArray;
import com.alessandrognoatto.basics.GaussLegendre;
import com.alessandrognoatto.matrixDecompositions.SchurDecomposition;

/**
 * @author Alessandro Gnoatto
 */
public class MatrixLogarithm {
	
	private static double[] xvals = {1.586970738772063e-005,
        								2.313807884242979e-003,
							        1.938179313533253e-002,
							        6.209171588994762e-002,
							        1.276404810806775e-001,
							        2.060962623452836e-001,
							        2.879093714241194e-001};
	private static int mmax = 7;
	
	private static int maxsqrt = 100;
	
	
	/**
	 * Computation of the logarithm of a Complex Matrix. 
	 * We assume that the matrix A has no nonpositive real eigenvalues.
	 * The computation is performed via inverse scaling and squaring method 
	 * with Pade approximation and a Schur decomposition.
	 * 
	 * Reference: A. H. Al-Mohy and N. J. Higham, Improved Inverse Scaling 
     * and Squaring Algorithms for the Matrix Logarithm, MIMS EPrint 2011.83,
     * The University of Manchester, October 2011.  
     * Name of corresponding algorithm in that paper: 
     * Algorithm 4.1/iss_schur_new.
     * 
     * @author Alessandro Gnoatto based on the Matlab code by Awad H. Al-Mohy and Nicholas J. Higham
	 * @param A
	 * @return the matrix logarithm of A 
	 */
	public static FieldMatrix<Complex> compute(FieldMatrix<Complex> A){
		int n = A.getColumnDimension();
		boolean foundm = false;
		
		SchurDecomposition mySchur = new SchurDecomposition(A);
		
		FieldMatrix<Complex> T = mySchur.getT();
		FieldMatrix<Complex> T0 = mySchur.getT();
		
		int p = 0;
		int s0 = optCost(T);
		int s = s0;
		
		for(int k = 1; k<= Math.min(s, maxsqrt);k++){
			T = sqrtmTri(T);
		}
		
		FieldMatrix<Complex> eye = MatrixUtils.createFieldIdentityMatrix(ComplexField.getInstance(), n);
		double d2 = pNorm(T.subtract(eye),2);
		double d3 = pNorm(T.subtract(eye),3);
		double alpha2 = Math.max(d2,	d3);
		
		int m = 0;
		
		if(alpha2 <= xvals[1]){
			m = FindInArray.findFirstGreaterThan(subArray(xvals,0,1), alpha2);
			foundm = true;
		}
		
		
		boolean more;
		
		while (!foundm){
			more = false;
			if(s > s0){
				d3 = pNorm(T.subtract(eye),3);
			}
			double d4 = pNorm(T.subtract(eye),4);
			double alpha3 = Math.max(d3,d4);
			
			if(alpha3 <= xvals[mmax-1]){
				int j = FindInArray.findFirstGreaterThan(subArray(xvals,2,mmax-1), alpha3)+2;
				
				if(j <= 5){
					m = j;
					break;
				} else{
					if(alpha3/2<= xvals[4] && p < 2){
						more = true;
						p = p + 1;
					}
				}
			}
			
			if(!more){
				double d5 = pNorm(T.subtract(eye),5);
				double alpha4 = Math.max(d4, d5);
				double eta = Math.min(alpha3, alpha4);
				
				if(eta<= xvals[mmax-1]){
					m = FindInArray.findFirstGreaterThan(subArray(xvals,5,mmax-1), eta)+5;
					break;
				}
			}
			
			if(s == maxsqrt){
				m = mmax;
				break;
			}
			
			T = sqrtmTri(T);
			s = s + 1;
		}
		
		//Computation of the superdiagonal of T^(1/2^s)
		for(int k = 0; k < n-1;k++){
			int[] myIndex = {k,k+1};
			FieldMatrix<Complex> utils = powerm2by2(T0.getSubMatrix(myIndex, myIndex),1/2^s);
			T.setEntry(k, k, utils.getEntry(k, k));
			T.setEntry(k, k+1, utils.getEntry(k, k+1));
			T.setEntry(k+1, k, utils.getEntry(k+1, k));
			T.setEntry(k+1, k+1, utils.getEntry(k+1, k+1));
		}
		
		//Computation of the diagonal of T^(1/2^s)-I
		//sqrtPower works on complex scalars so
		
		Complex[] d = new Complex[T0.getColumnDimension()];
		
		for(int k = 0; k < d.length-1; k++){
			d[k] = sqrtPower(T0.getEntry(k, k),s);
			T.setEntry(k, k, d[k]);
		}
		
		FieldMatrix<Complex> Y = logmPf(T,m);
		Complex factor = new Complex(2^s,0);
		FieldMatrix<Complex> X = Y.scalarMultiply(factor);
		
		//Compute diagonal and superdiagonal of log(T)
		for(int k = 0; k < n-1;k++){
			int[] myIndex = {k,k+1};
			FieldMatrix<Complex> utils = logm2by2(T0.getSubMatrix(myIndex, myIndex));
			X.setEntry(k, k, utils.getEntry(k, k));
			X.setEntry(k, k+1, utils.getEntry(k, k+1));
			X.setEntry(k+1, k, utils.getEntry(k+1, k));
			X.setEntry(k+1, k+1, utils.getEntry(k+1, k+1));
		}
		
		FieldMatrix<Complex> Uherm = (mySchur.getU()).transpose();
		for(int k = 0; k<n; k++){
			for(int l = 0; l<n;l++){
				Uherm.setEntry(k, l, Uherm.getEntry(k,l).conjugate());
			}
		}
		
		return ((mySchur.getU()).multiply(X)).multiply(Uherm); 
	}
	/**
	 * Extracts elements between firstIndex and lastIndex and returns the corresponding
	 * subarray
	 * 
	 * @param myArray
	 * @param firstIndex
	 * @param lastIndex
	 * @return
	 */
	private static double[] subArray(double[] myArray, int firstIndex, int lastIndex){
		
		int n = lastIndex - firstIndex +1;
		double[] myNewArray = new double[n];
		
		for(int i = firstIndex; i <= lastIndex; i++){
			myNewArray[i-firstIndex] = myArray[i];
		}
		return myNewArray;
	}
	
	/**
	 * Compute sum(abs(v).^p)^(1/p).
	 * @param myArray
	 * @param p
	 * @return
	 */
	private static double pNorm(FieldMatrix<Complex> A, int p){
		
		double mySum = 0;
		
		for(int i = 0; i<A.getRowDimension(); i++){
			for(int j = 0; j<A.getColumnDimension(); j++){
				mySum += Math.pow((A.getEntry(i, j)).abs(),p);
			}
			
		}
		
		return Math.pow(mySum, 1.0/p);
	}
	
	
	private static int optCost(FieldMatrix<Complex> A){
		
		int s = 0;
		
		Complex[] aDiag = new Complex[A.getColumnDimension()];
		
		for(int i = 0; i < aDiag.length;i++){
			aDiag[i] = A.getEntry(i, i);
		}
		
		while(infinityNormOfAMinusOne(aDiag)>xvals[mmax-1]){
			
			//update aDiag
			for(int i = 0; i < aDiag.length;i++){
				aDiag[i] = aDiag[i].sqrt();
			}
			//increment s
			s++;
		}	
		
		return s;
	}
	
	private static double infinityNormOfAMinusOne(Complex[] A){
		
		//norm(X,inf) = max(abs(X))
		double[] AbsDiagonal = new double[A.length];
		
		for(int i = 0; i<A.length; i++){
			AbsDiagonal[i] = (A[i].subtract(1.0)).abs();
		}
		//sort elements of the array
		Arrays.sort(AbsDiagonal);
		
		return AbsDiagonal[AbsDiagonal.length-1];
	}
	

	
	/**
	 * Pade approximation to matrix log by partial fraction expansion
	 * 
	 * @param A a complex Matrix
	 * @param m order of approximation
	 * @return the Matrix log
	 */
	private static FieldMatrix<Complex> logmPf(FieldMatrix<Complex> A, int m){
		
		GaussLegendre quad = new GaussLegendre(m,0.0,1.0);
		
		int n = A.getColumnDimension();
		FieldMatrix<Complex> S = MatrixUtils.createFieldMatrix(ComplexField.getInstance(),n,n);
		FieldMatrix<Complex> eye = MatrixUtils.createFieldIdentityMatrix(ComplexField.getInstance(), n);
		
		for(int j = 0; j<m;j++){
			
			//denominator
			Complex node = new Complex(quad.getNodes()[j],0);
			FieldMatrix<Complex> myMatrix = (A.scalarMultiply(node)).add(eye);
			
			//Compute the inverse of the denominator
			FieldMatrix<Complex> denominator = new FieldLUDecomposition<Complex>(myMatrix).getSolver().getInverse();
			
			Complex weight = new Complex(quad.getWeights()[j],0);
			S = S.add((A.multiply(denominator)).scalarMultiply(weight));
		}	
		
		return S;
	}

	/**
	 * Compute upper triangular square root R of T, a column at a time. The matrix T is in this context
	 * a full matrix i.e. we do not use a particular data structure for triangular matrices.
	 * 
	 * @param T
	 * @return The matrix square root or T.
	 */
	private static FieldMatrix<Complex> sqrtmTri(FieldMatrix<Complex> T){
		int n = T.getRowDimension();
		
		FieldMatrix<Complex> R = MatrixUtils.createFieldMatrix(ComplexField.getInstance(),n,n);
		
		for(int i = 0; i< R.getColumnDimension(); i++){
			R.setEntry(i, i, (T.getEntry(i, i)).sqrt());
		}
		
		for(int j = 1; j < R.getColumnDimension(); j++){
			for(int i = j-1; i >= 0; i--){
				
				Complex coeff = Complex.ZERO;
				
				for(int k = i+1; k <= j-1; k++){						
					coeff= coeff.add((R.getEntry(i, k)).multiply(R.getEntry(k,j)));
				}
				
				coeff = ((coeff.negate()).add(T.getEntry(i, j))).divide((R.getEntry(i, i)).add(R.getEntry(j,j)));
				
				R.setEntry(i,j,coeff);
			}
		}

		return R;
	}
	
	/**
	 * Computes the power of 2-by-2 upper triangular matrix.
	 * @param A (triangular matrix with complex entries)
	 * @return The p-th power of A.
	 */	
	private static FieldMatrix<Complex> powerm2by2(FieldMatrix<Complex> A, double p ){
		Complex a1 = A.getEntry(0, 0);
		Complex a2 = A.getEntry(1, 1);
		Complex a1p = a1.pow(p);
		Complex a2p = a2.pow(p);
		Complex loga1 = (a1).log();
		Complex loga2 = (a2).log();
		
		Complex X_12;
		
		if(a1 == a2){			
			X_12 = ((A.getEntry(0, 1)).multiply(p)).multiply(a1.pow(p-1));
		}else if(a1.abs()< 0.5*a2.abs() || a2.abs()< 0.5*a1.abs()){
			X_12 = (A.getEntry(0, 1)).multiply((a2p.subtract(a1p)).divide(a2.subtract(a1)));
		}else{
			Complex term1 = (((a2.subtract(a1)).divide(a2.add(a1))).atan());
			Complex term2 = (((Complex.I).multiply(Math.PI)).multiply(unwinding(loga2.subtract(loga1))));
			Complex term3 = (((loga1.add(loga2)).divide(2.0)).exp()).multiply(2.0);
			Complex w = term1.add(term2);
			X_12 = (term3.multiply((w.multiply(p)).sinh())).divide(a2.subtract(a1));			
		}
		
		Complex[][] arrayResult = {{a1p, X_12},{Complex.ZERO,a2p}};
		FieldMatrix<Complex> matrixResult = MatrixUtils.createFieldMatrix(arrayResult);
		return matrixResult;
		
	}
	
	/**
	 * Computation of a^(2^n)-1.
	 * @param Complex a
	 * @param int n
	 * @return a^(2^n)-1
	 */
	private static Complex sqrtPower(Complex a, int n){
		Complex result;
		
		if(n == 0){
			result = a.subtract(1);
		}else{
			int n0 = n;
			if(a.getArgument()>Math.PI/2){
				a = a.sqrt();
				n0 = n-1;
			}
			Complex z0 = a.subtract(1.0);
			a = a.sqrt();
			result = a.add(1.0);
			
			for(int i = 1; i <= n0-1;i++){
				a = a.sqrt();
				result = result.multiply(a.add(1.0));
			}
			
			result = z0.divide(result);
		}
			
		return result;
	}
	
	
	/**
	 * Computes the matrix logarithm of a 2x2 upper triangular matrix
	 * @param A (triangular matrix with complex entries)
	 * @return The matrix log of A.
	 */	
	private static FieldMatrix<Complex> logm2by2(FieldMatrix<Complex> A){
		Complex a1 = A.getEntry(0, 0);
		Complex a2 = A.getEntry(1, 1);
		Complex loga1 = (a1).log();//A(1,1);
		Complex loga2 = (a2).log();//A(2,2);
		
		Complex X_12;
		
		if(a1 == a2){
			X_12 = (A.getEntry(0,1)).divide(a1);
		}else if(a1.abs()< 0.5*a2.abs() || a2.abs()< 0.5*a1.abs()){
			X_12 = ((A.getEntry(0,1)).multiply(loga2.subtract(loga1))).divide(a2.subtract(a1));
		}else{
			Complex term1 = (((a2.subtract(a1)).divide(a2.add(a1))).atan()).multiply(2.0);
			Complex term2 = (((Complex.I).multiply(2*Math.PI)).multiply(unwinding(loga2.subtract(loga1))));
			X_12 = (A.getEntry(0,1)).multiply((term1.add(term2)).divide(a2.subtract(a1)));
		}
		
		Complex[][] arrayResult = {{loga1, X_12},{Complex.ZERO,loga2}};
		FieldMatrix<Complex> matrixResult = MatrixUtils.createFieldMatrix(arrayResult);
		return matrixResult;
	}
	/**
	 * 
	 * @param z
	 * @return the unwinding number of the complex number z
	 */
	private static Complex unwinding(Complex z){
		
		double u = Math.ceil((z.getImaginary()-Math.PI)/(2*Math.PI));	
		Complex result = new Complex(u,0);
		return result; 
	}
	
	public static void main(String[] args){
		Complex A11 = new Complex(0.4854,0.9157);
		Complex A12 = new Complex(0.1419,0.9595);
		Complex A21 = new Complex(0.8003,0.7922);
		Complex A22 = new Complex(0.4218,0.6557);
		
		Complex[][] data2 = {{A11,A12},{A21,A22}};
		
		FieldMatrix<Complex> myMatrix2 = MatrixUtils.createFieldMatrix(data2);
		
		System.out.println("Matrix logarithm via Java");
		long startTime = System.currentTimeMillis();
		System.out.println(MatrixLogarithm.compute(myMatrix2));
		long endTime = System.currentTimeMillis();
		long time = endTime-startTime;
		System.out.println("The computation required "+time +" milliseconds.");
		
		System.out.println("Matrix logarithm via Matlab");
		System.out.println("-0.474206446640814 - 0.098186165925212i  0.707064374568391 + 1.575697818629481i");
		System.out.println(" 1.754484780944842 + 0.970678044155572i -0.708900389443337 - 0.513001391672771i");
		
	}
}
