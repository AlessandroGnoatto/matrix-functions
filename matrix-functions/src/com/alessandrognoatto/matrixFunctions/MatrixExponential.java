/**
 * 
 */
package com.alessandrognoatto.matrixFunctions;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.complex.*;

/**
 * This class computes the matrix exponential for a square matrix with complex entries.
 * 
 * The matrix exponential is computed by using a scaling and squaring algorithm with a Pade approximation.
 * 
 * Reference:
 *  N. J. Higham, The scaling and squaring method for the matrix
 *  exponential revisited. SIAM J. Matrix Anal. Appl.,
 *  26(4) (2005), pp. 1179-1193.
 *  
 *
 * @author Alessandro Gnoatto
 * @param A complex Matrix i.e. FieldMatrix<Complex>
 * @return The matrix exponential of the complex Matrix
 */
public class MatrixExponential {
	
	public static FieldMatrix<Complex> compute(FieldMatrix<Complex> A){
		
		int n = A.getRowDimension();
		//line 9 of algorithm 2.3 in Higham (2005)
		int[] m_vals = {3, 5, 7, 9, 13};
		
		//See Table 2.3 In Higham (2005) page 1186
		double[] theta  = {1.495585217958292E-002, //  m_vals = 3
		                   2.539398330063230E-001,  // m_vals = 5
		                   9.504178996162932E-001,  // m_vals = 7
		                   2.097847961257068E+000,  // m_vals = 9
		                   5.371920351148152E+000};// m_vals = 13
		
		//compute the 1-norm of the matrix A
		double normA = computeOneNormComplexMatrix(A);

		FieldMatrix<Complex> F = MatrixUtils.createFieldMatrix(ComplexField.getInstance(),n,n);
		
		if(normA<=theta[theta.length-1]){
			 //no scaling and squaring is required.
			for(int i = 0; i<m_vals.length;i++){
				
				if(normA<=theta[i]){
					
					F = getPadeApproximantOfDegree(m_vals[i],A);
					
					break;
				}
			}
		} else {
			double t = frexp(normA/theta[theta.length-1]).mantissa;
			int s = frexp(normA/theta[theta.length-1]).exponent;
			
			if(t == 0.5){
				s = s-1;
			}
			System.out.println(s);
			Complex scalingFactor = new Complex(Math.pow(2, -s),0);
			A = A.scalarMultiply(scalingFactor);
			System.out.println(A);
			F = getPadeApproximantOfDegree(m_vals[m_vals.length-1],A);
			System.out.println("Sono qui");
			System.out.println(F);
			for(int k = 1; k <= s;k++){
				F = F.multiply(F);
			}
		}
		
		return F;
	}
	
	
			
	/**
	 * @param A matrix of complex number FieldMatrix<Complex>
	 * @return This method returns the one-norm of a complex matrix. For a matrix A, in Matlab this is given by norm(A,1).
	 * The function norm(A,1) is equivalent to the matlab instruction max(sum(abs(A))), where sum computes the sum over
	 * all columns of A i.e. take the absolute value of each entry of A, create a row vector containing the sum of all
	 * entries for each column in absolute value and finally take the biggest value in the resulting row vector.
	 */
	private static double computeOneNormComplexMatrix(FieldMatrix<Complex> A){
		// get the number of columns of the matrix
		int n = A.getColumnDimension();
		int m = A.getRowDimension();
		
		double[] partialSums = new double[n];
		
		for(int i=0; i<n;i++){
			
			double sum = 0;
			
			for(int j = 0; j<m;j++){
				sum += ((A.copy()).getEntry(j, i)).abs();
			}
			
			partialSums[i] = sum;
		}
		
		//sort element of the array
		Arrays.sort(partialSums);
		//return the largest element which is stored at the end
		return partialSums[partialSums.length-1];
	}
	



	private static FieldMatrix<Complex> getPadeApproximantOfDegree(int m, FieldMatrix<Complex> A){
		
		int n = A.getColumnDimension();
		Complex[] c = getPadeCoefficients(m);
		
		FieldMatrix<Complex> eye = MatrixUtils.createFieldIdentityMatrix(ComplexField.getInstance(), n);
		
		//initialize the result
		FieldMatrix<Complex> F = MatrixUtils.createFieldMatrix(ComplexField.getInstance(),n,n);
		
		//Initialize two matrices to zero
		FieldMatrix<Complex> U = MatrixUtils.createFieldMatrix(ComplexField.getInstance(),n,n);
		FieldMatrix<Complex> V = MatrixUtils.createFieldMatrix(ComplexField.getInstance(),n,n);
		
		FieldMatrix<Complex> negateU;
		FieldMatrix<Complex> firstTerm;
		FieldMatrix<Complex> invertedFirstTerm;
		
		switch (m){
			case 3:
			case 5:
			case 7:
			case 9:
				
				int index = (int) (Math.ceil((m+1)/2.0));
				
				ArrayList<FieldMatrix<Complex>> myList = new ArrayList<FieldMatrix<Complex>>();
				
				myList.add(eye);
				
				FieldMatrix<Complex> powerA = A;
				myList.add(powerA.multiply(A));
				
				/* Corresponding Matlab code
				 * 
				 * for j = 3:ceil((m+1)/2)
				 * 		Apowers{j} = Apowers{j-1}*Apowers{2};
				 * end
				 */
				for(int i = 2; i < index; i++){
					myList.add((myList.get(i-1)).multiply(myList.get(1)));
				}
				
				
				/* Corresponding Matlab code
				 * for j = m+1:-2:2
				 * 		U = U + c(j)*Apowers{j/2};
                	 * end
				 */
				for(int j = m+1; j>=2; j = j - 2){
					
					U = U.add((myList.get(j/2 - 1)).scalarMultiply(c[j-1]));
					
				}
				U = U.multiply(A);
				/* Corresponding Matlab code
				 *  
				 *  for j = m:-2:1
				 * 		V = V + c(j)*Apowers{(j+1)/2};
				 *  end
				 */			
				for(int j = m; j>= 1; j = j - 2){
					V = V.add((myList.get((j+1)/2-1).scalarMultiply(c[j-1])));
				}
				
				
				// F = (-U+V)\(U+V);
				negateU = MatrixUtils.createFieldMatrix(ComplexField.getInstance(),n,n);
				negateU = negateU.subtract(U);					
				firstTerm = (negateU.add(V));
				invertedFirstTerm = new FieldLUDecomposition<Complex>(firstTerm).getSolver().getInverse();
				F = invertedFirstTerm.multiply((U.add(V)));				
				break;
			case 13:	
				// A2 = A*A;
				FieldMatrix<Complex> A2 = (A.copy()).multiply(A.copy());
				
				//A4 = A2*A2;
				FieldMatrix<Complex> A4 = (A2.copy()).multiply(A2.copy());
				
				//A6 = A2*A4;
				FieldMatrix<Complex> A6 = (A4.copy()).multiply(A2.copy());
				
				/*
				 *	U = A * (A6*(c(14)*A6 + c(12)*A4 + c(10)*A2) ...
				 *		 + c(8)*A6 + c(6)*A4 + c(4)*A2 + c(2)*eye(n) ); 
				 */				
				FieldMatrix<Complex> term1 = (A6.copy()).multiply((((A6.copy()).scalarMultiply(c[13])).add((A4.copy()).scalarMultiply(c[11]))).add((A2.copy()).scalarMultiply(c[9])));
				FieldMatrix<Complex> term2 = ((((A6.copy()).scalarMultiply(c[7])).add((A4.copy()).scalarMultiply(c[5]))).add((A2.copy()).scalarMultiply(c[3]))).add((eye.copy()).scalarMultiply(c[1]));				
				U = (A.copy()).multiply((term1).add(term2));
				
				/*
				 *  V = A6*(c(13)*A6 + c(11)*A4 + c(9)*A2) ...
				 *  		 + c(7)*A6 + c(5)*A4 + c(3)*A2 + c(1)*eye(n);
				 */
				FieldMatrix<Complex> term3 = (((A6.copy()).scalarMultiply(c[12])).add((A4.copy()).scalarMultiply(c[10]))).add((A2.copy()).scalarMultiply(c[8]));
				FieldMatrix<Complex> term4 = ((((A6.copy()).scalarMultiply(c[6])).add((A4.copy()).scalarMultiply(c[4]))).add((A2.copy()).scalarMultiply(c[2]))).add((eye.copy()).scalarMultiply(c[0]));(((A6.copy()).scalarMultiply(c[6])).add((A4.copy()).scalarMultiply(c[4]))).add((eye.copy()).scalarMultiply(c[0]));
				V = ((A6.copy()).multiply(term3)).add(term4);
				
				// F = (-U+V)\(U+V);
				negateU = MatrixUtils.createFieldMatrix(ComplexField.getInstance(),n,n);
				negateU = negateU.subtract(U);					
				firstTerm = (negateU.add(V));
				invertedFirstTerm = new FieldLUDecomposition<Complex>(firstTerm).getSolver().getInverse();
				F = invertedFirstTerm.multiply((U.add(V)));		
				break;
		}
		

		return F;
	}
	
	
	private static Complex[] getPadeCoefficients(int m){
		ArrayList<Double> al = new ArrayList<Double>();
		
		switch (m){
			case 3:
				Double[] array3 = {120.0, 60.0, 12.0, 1.0};
				al.addAll(Arrays.asList(array3));
				break;
			case 5:
				Double[] array5 = {30240.0, 15120.0, 3360.0, 420.0, 30.0, 1.0};
				al.addAll(Arrays.asList(array5));
				break;
			case 7:
				Double[] array7 = {17297280.0, 8648640.0, 1995840.0, 277200.0, 25200.0, 1512.0, 56.0, 1.0};
				al.addAll(Arrays.asList(array7));
				break;
			case 9:
				Double[] array9 = {17643225600.0, 8821612800.0, 2075673600.0, 302702400.0, 30270240.0,
							2162160.0, 110880.0, 3960.0, 90.0, 1.0};
				al.addAll(Arrays.asList(array9));
				break;
			case 13:
				Double[] array13 = {64764752532480000.0, 32382376266240000.0, 7771770303897600.0,
                        1187353796428800.0,  129060195264000.0,   10559470521600.0,
                        670442572800.0,      33522128640.0,       1323241920.0,
                        40840800.0,          960960.0,            16380.0,  182.0,  1.0};	
				al.addAll(Arrays.asList(array13));
				break;
		}
		
		Complex resu[] = new Complex[al.size()];
				
		for(int i = 0; i< al.size(); i++){
			Complex c = new Complex(al.get(i),0);
			resu[i] = c;
		}
		return resu;
	}
	
	
	//just some tests
	public static void main(String[] Args){
		
		//FieldMatrix<Complex> eye = MatrixUtils.createFieldIdentityMatrix(ComplexField.getInstance(), 2);
		//System.out.println(eye.getEntry(0, 0));
		
		//System.out.println("Matrix Exponential");
		//System.out.println(MatrixExponential.compute(eye));
		
		Complex A11 = new Complex(0.4854,0.9157);
		Complex A12 = new Complex(0.1419,0.9595);
		Complex A21 = new Complex(0.8003,0.7922);
		Complex A22 = new Complex(0.4218,0.6557);
		
		Complex[][] data = {{A11,A12},{A21,A22}};
		FieldMatrix<Complex> myMatrix = MatrixUtils.createFieldMatrix(data);
		System.out.println("The complex matrix A is");
		System.out.println(myMatrix);
		
		System.out.println("The matrix exponential of A is");
		
		long startTime = System.currentTimeMillis();
		System.out.println(MatrixExponential.compute(myMatrix));
		long endTime = System.currentTimeMillis();
		
		System.out.println("That took " + (endTime - startTime) + " milliseconds");
		
		System.out.println("Matlab answer for the same matrix is");
		System.out.println(" 0.167128030866949 + 1.315269960148860i -0.976979580583788 + 0.960988622386044i");
		System.out.println("-0.237742415131101 + 1.573148823645988i  0.410921842402859 + 1.026162507008104i");
		
	}
	
	/**
	 * [F,E] = log2(X) for each element of the real array X, returns an
	 * array F of real numbers, usually in the range 0.5 <= abs(F) < 1,
	 * and an array E of integers, so that X = F .* 2.^E.  Any zeros in X
	 * produce F = 0 and E = 0.  This corresponds to the ANSI C function
	 *  frexp()
	 *  
	 *  Example [t, s] = log2(10)
	 *  then 
	 *  t = 0.6250
	 *  s = 4
	 *  and
	 *  
	 *  t*2^s = 10
	 *  
	 *  The following code is taken from stackoverflow.com 
	 *  Question url: http://stackoverflow.com/questions/1552738/is-there-a-java-equivalent-of-frexp
	 *  Author of solution: Jay R.
	 *  Author url: http://stackoverflow.com/users/5074/jay-r
	 *  
	 */
	public static class FRexpResult
	{
	   public int exponent = 0;
	   public double mantissa = 0.;
	}

	public static FRexpResult frexp(double value)
	{
	   final FRexpResult result = new FRexpResult();
	   long bits = Double.doubleToLongBits(value);
	   double realMant = 1.;

	   // Test for NaN, infinity, and zero.
	   if (Double.isNaN(value) || 
	       value + value == value || 
	       Double.isInfinite(value))
	   {
	      result.exponent = 0;
	      result.mantissa = value;
	   }
	   else
	   {

	      boolean neg = (bits < 0);
	      int exponent = (int)((bits >> 52) & 0x7ffL);
	      long mantissa = bits & 0xfffffffffffffL;

	      if(exponent == 0)
	      {
	         exponent++;
	      }
	      else
	      {
	         mantissa = mantissa | (1L<<52);
	      }

	      // bias the exponent - actually biased by 1023.
	      // we are treating the mantissa as m.0 instead of 0.m
	      //  so subtract another 52.
	      exponent -= 1075;
	      realMant = mantissa;

	      // normalize
	      while(realMant > 1.0) 
	      {
	         mantissa >>= 1;
	         realMant /= 2.;
	         exponent++;
	      }

	      if(neg)
	      {
	         realMant = realMant * -1;
	      }

	      result.exponent = exponent;
	      result.mantissa = realMant;
	   }
	   return result;
	}

}
