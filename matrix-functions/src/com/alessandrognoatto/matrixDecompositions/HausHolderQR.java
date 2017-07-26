package com.alessandrognoatto.matrixDecompositions;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexField;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.linear.FieldVector;
import org.apache.commons.math3.linear.MatrixUtils;

public class HausHolderQR {
	FieldMatrix<Complex> Q;
	FieldMatrix<Complex> R;
	
	/**
	 * This class performs the QR decomposition of a Matrix by means of Hausholder reflections.
	 * 
	 * See e.g. http://en.wikipedia.org/wiki/QR_decomposition#Using_Householder_reflections
	 * 
	 * @param A
	 */	
	public HausHolderQR(FieldMatrix<Complex> A){
		int m = A.getRowDimension();
		int n = A.getColumnDimension();
		
		FieldMatrix<Complex> eye = MatrixUtils.createFieldIdentityMatrix(ComplexField.getInstance(), m);
		
		FieldMatrix<Complex> Ak = A;
		FieldMatrix<Complex> Qt = eye;
		
		for(int i = 0; i<n; i++){
			
			//calculate the zapping vector for column i
			FieldVector<Complex> e = MatrixUtils.createFieldVector(eye.getColumn(i));
			FieldVector<Complex> a = MatrixUtils.createFieldVector(Ak.getColumn(i));
			
			for(int k = 0; k < i;k++){
				a.setEntry(k, Complex.ZERO);
			}
			
			
			//v = a + exp(1i*angle(a(i))) * norm(a) * e ;
			Complex[] arrayA = a.toArray();
			double absValueA = 0;
			for(int j = 0; j < arrayA.length; j++){
				absValueA += Math.pow(arrayA[j].abs(),2);
			}
			double normA = Math.sqrt(absValueA);
						
			Complex vv = (((Complex.I).multiply((a.getEntry(i).getArgument()))).exp()).multiply(normA);
			
			FieldVector<Complex> v = a.add(e.mapMultiply(vv));
			
			//explicit calculation and application of H
			Complex[] arrayV = v.toArray();
			double absValueV = 0;
			for(int j = 0; j < arrayV.length; j++){
				absValueV += Math.pow(arrayV[j].abs(),2);
			}
			Complex normV = new Complex(Math.sqrt(absValueV),0);
			
			// v = v./norm(v); mapDivideToSelf substitutes this with this/argument
			v.mapDivideToSelf(normV);
			
			
			Complex wnum = Complex.ZERO;
			Complex wden = Complex.ZERO;
			//w = (a'* v) / (v' * a); 
			for(int j = 0; j< a.getDimension(); j++){
				wnum = wnum.add((a.getEntry(j).conjugate()).multiply(v.getEntry(j)));
				wden = wden.add((v.getEntry(j).conjugate()).multiply(a.getEntry(j)));
			}
			
			Complex w = wnum.divide(wden);
			
			//I have to conjugate v element-wise by hand
			Complex[] dataConj = new Complex[v.getDimension()];
			for(int j = 0; j < dataConj.length; j++){
				dataConj[j] = (v.getEntry(j)).conjugate();
			}
			FieldVector<Complex> conjV = MatrixUtils.createFieldVector(dataConj);
		
			//H = I - (1 + w) * (v * v');
			FieldMatrix<Complex> H = eye.subtract((v.outerProduct(conjV)).scalarMultiply(w.add(1.0)));
			
			//Ak = H * Ak ;
			Ak = H.multiply(Ak);
			//System.out.println("Ak = "+Ak);
			//Qt = H * Qt;
			Qt = H.multiply(Qt);
			//System.out.println("Qt = "+ Qt);			
		}
		
		Qt = Qt.transpose();
		for(int j = 0; j < Qt.getColumnDimension();j++){
			for(int k = 0; k< Qt.getRowDimension();k++){
				Qt.setEntry(k, j, (Qt.getEntry(k, j)).conjugate());
			}
		}
		
		this.R = Ak;
		this.Q = Qt;		
	}
	
	public FieldMatrix<Complex> getQ(){
		return this.Q;
	}
	
	public FieldMatrix<Complex> getR(){
		return this.R;
	}
	
	/*public static void main(String[] args){
		
		Complex A11 = new Complex(0.4854,0.9157);
		Complex A12 = new Complex(0.1419,0.9595);
		Complex A21 = new Complex(0.8003,0.7922);
		Complex A22 = new Complex(0.4218,0.6557);
		
		Complex[][] data2 = {{A11,A12},{A21,A22}};
		FieldMatrix<Complex> myMatrix2 = MatrixUtils.createFieldMatrix(data2);
		
		HausHolderQR myQR = new HausHolderQR(myMatrix2);
		
		System.out.println("The matrix Q is");
		System.out.println(myQR.getQ());
		System.out.println("Matlab Result for Q is");
		System.out.println("-0.677198742456113 + 0.000000000000000i  0.702268388065559 + 0.219592746555210i");
		System.out.println("-0.702268388065558 + 0.219592746555210i -0.677198742456114 - 0.000000000000000i");
		
		System.out.println("The matrix R is");
		System.out.println(myQR.getR());
		
		System.out.println("Matlab Result for R is");
		System.out.println("-0.716776286736026 - 1.352187980560733i -0.248324343724324 - 1.202873795938215i");
		System.out.println("-0.000000000000000 - 0.000000000000000i  0.024708695018238 + 0.198627092184245i");
	}*/
}
