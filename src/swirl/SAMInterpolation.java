/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package swirl;

import org.apache.commons.math3.exception.MathArithmeticException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

/**
 * Implements the SAM-Interpolation
 * @author Sebastian Weiss
 */
public class SAMInterpolation implements Interpolation {
	private Frame startFrame;
	private Frame endFrame;
	private Matrix4f A0 = new Matrix4f();
	private Matrix4f A1 = new Matrix4f();
	private Matrix4f deltaA = new Matrix4f();
	private Matrix4f V = new Matrix4f();
	private Matrix4f Vt = new Matrix4f();
	private Matrix4f D = new Matrix4f();
	private boolean solved = false;

	@Override
	public void setStartEnd(Frame start, Frame end) {
		startFrame = start;
		endFrame = end;
		A0.fromFrame(start.P, start.J, start.K, start.I.negate());
		A1.fromFrame(end.P, end.J, end.K, end.I.negate());
		
		//Calculate transition matrix and diagonalize it
		A0.invert(deltaA);
		deltaA.multLocal(A1);
		//diagonalize
		Array2DRowRealMatrix matrix = new Array2DRowRealMatrix(4, 4);
		copy(deltaA, matrix);
		EigenDecomposition eigen = null;
		try {
			eigen = new EigenDecomposition(matrix);
		} catch (MathArithmeticException mathArithmeticException) {
			solved = false;
			return;
		}
		solved = true;
		copy(eigen.getV(), V);
		copy(eigen.getVT(), Vt);
		copy(eigen.getD(), D);
		
		System.out.println("A0: ");
		System.out.println(A0);
		System.out.println("A1: ");
		System.out.println(A1);
		System.out.println("deltaA: ");
		System.out.println(deltaA);
		System.out.println("V: ");
		System.out.println(V);
		System.out.println("D: ");
		System.out.println(D);
		
		Matrix4f test = V.mult(D).multLocal(Vt);
		assert (test.almostEqual(deltaA));
	}
	private void copy(Matrix4f source, RealMatrix target) {
		for (int row=0; row<4; ++row) {
			for (int column=0; column<4; ++column) {
				target.setEntry(row, column, source.get(row, column));
			}
		}
	}
	private void copy(RealMatrix source, Matrix4f target) {
		for (int row=0; row<4; ++row) {
			for (int column=0; column<4; ++column) {
				target.set(row, column, (float) source.getEntry(row, column));
			}
		}
	}
	
	private Matrix4f getPowerMatrix(float power) {
		Matrix4f mat = D.power(power);
		mat = V.mult(mat);
		mat.multLocal(Vt);
		return mat;
	}

	@Override
	public void interpolate(float t, Frame toSet) {
		if (!solved) {
			return;
		}
		Matrix4f M = getPowerMatrix(t);
		M.mult(startFrame.P, toSet.P);
		M.mult(startFrame.P.add(startFrame.I), toSet.I);
		M.mult(startFrame.P.add(startFrame.J), toSet.J);
		M.mult(startFrame.P.add(startFrame.K), toSet.K);
		toSet.I.subtractLocal(toSet.P);
		toSet.J.subtractLocal(toSet.P);
		toSet.K.subtractLocal(toSet.P);
//		System.out.println("Frame at t="+t+": "+toSet);
	}

	@Override
	public void debugDraw(Swirl swirl) {
	}

	@Override
	public String debugString() {
		StringBuilder str = new StringBuilder();
		str.append("SAM-Interpolation\n");
		str.append("matrix between start and end frame:\n");
		str.append(deltaA).append("\n");
		if (solved) {
			str.append("diagonalized matrix:\n");
			str.append(D);
		} else {
			str.append("unable to diagonalize matrix!");
		}
		return str.toString();
	}
	
}
