/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package swirl;

public class SwirlInterpolation implements Interpolation {

	private static final double EPSILON = 0.00001;

	private Frame start;
	private Frame end;

	private float m;
	Vector3f N = new Vector3f();
	private double alpha;
	Vector3f F;

	@Override
	public void setStartEnd(Frame start, Frame end) {
		this.start = start;
		this.end = end;
		System.out.println("Calculate Swirl P(t)=F+m^t*(P0-F)^(alpha*t, N)");
		System.out.println(" Start frame: " + start);
		System.out.println(" End frame: " + end);

		//Calculate m
		this.m = end.I.length() / start.I.length();

		//Calculate N and alpha
		calcNAlpha();

		System.out.println(" scaling m=" + m);
		System.out.println(" rotation axis N=" + N);
		System.out.println(" rotation angle alpha=" + alpha + "  ( " + (alpha * 180 / Math.PI) + "° )");

		//Calculate F
		Vector3f lhs = new Vector3f(end.P);
		Vector3f P0mP0N = new Vector3f(start.P);
		P0mP0N.subtractLocal(project(start.P, N));
		Vector3f lhsA = project(start.P, N);
		Vector3f lhsB = new Vector3f(P0mP0N);
		lhsB.multLocal((float) Math.cos(alpha));
		Vector3f lhsC = P0mP0N.cross(N);
		lhsC.multLocal((float) Math.sin(alpha));
		lhsA.addLocal(lhsB);
		lhsA.addLocal(lhsC);
		lhsA.multLocal(m);
		lhs.subtractLocal(lhsA);

		double a = 1 - m * Math.cos(alpha);		// It's a kind of magic.
		double b = m * (1 + Math.cos(alpha));	// It's a kind of magiiiic,
		double c = -m * Math.sin(alpha);		// magic,
		Matrix3f M = new Matrix3f();		// magic,
		M.m00 = (float) (a + b * N.x * N.x);			//  M
		M.m01 = (float) (b * N.x * N.y - c * N.z);		//  A
		M.m02 = (float) (b * N.x * N.z + c * N.y);		//  A
		M.m10 = (float) (b * N.y * N.x + c * N.z);		//  G
		M.m11 = (float) (a + b * N.y * N.y);				//  I
		M.m12 = (float) (b * N.y * N.z + c * N.y);		//  I
		M.m20 = (float) (b * N.z * N.x - c * N.y);		//  I
		M.m21 = (float) (b * N.z * N.y + c * N.x);		//  C
		M.m22 = (float) (a + b * N.z * N.z);				//  !
		F = solve(M, lhs);					// Guitar solo ...

		//System.out.println(" lhs="+lhs);
		//System.out.println(" M="+M);
		System.out.println(" rotation center F=" + F);
	}

	private boolean calcNAlpha() {
		Vector3f I0 = start.I.normalize();
		Vector3f J0 = start.J.normalize();
		Vector3f K0 = start.K.normalize();
		Vector3f I1 = end.I.normalize();
		Vector3f J1 = end.J.normalize();
		Vector3f K1 = end.K.normalize();

		Vector3f dI = I1.subtract(I0);
		Vector3f dJ = J1.subtract(J0);
		Vector3f dK = K1.subtract(K0);

		Vector3f N1 = dI.cross(dJ);
		Vector3f N2 = dI.cross(dK);
		Vector3f N3 = dJ.cross(dK);
		float l1 = N1.lengthSquared();
		float l2 = N2.lengthSquared();
		float l3 = N3.lengthSquared();
		System.out.println("N1=" + N1 + "  N2=" + N2 + "  N3=" + N3);
		System.out.println("N1=" + N1.normalize() + "  N2=" + N2.normalize() + "  N3=" + N3.normalize() + "  (vectors normalized)");

		//use vector with the greatest magnitude
		if (l1 < EPSILON && l2 < EPSILON && l3 < EPSILON) {
			return false; //no rotation
		}
		if (l1 >= l2 && l1 >= l3) {
			N = N1;
		} else if (l2 >= l1 && l2 >= l3) {
			N = N2;
		} else {
			N = N3;
		}
		N.normalizeLocal();

		//project vectors in the plane of N
		Vector3f pI0 = projectOnPlane(I0, N).normalizeLocal();
		Vector3f pJ0 = projectOnPlane(J0, N).normalizeLocal();
		Vector3f pK0 = projectOnPlane(K0, N).normalizeLocal();
		Vector3f pI1 = projectOnPlane(I1, N).normalizeLocal();
		Vector3f pJ1 = projectOnPlane(J1, N).normalizeLocal();
		Vector3f pK1 = projectOnPlane(K1, N).normalizeLocal();
		float angle1 = pI0.angleBetween(pI1);
		float angle2 = pJ0.angleBetween(pJ1);
		float angle3 = pK0.angleBetween(pK1);
		System.out.println("a1= " + angle1 + "  a2=" + angle2 + "  a3=" + angle3);
		//choose greatest angle
		if (Math.abs(angle1) >= Math.abs(angle2) && Math.abs(angle1) >= Math.abs(angle3)) {
			alpha = angle1;
		} else if (Math.abs(angle2) >= Math.abs(angle1) && Math.abs(angle2) >= Math.abs(angle3)) {
			alpha = angle2;
		} else {
			alpha = angle3;
		}

		return true;
	}

	private Vector3f project(Vector3f X, Vector3f N) {
		Vector3f v = new Vector3f(N);
		v.multLocal(X.x * N.x + X.y * N.y + X.z * N.z);
		return v;
	}

	private Vector3f projectOnPlane(Vector3f X, Vector3f N) {
		Vector3f v = N.mult(X.dot(N));
		v.set(X.x - v.x, X.y - v.y, X.z - v.z);
		return v;
	}

	/**
	 * Solves Ax=b
	 *
	 * @param A the matrix
	 * @param b the right hand side of the equation
	 * @return the solution x or {@code null} if not soveable
	 */
	private Vector3f solve(Matrix3f A, Vector3f b) {
		//Solve using cramers rule
		Vector3f result = new Vector3f();
		float detA = A.determinant();
		if (Math.abs(detA) < EPSILON) {
			System.err.println("determinant is " + detA + ", not solveable");
		}
		Matrix3f A1 = new Matrix3f(A);
		A1.setColumn(0, b);
		result.x = A1.determinant() / detA;
		Matrix3f A2 = new Matrix3f(A);
		A2.setColumn(1, b);
		result.y = (A2.determinant() / detA);
		Matrix3f A3 = new Matrix3f(A);
		A3.setColumn(2, b);
		result.z = A3.determinant() / detA;
		return result;
	}

	@Override
	public void interpolate(float t, Frame toSet) {
		interpolate(t, start.P, toSet.P);

		Vector3f tmp = start.P.add(start.I);
		interpolate(t, tmp, toSet.I);
		toSet.I.subtractLocal(toSet.P);

		tmp = start.P.add(start.J);
		interpolate(t, tmp, toSet.J);
		toSet.J.subtractLocal(toSet.P);
		
		tmp = start.P.add(start.K);
		interpolate(t, tmp, toSet.K);
		toSet.K.subtractLocal(toSet.P);

//		System.out.println("Frame at t="+t+": "+toSet);
	}

	private void interpolate(double t, Vector3f P0, Vector3f store) {
		store.set(F);
		Vector3f FP0 = P0.subtract(F);
		Vector3f rot = rotate(t * alpha, FP0);
		store.addScaleLocal(rot, (float) Math.pow(m, t));
	}

	private Vector3f rotate(double angle, Vector3f X) {
		Vector3f W = project(X, N);
		Vector3f U = X.subtract(W);
		Vector3f result = new Vector3f(W);
		result.addScaleLocal(U, (float) Math.cos(angle));
		result.addScaleLocal(N.cross(U), -(float) Math.sin(angle));
		return result;
	}

	@Override
	public void debugDraw(Swirl swirl) {
		swirl.pushMatrix();
		swirl.translate(F.x * Swirl.SCALE, F.y * Swirl.SCALE, F.z * Swirl.SCALE);
		swirl.noStroke();
		swirl.fill(Swirl.grey);
		swirl.sphere(0.15f * Swirl.SCALE);
		swirl.collar(Pv3D.P(), Pv3D.V(N.x * 50 * Swirl.SCALE, N.y * 50 * Swirl.SCALE, N.z * 50 * Swirl.SCALE),
				0.05f * Swirl.SCALE, 0.05f * Swirl.SCALE);
		swirl.collar(Pv3D.P(), Pv3D.V(N.x * -50 * Swirl.SCALE, (N.y * -50 * Swirl.SCALE), (N.z * -50 * Swirl.SCALE)),
				0.05f * Swirl.SCALE, 0.05f * Swirl.SCALE);
		swirl.popMatrix();
	}

}
