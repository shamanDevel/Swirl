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
	private float alpha;
	Vector3f F;

	@Override
	public void setStartEnd(Frame start, Frame end) {
		this.start = start;
		this.end = end;
		System.out.println("Calculate Swirl P(t)=F+m^t*(P0-F)^'(alpha*t, N)");
		System.out.println(" Start frame: " + start);
		System.out.println(" End frame: " + end);

		//Calculate m
		this.m = end.I.length() / start.I.length();

		//Calculate N and alpha
		calcNAlpha();
		//calculate F
		calcF();
		//Check if the rotation axis N is flipped
		if (checkForAxisFlipped()) {
			System.out.println("flip N");
			N.negateLocal();
			calcF();
		}

		System.out.println(" scaling m=" + m);
		System.out.println(" rotation axis N=" + N);
		System.out.println(" rotation angle alpha=" + alpha + "  ( " + (alpha * 180 / Math.PI) + "° )");
		
		//Test
		Vector3f test = new Vector3f();
		test.addLocal(F);
		test.addScaleLocal(project(start.P.subtract(F), N), -m);
		test.addScaleLocal(start.P.subtract(F).subtract(project(start.P.subtract(F), N)), (float) (m*Math.cos(alpha)));
		test.addScaleLocal(N.cross(start.P.subtract(F).subtract(project(start.P.subtract(F), N))), (float) (-m*Math.sin(alpha)));
		System.out.println(" calculated P1=" + test+", expected P1="+end.P);

		//System.out.println(" lhs="+lhs);
		//System.out.println(" M="+M);
		System.out.println(" rotation center F=" + F);
	}
	
	private void calcF() {
		double nx = N.x;
		double ny = N.y;
		double nz = N.z;
		double norm = Math.sqrt(nx*nx + ny*ny + nz*nz);
		nx/=norm;
		ny/=norm;
		nz/=norm;
		double ox = start.P.x;
		double oy = start.P.y;
		double oz = start.P.z;
		double px = end.P.x;
		double py = end.P.y;
		double pz = end.P.z;
		double m2 = m*m;
		double m3 = m2*m;
		double nx2 = nx*nx;
		double nx4 = nx2*nx2;
		double ny2 = ny*ny;
		double ny4 = ny2*ny2;
		double nz2 = nz*nz;
		double nz4 = nz2*nz2;
		double Nsq = nx2 + ny2 + nz2; //should be one
		double Nsqsq = -nx2+nx4-ny2+ny4-nz2+nz4;
		double cos = Math.cos(alpha);
		double sin = Math.sin(alpha);
		double s2a = Math.sin(alpha*2);
		double cos2 = cos*cos;
		double cos3 = cos2*cos;
		double sin2 = sin*sin;
		
		double D = (1 - m*Nsq + m*(-1 + Nsq)*cos)*(1 + cos2*m2 - 2*m*cos + m2*Nsq*sin2);
		
		double fx = -((px - m*ox*Math.cos(alpha) - 2*m*px*Math.cos(alpha) + 2*Math.pow(m,2)*ox*Math.pow(Math.cos(alpha),2) + Math.pow(m,2)*px*Math.pow(Math.cos(alpha),2) - Math.pow(m,3)*ox*Math.pow(Math.cos(alpha),3) + 
       2*m*(1 + m)*Math.pow(ny,2)*(ox - px)*Math.pow(Math.sin(alpha/2.),2) + 2*m*(1 + m)*Math.pow(nz,2)*(ox - px)*Math.pow(Math.sin(alpha/2.),2) - 2*m*(1 + m)*nx*ny*(oy - py)*Math.pow(Math.sin(alpha/2.),2) - 
       2*m*(1 + m)*nx*nz*(oz - pz)*Math.pow(Math.sin(alpha/2.),2) - m*nz*(oy - py)*(-1 + m*Math.cos(alpha))*Math.sin(alpha) + m*ny*(oz - pz)*(-1 + m*Math.cos(alpha))*Math.sin(alpha))/Math.pow(-1 + m*Math.cos(alpha),3));
		
		double fy = (-py + m*oy*Math.cos(alpha) + 2*m*py*Math.cos(alpha) - 2*Math.pow(m,2)*oy*Math.pow(Math.cos(alpha),2) - Math.pow(m,2)*py*Math.pow(Math.cos(alpha),2) + Math.pow(m,3)*oy*Math.pow(Math.cos(alpha),3) + 
     2*m*(1 + m)*nx*ny*(ox - px)*Math.pow(Math.sin(alpha/2.),2) + 2*m*(1 + m)*Math.pow(ny,2)*(oy - py)*Math.pow(Math.sin(alpha/2.),2) + 2*m*(1 + m)*ny*nz*(oz - pz)*Math.pow(Math.sin(alpha/2.),2) - 
     m*nz*(ox - px)*(-1 + m*Math.cos(alpha))*Math.sin(alpha) + m*nx*(oz - pz)*(-1 + m*Math.cos(alpha))*Math.sin(alpha))/Math.pow(-1 + m*Math.cos(alpha),3);
		
		double fz = -((pz + Math.pow(m,2)*(2*oz + pz)*Math.pow(Math.cos(alpha),2) - Math.pow(m,3)*oz*Math.pow(Math.cos(alpha),3) - 2*m*(1 + m)*nx*nz*(ox - px)*Math.pow(Math.sin(alpha/2.),2) - 2*m*(1 + m)*ny*nz*(oy - py)*Math.pow(Math.sin(alpha/2.),2) - 
       2*m*(1 + m)*Math.pow(nz,2)*(oz - pz)*Math.pow(Math.sin(alpha/2.),2) + m*ny*(ox - px)*Math.sin(alpha) + m*nx*(-oy + py)*Math.sin(alpha) - m*Math.cos(alpha)*(oz + 2*pz + m*(ny*(ox - px) + nx*(-oy + py))*Math.sin(alpha)))/
     Math.pow(-1 + m*Math.cos(alpha),3));
		
		F = new Vector3f((float) fx, (float) fy, (float) fz);
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
		
//		//normalize signs
//		int sign1 = ((int) Math.signum(N1.x)) + ((int) Math.signum(N1.y)) + ((int) Math.signum(N1.z));
//		int sign2 = ((int) Math.signum(N2.x)) + ((int) Math.signum(N2.y)) + ((int) Math.signum(N2.z));
//		int sign3 = ((int) Math.signum(N3.x)) + ((int) Math.signum(N3.y)) + ((int) Math.signum(N3.z));
//		int sign = sign1 + sign2 + sign3;
//		if ((sign>0 && sign1<0) || (sign<0 && sign1>0)) {
//			N1.negateLocal();
//		}
//		if ((sign>0 && sign2<0) || (sign<0 && sign2>0)) {
//			N2.negateLocal();
//		}
//		if ((sign>0 && sign3<0) || (sign<0 && sign3>0)) {
//			N3.negateLocal();
//		}

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
		v.multLocal(X.dot(N));
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
		
//		Vector3f test = A.mult(result);
//		assert (test.equals(b));
		
		return result;
	}

	/**
	 * It can happen that the normal axis gets reversed and therefore the whole
	 * frame is reversed as well. Check this.
	 * @return {@code true} if the normal axis is flipped
	 */
	private boolean checkForAxisFlipped() {
		Vector3f v = new Vector3f();
		
		Vector3f refP = new Vector3f();
		interpolate(1, start.P, refP);
		
		Vector3f tmp = start.P.add(start.I);
		interpolate(1, tmp, v);
		v.subtractLocal(refP);
		if (!v.equals(end.I)) { //includes an epsilon
			return true;
		}

		tmp = start.P.add(start.J);
		interpolate(1, tmp, v);
		v.subtractLocal(refP);
		if (!v.equals(end.J)) { //includes an epsilon
			return true;
		}
		
		tmp = start.P.add(start.K);
		interpolate(1, tmp, v);
		v.subtractLocal(refP);
		if (!v.equals(end.K)) { //includes an epsilon
			return true;
		}
		
		return false;
	}
	
	@Override
	public void interpolate(float t, Frame toSet) {
		interpolate(t, start.P, toSet.P);
		Vector3f refP = toSet.P;

		Vector3f tmp = start.P.add(start.I);
		interpolate(t, tmp, toSet.I);
		toSet.I.subtractLocal(refP);

		tmp = start.P.add(start.J);
		interpolate(t, tmp, toSet.J);
		toSet.J.subtractLocal(refP);
		
		tmp = start.P.add(start.K);
		interpolate(t, tmp, toSet.K);
		toSet.K.subtractLocal(refP);

//		System.out.println("Frame at t="+t+": "+toSet);
	}

	private void interpolate(float t, Vector3f P0, Vector3f store) {
		store.set(F);
		Vector3f FP0 = P0.subtract(F);
		Vector3f rot = rotate(alpha*t, FP0);
		store.addScaleLocal(rot, (float) Math.pow(m, t));
	}

	/**
	 * Traditional rotation
	 * @param angle
	 * @param X
	 * @return 
	 */
	private Vector3f rotate(float angle, Vector3f X) {
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

	@Override
	public String debugString() {
		StringBuilder str = new StringBuilder();
		str.append("Swirl Interpolation P(t)=F+(m^t)*FPo^(alpha*t,N)\n");
		str.append("scaling m=").append(m).append("\n");
		str.append("rotation axis N=").append(N).append("\n");
		str.append("rotation angle alpha=").append(alpha).append(" (").append(alpha * 180 / Math.PI).append("°)\n");
		str.append("rotation center F=").append(F);
		return str.toString();
	}
}
