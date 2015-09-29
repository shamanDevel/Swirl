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
		
		double fx = (cos3*m3*(-1 + Nsq)*ox - m*nx2*ox - m*nx*ny*oy - m*nx*nz*oz + px - m*ny2*px - m*nz*px + m*nx*ny*py - 
     cos2*m2*((-2 + (2 + m)*nx2 + (1 + m)*ny2 + nz + m*nz)*ox + (-1 + ny2 + nz)*px + nx*(ny*(oy - py) + nz*(oz - pz))) + m*nx*nz*pz - 
     m2*(m*Nsq*Nsq*ox - ny2*ox - nz*ox + nx*nz*oz - nx2*px + nx*ny*(oy - py) - nx*nz*pz)*sin2 + 
     m*cos*((-1 + (1 + 2*m)*nx2 + m*(ny2 + nz))*ox + nx*ny*oy + m*nx*ny*oy + nx*nz*oz + m*nx*nz*oz - 2*px + ny2*px + m*ny2*px + 
        nz*px + m*nz*px - nx*ny*py - m*nx*ny*py - nx*nz*pz - m*nx*nz*pz + m2*Nsqsq*ox*sin2 + 
        m*(-1 + Nsq)*(nz*(-oy + py) + ny*(oz - pz))*sin) + 
     m*sin*((-1 + m*Nsq)*(nz*(oy - py) + ny*(-oz + pz)) + m2*(ny2*nz + nx2*(ny2 + nz))*ox*s2a));
		fx /= D;
		
		double fy = (-(m*nx*ny*ox) + cos3*m3*(-1 + Nsq)*oy - m*ny2*oy - m*ny*nz*oz + m*nx*ny*px + py - m*nx2*py - m*nz*py + m*ny*nz*pz - 
     cos2*m2*((-2 + (2 + m)*ny2 + (1 + m)*nz)*oy + ny*nz*oz + nx*ny*(ox - px) - py + nz*py + nx2*(oy + m*oy + py) - ny*nz*pz) - 
     m2*(m*nx4*oy - nz*oy + m*(ny2 + nz)*(ny2 + nz)*oy + nx2*(-1 + 2*m*(ny2 + nz))*oy + ny*nz*oz + nx*ny*(ox - px) - ny2*py - ny*nz*pz)*
      sin2 + m*cos*(-oy + ny2*oy + 2*m*ny2*oy + m*nz*oy + ny*nz*oz + m*ny*nz*oz + (1 + m)*nx*ny*(ox - px) - 2*py + nz*py + 
        m*nz*py + nx2*(py + m*(oy + py)) - ny*nz*pz - m*ny*nz*pz + m2*Nsqsq*oy*sin2 + 
        m*(-1 + Nsq)*(nz*(ox - px) + nx*(-oz + pz))*sin) + 
     m*sin*((-1 + m*Nsq)*(nz*(-ox + px) + nx*(oz - pz)) + m2*(ny2*nz + nx2*(ny2 + nz))*oy*s2a));
		fy /= D;
		
		double fz = (-(m*nx*nz*ox) - m*ny*nz*oy + cos3*m3*(-1 + Nsq)*oz - m*nz*oz + m*nx*nz*px + m*ny*nz*py + pz - m*nx2*pz - m*ny2*pz - 
     cos2*m2*(-2*oz + 2*nz*oz + m*nz*oz + nx*nz*(ox - px) + ny*nz*(oy - py) - pz + nx2*(oz + m*oz + pz) + ny2*(oz + m*oz + pz)) - 
     m2*(m*nx4*oz + m*ny4*oz + ny2*(-1 + 2*m*nz)*oz + nx2*(-1 + 2*m*(ny2 + nz))*oz + m*nz4*oz + nx*nz*(ox - px) + ny*nz*(oy - py) - 
        nz*pz)*sin2 + m*cos*(ny*nz*oy + m*ny*nz*oy - oz + m*ny2*oz + nz*oz + 2*m*nz*oz + (1 + m)*nx*nz*(ox - px) - ny*nz*py - 
        m*ny*nz*py - 2*pz + ny2*pz + m*ny2*pz + nx2*(pz + m*(oz + pz)) + m2*Nsqsq*oz*sin2 + 
        m*(-1 + Nsq)*(ny*(-ox + px) + nx*(oy - py))*sin) + 
     m*sin*((-1 + m*Nsq)*(ny*(ox - px) + nx*(-oy + py)) + m2*(ny2*nz + nx2*(ny2 + nz))*oz*s2a));
		fz /= D;
		
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
