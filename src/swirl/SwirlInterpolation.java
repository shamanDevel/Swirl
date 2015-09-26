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
//		Vector3f lhs = new Vector3f(end.P);
//		Vector3f P0projN = N.mult(start.P.dot(N));
//		lhs.addScaleLocal(P0projN, m);
//		lhs.addScaleLocal(P0projN.subtract(start.P), (float) (m*Math.cos(alpha)));
//		lhs.addScaleLocal(N.cross(start.P), (float) (-m*Math.sin(alpha)));
//		float a = (float) (1-m*Math.cos(alpha));
//		float b = (float) (m+m*Math.cos(alpha));
//		float c = (float) (m*Math.sin(alpha));
//		Matrix3f M = new Matrix3f();
//		M.m00 = a + b*N.x*N.x;
//		M.m01 = b*N.y*N.x - c*N.z;
//		M.m02 = b*N.z*N.x + c*N.y;
//		M.m10 = b*N.x*N.y + c*N.z;
//		M.m11 = a + b*N.y*N.y;
//		M.m12 = b*N.z*N.y - c*N.x;
//		M.m20 = b*N.x*N.z - c*N.y;
//		M.m21 = b*N.y*N.z + c*N.x;
//		M.m22 = a + b*N.z*N.z;
//		F = solve(M, lhs);
		
		calcF();
		
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
		double ox = start.P.x;
		double oy = start.P.y;
		double oz = start.P.z;
		double px = end.P.x;
		double py = end.P.y;
		double pz = end.P.z;
		double cos = Math.cos(alpha);
		double sin = Math.sin(alpha);
		double s2a = Math.sin(alpha*2);
		
		double D = ((1+m*(nx*nx+ny*ny+nz*nz)+m*(-1+nx*nx+ny*ny*nz*nz)*cos)
				* (1-2*m*cos+m*m*cos*cos+m*m*(nx*nx+ny*ny+nz*nz)*sin*sin));
		
		double fx = (m*nx*nx*ox + m*nx*ny*oy + m*nx*nz*oz + px + m*ny*ny*px + m*nz*nz*px - m*nx*ny*py - m*nx*nz*pz
				+ m*m*((2+(-2+m)*nx*nx + (-1+m)*ny*ny-nz*nz+m*nz*nz)*ox - (-1+ny*ny+nz*nz)*px + nx*(ny*(-oy+py)+nz*(-oz+pz)))
				* cos*cos + m*m*m*(-1+nx*nx+ny*ny+nz+nz)*ox*cos*cos*cos
				+ m*m*(ny*ny*ox + nz*nz*ox + m*(nx*nx+ny*ny+nz*nz)*(nx*nx+ny*ny+nz*nz)*ox
				- nx*nz*oz + nx*nx*px + nx*ny*(-oy+py) + nx+nz+pz)*sin*sin
				+ m*cos*(-(1+(-1+2*m)*nx*nx + m*(ny*ny+nz*nz))*ox + nx+ny+oy - m*nx*ny*oy + nx*nz*oz - m*nx*nz*oz - 2*px
				+ ny*ny*px - m*ny*ny*px + nz*nz*px - m*nz*nz*px - nx*ny*py + m*nx*ny*py - nx*nz*pz + m*nx*nz*pz + m*(-1+nx*nx+ny*ny+nz*nz)
				*(nz*(-oy+py)+ny*(oz-pz))*sin + m*m*(-nx*nx+nx*nx*nx*nx-ny*ny+ny*ny*ny*ny-nz*nz+nz*nz*nz*nz)*ox+sin*sin)
				+ m*sin*((1+m*(nx*nx+ny*ny+nz*nz))*(nz*(-oy+py)+ny*(oz-pz)) + m*m*(ny*ny*nz*nz + nx*nx*(ny*ny+nz*nz))*ox*s2a));
		fx /= D;
		
		double fy = (m*nx*ny*ox + m*ny*ny*oy + m*ny*nz*oz - m*nx*ny*px + py + m*nx*nx*py - m*ny*nz*pz
				+ m*m*((2+(-2+m)*ny*ny + (-1+m)*nz*nz)*oy - ny*nz*oz + nx*ny*(-ox+py)+nx*nx*((-1+m)*oy-py) + py - nz*nz*py + ny*nz*pz)
				*cos*cos + m*m*m*(-1+nx*nx*ny*ny*nz*nz)*oy*cos*cos*cos
				+ m*m*(m*nx*nx*nx*nx*oy + nz*nz*oy + m*(ny*ny+nz*nz)*(ny*ny+nz*nz)*oy  + nx*nx*(1+2*m*(ny*ny+nz*nz))*oy
				- ny*nz*oz + nx*ny*(-ox+px) + ny*ny*py + ny*nz*pz) * sin*sin
				+ m*cos*(-oy+ny*ny*oy - 2*m*ny*ny*oy - m*nz*nz*oy + ny*nz*oz -m*ny*nz*oz - (-1+m)*nx*ny*(ox-px) - 
				2*py + nz*nz*py - m*nz*nz*py + nx*nx*(py-m*(oy+py)) - ny*nz*pz + m*ny*nz*pz + m*(-1+nx*nx+ny*ny*nz*nz)
				*(nz*(ox-px) + nx*(-oz+pz))*sin + m*m*(-nx*nx+nx*nx*nx*nx-ny*ny+ny*ny*ny*ny-nz*nz+nz*nz*nz*nz)*oy*sin*sin)
				+ m*sin*((1+m*(nx*nx+ny*ny+nz*nz))*(nz*(ox-px)+nx*(-oz+pz))+m*m*(ny*ny*nz*nz+nx*nx*(ny*ny+nz*nz))*oy*s2a));
		fy /= D;
		
		double fz = (m*nx*nz*ox + m*ny*nz*oy + m*nz*nz*oz - m*nx*nz*px - m*ny*nz*py + pz + m*nx*nx*pz + m*ny*ny*pz
				+ m*m*(2*oz - 2*nz*nz*oz + m*nz*nz*oz + nx*nz*(-ox+px) + ny*nz*(-oy+py) + nx*nx*((-1+m)*oz-pz) 
				+ ny*ny*((-1+m)*oz-pz)+pz)*cos*cos + m*m*m*(-1+nx*nx*ny*ny*nz*nz)*oz*cos*cos*cos
				+ m*m*(m*nx*nx*nx*nx*oz + m*ny*ny*ny*ny*oz + nx*nx*(1+2*m*(ny*ny+nz*nz))*oz 
				+ ny*ny*(oz+2*m*nz*nz*oz) + nx*nz*(-ox+px) + ny*nz*(-oy+py)+nz*nz*(m*nz*nz*oz+pz))*sin*sin
				+ m*cos*(ny*nz*oy - m*ny*nz*oy - oz - m*ny*ny*oz + nz*nz*oz - 2*m*nz*nz*oz - (-1+m)*nx*nz*(ox-py)
				- ny*nz*py + m*ny*nz*py - 2*pz + ny*ny*pz - m*ny*ny*pz + nx*nx*(pz-m*(oz+pz)) + m*(-1+nx*nx+ny*ny+nz*nz)
				*(ny*(-ox+px)+nx*(oy-py))*sin + m*m*(-nx*nx + nx*nx*nx*nx - ny*ny + ny*ny*ny*ny - nz*nz + nz*nz*nz*nz)*oz*sin*sin)
				+ m*sin*((1+m*(nx*nx+ny*ny+nz*nz))*(ny*(-ox+px)+nx*(oy-py)) + m*m*(ny*ny*nz*nz + nx*nx*(ny*ny+nz*nz))*oz*s2a));
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

	@Override
	public void interpolate(float t, Frame toSet) {
		interpolate(t, start.P, toSet.P, true);
		Vector3f refP = new Vector3f();
		interpolate(t, start.P, refP, false);

		Vector3f tmp = start.P.add(start.I);
		interpolate(t, tmp, toSet.I, false);
		toSet.I.subtractLocal(refP);

		tmp = start.P.add(start.J);
		interpolate(t, tmp, toSet.J, false);
		toSet.J.subtractLocal(refP);
		
		tmp = start.P.add(start.K);
		interpolate(t, tmp, toSet.K, false);
		toSet.K.subtractLocal(refP);

//		System.out.println("Frame at t="+t+": "+toSet);
	}

	private void interpolate(float t, Vector3f P0, Vector3f store, boolean lerpZ) {
		store.set(F);
		Vector3f FP0 = P0.subtract(F);
		Vector3f rot = rotate2(alpha, t, FP0, lerpZ);
//		Vector3f rot = rotate(alpha*t, FP0);
		store.addScaleLocal(rot, (float) Math.pow(m, t));
	}
	
	private Vector3f lerp(Vector3f A, Vector3f B, float t) {
		Vector3f result = new Vector3f();
		result.addScaleLocal(A, 1-t);
		result.addScaleLocal(B, t);
		return result;
	}

	private Vector3f rotate(float angle, Vector3f X) {
		Vector3f W = project(X, N);
		Vector3f U = X.subtract(W);
		Vector3f result = new Vector3f(W);
		result.addScaleLocal(U, (float) Math.cos(angle));
		result.addScaleLocal(N.cross(U), -(float) Math.sin(angle));
		return result;
	}
	
	private Vector3f rotate2(float angle, float t, Vector3f X, boolean lerpZ) {
		Vector3f W = project(X, N);
		Vector3f U = X.subtract(W);
		Vector3f result = lerp(W, W.negate(), lerpZ ? t : 0);
		result.addScaleLocal(U, (float) Math.cos(angle*t));
		result.addScaleLocal(N.cross(U), -(float) Math.sin(angle*t));
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
