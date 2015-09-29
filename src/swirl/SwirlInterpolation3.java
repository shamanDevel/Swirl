/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package swirl;

/**
 *
 * @author Sebastian Weiss
 */
public class SwirlInterpolation3 implements Interpolation {
	private static final double EPSILON = 0.00001;
	
	private Frame start;
	private Frame end;
	
	private Vector3f N = new Vector3f();
	private float spiralAngle;
	private float spiralScale;
	private Vector2f spiralCenter;
	private Vector3f projI;
	private Vector3f projJ;
	private Vector2f A;
	private float zTranslation;
	private float translationScale;
	private Quaternion startQuat = new Quaternion();
	private Quaternion endQuat = new Quaternion();
	
	private Pv3D.pt testA;
	private Pv3D.vec testAB;
	private Pv3D.pt testC;
	private Pv3D.vec testCD;
	private Pv3D.pt testF;

	@Override
	public void setStartEnd(Frame start, Frame end) {
		this.start = start;
		this.end = end;
		System.out.println("Calculate Swirl 3");
		System.out.println(" Start frame: "+start);
		System.out.println(" End frame: "+end);
		
		//Calculate N and alpha
		Vector3f dI = new Vector3f(end.I.normalize());
		 dI.subtractLocal(start.I.normalize());
		Vector3f dJ = new Vector3f(end.J.normalize()); 
		 dJ.subtractLocal(start.J.normalize());
		Vector3f dK = new Vector3f(end.K.normalize());
		 dK.subtractLocal(start.K.normalize());
		Vector3f dir1, dir2;
		dI.cross(dI, N);
		if (N.lengthSquared() < EPSILON) {
			dI.cross(dK, N);
			if (N.lengthSquared() < EPSILON) {
				dJ.cross(dK, N);
				if (N.lengthSquared() < EPSILON) {
					System.out.println("Special case: no rotation");
					//TODO: implement
					return;
				} else {
					dir1 = start.J; dir2 = end.J;
				}
			} else {
				dir1 = start.I; dir2 = end.I;
			}
		} else {
			dir1 = start.I; dir2 = end.I;
		}
		N.normalizeLocal();
		System.out.println(" rotation axis N="+N);
		System.out.println(" dir1="+dir1+"  dir2="+dir2);
		
		//compute the swirl in 2d
		projI = projectOnPlane(dir1, N);
		projI.normalizeLocal();
		projJ = N.cross(projI);
		System.out.println(" projI="+projI+"  projJ="+projJ);
		A = new Vector2f(start.P.dot(projI), start.P.dot(projJ));
		float projILength = projI.length();
		Vector2f B = new Vector2f(A.x + projILength, A.y);
		Vector2f C = new Vector2f(end.P.dot(projI), end.P.dot(projJ));
		Vector3f tmp = end.P.add(dir2);
		Vector2f D = new Vector2f(tmp.dot(projI), tmp.dot(projJ));
		System.out.println(" A="+A+"  B="+B+"  C="+C+"  D="+D);
		Vector2f AB = new Vector2f(B); AB.subtractLocal(A);
		Vector2f CD = new Vector2f(D); CD.subtractLocal(C);
		System.out.println(" AB="+AB+"  CD="+CD);
		spiralAngle = spiralAngle(AB, CD);
		spiralScale = spiralScale(AB, CD);
		spiralCenter = spiralCenter(spiralAngle, spiralScale, A, C);
		System.out.println(" spiral angle alpha="+spiralAngle+"  ( "+(spiralAngle*180/Math.PI)+"° )");
		System.out.println(" spiral scale m="+spiralScale);
		System.out.println(" spiral center F="+spiralCenter);
		
		zTranslation = end.P.dot(N) - start.P.dot(N);
		System.out.println(" z-translation="+zTranslation);
		translationScale = end.I.length() / start.I.length();
		System.out.println(" translation scale="+translationScale);
		
		startQuat.fromAxes(start.I.normalize(), start.J.normalize(), start.K.normalize());
		endQuat.fromAxes(end.I.normalize(), end.J.normalize(), end.K.normalize());
		
		testA = Pv3D.P(A.x, Pv3D.P(projI.x, projI.y, projI.z), A.y, Pv3D.P(projJ.x, projJ.y, projJ.z));
		Pv3D.pt testB = Pv3D.P(B.x, Pv3D.P(projI.x, projI.y, projI.z), B.y, Pv3D.P(projJ.x, projJ.y, projJ.z));
		testAB = Pv3D.V(testA, testB);
		testC = Pv3D.P(C.x, Pv3D.P(projI.x, projI.y, projI.z), C.y, Pv3D.P(projJ.x, projJ.y, projJ.z));
		Pv3D.pt testD = Pv3D.P(D.x, Pv3D.P(projI.x, projI.y, projI.z), D.y, Pv3D.P(projJ.x, projJ.y, projJ.z));
		testCD = Pv3D.V(testC, testD);
		testF = Pv3D.P(spiralCenter.x, Pv3D.P(projI.x, projI.y, projI.z), spiralCenter.y, Pv3D.P(projJ.x, projJ.y, projJ.z));
		
		//System.out.println(" scaling m="+m);
		//System.out.println(" rotation angle alpha="+alpha+"  ( "+(alpha*180/Math.PI)+"° )");
	}
	
	private Vector3f project(Vector3f X, Vector3f N) {
		return N.mult(X.dot(N));
	}
	
	private Vector3f projectOnPlane(Vector3f X, Vector3f N) {
		Vector3f v = N.mult(X.dot(N));
		v.set(X.x - v.x, X.y - v.y, X.z - v.z);
		return v;
	}
	
	private void rotate(float angle, float scale, Vector2f Q, Vector2f origin, Vector2f store) {
		float dx = Q.x - origin.x;
		float dy = Q.y - origin.y;
		float c = (float) Math.cos(angle);
		float s = (float) Math.sin(angle);
		store.set(origin.x + scale * (c * dx - s * dy), origin.y + scale * (s * dx + c * dy));
	}
	
	private Vector2f spiralCenter(Vector2f A, Vector2f B, Vector2f C, Vector2f D) { // computes center of spiral that takes A to C and B to D
		Vector2f AB = new Vector2f(B); AB.subtractLocal(A);
		Vector2f CD = new Vector2f(D); AB.subtractLocal(C);
		float a = spiralAngle(AB, CD);
		float z = spiralScale(AB, CD);
		return spiralCenter(a, z, A, C);
	}

	private float spiralAngle(Vector2f AB, Vector2f CD) {
		return AB.angleBetween(CD);
	}

	private float spiralScale(Vector2f AB, Vector2f CD) {
		return CD.length() / AB.length();
	}

	private static final float sq(float x) {return x*x;}
	private Vector2f spiralCenter(float a, float z, Vector2f A, Vector2f C) {
		float c = (float) Math.cos(a), s = (float) Math.sin(a);
		float D = sq(c * z - 1) + sq(s * z);
		float ex = c * z * A.x - C.x - s * z * A.y;
		float ey = c * z * A.y - C.y + s * z * A.x;
		float x = (ex * (c * z - 1) + ey * s * z) / D;
		float y = (ey * (c * z - 1) - ex * s * z) / D;
		return new Vector2f(x, y);
	}
	
	private Vector3f project(float x, float y, float z,
			Vector3f basisI, Vector3f basisJ,  Vector3f basisK) {
		return new Vector3f(
				x*basisI.x + y*basisJ.x + z*basisK.x, 
				x*basisI.y + y*basisJ.y + z*basisK.y, 
				x*basisI.z + y*basisJ.z + z*basisK.z);
	}
	
	@Override
	public void interpolate(float t, Frame toSet) {
		//Interpolate P
		toSet.P.set(interpolate(t, new Vector3f(start.P)));
		
		//Interpolate I,J,K
		Quaternion q = new Quaternion();
		q.slerp(startQuat, endQuat, t);
		q.toAxes(new Vector3f[]{toSet.I, toSet.J, toSet.K});
		float s = (float) (start.I.length() * Math.pow(translationScale, t));
		toSet.I.multLocal(s);
		toSet.J.multLocal(s);
		toSet.K.multLocal(s);
		
//		Vector3f PI = new Vector3f(start.P); PI.addLocal(start.I);
//		PI = interpolate(t, PI);
//		PI.subtractLocal(toSet.P);
//		PI.multLocal((float) Math.pow(translationScale, t));
//		toSet.I.set(PI);
//		Vector3f PJ = new Vector3f(start.P); PJ.addLocal(start.J);
//		PJ = interpolate(t, PJ);
//		PJ.subtractLocal(toSet.P);
//		PJ.multLocal((float) Math.pow(translationScale, t));
//		toSet.J.set(PJ);
//		Vector3f PK = new Vector3f(start.P); PK.addLocal(start.K);
//		PK = interpolate(t, PK);
//		PK.subtractLocal(toSet.P);
//		PK.multLocal((float) Math.pow(translationScale, t));
//		toSet.K.set(PK);
	}
	private Vector3f interpolate(float t, Vector3f point) {
		//Project in rotation plane
		float x = projI.dot(point);
		float y = projJ.dot(point);
		float z = N.dot(point);
		//Rotate in 2D
		Vector2f P = new Vector2f();
		float angle = spiralAngle * t;
		rotate(angle, (float) Math.pow(spiralScale, t), new Vector2f(x, y), spiralCenter, P);
		//project back
		Vector3f projected = project(P.x, P.y, z, projI, projJ, N);
		//move point
		Vector3f translation = new Vector3f(N);
		double dz = z - N.dot(new Vector3f(start.P));
		translation.multLocal((float) (t * zTranslation + (Math.pow(translationScale, t) - 1) * dz)); //unsteady
//		translation.multLocal(t * zTranslation); //unsteady
		projected.addLocal(translation);
		
		return projected;
	}
	
	@Override
	public void debugDraw(Swirl applet) {
		applet.noStroke();
		applet.fill(Swirl.black);
		applet.pushMatrix();
		 applet.translate(testA.x * Swirl.SCALE, testA.y * Swirl.SCALE, testA.z * Swirl.SCALE);
		 applet.sphere(0.1f*Swirl.SCALE);
		applet.popMatrix();
		applet.arrow(Pv3D.P(testA.x * Swirl.SCALE, testA.y * Swirl.SCALE, testA.z * Swirl.SCALE),
				Pv3D.V(Swirl.SCALE, testAB), 0.1f*Swirl.SCALE);
		applet.pushMatrix();
		 applet.translate(testC.x * Swirl.SCALE, testC.y * Swirl.SCALE, testC.z * Swirl.SCALE);
		 applet.sphere(0.1f*Swirl.SCALE);
		applet.popMatrix();
		applet.arrow(Pv3D.P(testC.x * Swirl.SCALE, testC.y * Swirl.SCALE, testC.z * Swirl.SCALE),
				Pv3D.V(Swirl.SCALE, testCD), 0.1f*Swirl.SCALE);
		
		applet.fill(Swirl.grey);
		applet.pushMatrix();
		 applet.translate(testF.x * Swirl.SCALE, testF.y * Swirl.SCALE, testF.z * Swirl.SCALE);
		 applet.sphere(0.1f*Swirl.SCALE);
		 applet.collar(Pv3D.P(), Pv3D.V(N.x * 50 * Swirl.SCALE, N.y * 50 * Swirl.SCALE, N.z * 50 * Swirl.SCALE), 
				 0.05f*Swirl.SCALE, 0.05f*Swirl.SCALE);
		 applet.collar(Pv3D.P(), Pv3D.V(N.x *-50 * Swirl.SCALE, (N.y *-50 *Swirl. SCALE), (N.z *-50 * Swirl.SCALE)), 
				 0.05f*Swirl.SCALE, 0.05f*Swirl.SCALE);
		applet.popMatrix();
	}

	@Override
	public String debugString() {
		StringBuilder str = new StringBuilder();
		str.append("Modified Swirl Interpolation: 2d log-spiral in the plane orthogonal to N + exponential translation along N\n");
		str.append("rotation axis N=").append(N).append("\n");
		str.append("spiral angle=").append(spiralAngle).append(" (").append(spiralAngle * 180 / Math.PI)
				.append("°), spiral scale=").append(spiralScale).append(", spiral center=").append(spiralCenter).append("\n");
		str.append("translation along N=").append(zTranslation).append(", translation scale=").append(translationScale);
		return str.toString();
	}
	
}
