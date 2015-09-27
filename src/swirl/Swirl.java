/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package swirl;

import processing.core.*;
import processing.event.*;

import static swirl.Pv3D.*;

/**
 *
 * @author Sebastian Weiss
 */
public class Swirl extends PApplet {
	public static final float SCALE = 50;
	private static final float PICK_TOLERANCE = 20;
	public static final float DEG_TO_RAD = PI / 180.0f;

	float dz = 0; // distance to camera. Manipulated with wheel or when 
//float rx=-0.06*TWO_PI, ry=-0.04*TWO_PI;    // view angles manipulated when space pressed but not mouse
	float rx = 0, ry = 0;    // view angles manipulated when space pressed but not mouse
	Boolean twistFree = false, animating = true, tracking = false, center = true, gouraud = true, showControlPolygon = false, showNormals = false;
	float t = 0, s = 0;
	boolean viewpoint = false;
	PImage myFace;
	boolean filming = false;
	boolean change = false;
	int frameCounter = 0;
	
	private final Frame startFrame = new Frame();
	private final Frame endFrame = new Frame();
	private final int interpolatingFramesCount = 10;
	private final Frame[] interpolatingFrames = new Frame[interpolatingFramesCount];
	private final int extrapolatingFramesCount = 50;
	private final float extrapolatingFramesLength = 5;
	private final Frame[] extrapolatingFrames = new Frame[extrapolatingFramesCount];
	private boolean recalculateFrames = false;
	private int selectedInterpolation;
	private Interpolation[] interpolations;
	
	private int pickedPoint = 0;

	String title = "6491 P2 2015: 3D swirl", name = "Sebastian Weiß",
			menu = "1-4: change interpolation, !:picture, ~:(start/stop)capture, space:rotate, s/wheel:closer, a:anim, #:quit",
			guide = "click'n'drag center of frames or arrow tips to change the start frame (green) and end frame (red)"; // user's guide

	public void settings() {
		size(900, 900, P3D); // p3D means that we will do 3D graphics
		noSmooth();
	}

	public void setup() {
		myFace = loadImage("data/pic.jpg");  // load image from file pic.jpg in folder data *** replace that file with your pic of your own face
		textureMode(NORMAL);
		
		startFrame.P.set(-4, -1, -1);
		startFrame.I.set(1, 0, 0);
		startFrame.J.set(0, 1, 0);
		startFrame.K.set(0, 0, 1);
		Quaternion quat = new Quaternion();
		quat.fromAngles(DEG_TO_RAD * 20, DEG_TO_RAD * -50, DEG_TO_RAD * 120);
		quat.multLocal(startFrame.I);
		quat.multLocal(startFrame.J);
		quat.multLocal(startFrame.K);
		
		endFrame.P.set(4, 1, 1);
		endFrame.I.set(1.1f, 0, 0);
		endFrame.J.set(0, 1.1f, 0);
		endFrame.K.set(0, 0, 1.1f);
		quat.fromAngles(DEG_TO_RAD * -60, DEG_TO_RAD * 15, DEG_TO_RAD * 30);
		quat.multLocal(endFrame.I);
		quat.multLocal(endFrame.J);
		quat.multLocal(endFrame.K);
		
		System.out.println("Start frame: "+startFrame);
		System.out.println("End frame: "+endFrame);
		
		for (int i=0; i<interpolatingFramesCount; ++i) {
			interpolatingFrames[i] = new Frame();
		}
		for (int i=0; i<extrapolatingFramesCount; ++i) {
			extrapolatingFrames[i] = new Frame();
		}
		
		interpolations = new Interpolation[]{
			new LinearInterpolation(),
			new SwirlInterpolation(),
			new SwirlInterpolation2(),
			new SAMInterpolation()
		};
		selectedInterpolation = 1;
		recalculateFrames = true;
		
		rx+=0.0001f;
		ry+=0.0001f;
	}

	public void draw() {
		background(255);
		pushMatrix();   // to ensure that we can restore the standard view before writing on the canvas

		float fov = PI / 3.0f;
		float cameraZ = (height / 2.0f) / tan(fov / 2.0f);
		camera(0, 0, cameraZ, 0, 0, 0, 0, 1, 0);       // sets a standard perspective
		perspective(fov, 1.0f, 0.1f, 10000);

		translate(0, 0, dz); // puts origin of model at screen center and moves forward/away by dz
		lights();  // turns on view-dependent lighting
		rotateX(rx);
		rotateY(ry); // rotates the model around the new origin (center of screen)
		rotateX(PI / 2); // rotates frame around X to make X and Y basis vectors parallel to the floor
		noStroke(); // if you use stroke, the weight (width) of it will be scaled with you scaleing factor
		showFrame(50); // X-red, Y-green, Z-blue arrows

		computeProjectedVectors(); // computes screen projections I, J, K of basis vectors (see bottom of pv3D): used for dragging in viewer's frame    
		
		calculateAndShowFrames();
		interpolations[selectedInterpolation].debugDraw(this);
		
		popMatrix(); // done with 3D drawing. Restore front view for writing text on canvas

      // for demos: shows the mouse and the key pressed (but it may be hidden by the 3D model)
		//  if(keyPressed) {stroke(red); fill(white); ellipse(mouseX,mouseY,26,26); fill(red); text(key,mouseX-5,mouseY+4);}
		fill(0xff000000);
		displayHeader(); // dispalys header on canvas, including my face
		if (!filming) {
			displayFooter(); // shows menu at bottom, only if not filming
		}
		if (animating) {
			t += PI / 180 / 2;
			if (t >= TWO_PI) {
				t = 0;
			}
			s = (float) ((cos(t) + 1.) / 2);
		} // periodic change of time 
		if (filming && (animating || change)) {
			saveFrame("FRAMES/F" + nf(frameCounter++, 4) + ".tif");  // save next frame to make a movie
		}
		change = false; // to avoid capturing frames when nothing happens (change is set uppn action)
	}
	
	private void showFrame(Frame frame, int color) {
		noStroke();
		fill(color);
		Pv3D.pt P = new Pv3D.pt(frame.P.x*SCALE, frame.P.y*SCALE, frame.P.z*SCALE);
		Pv3D.vec I = V(frame.I.x*SCALE, frame.I.y*SCALE, frame.I.z*SCALE);
		Pv3D.vec J = V(frame.J.x*SCALE, frame.J.y*SCALE, frame.J.z*SCALE);
		Pv3D.vec K = V(frame.K.x*SCALE, frame.K.y*SCALE, frame.K.z*SCALE);
		pushMatrix();
		translate(P.x, P.y, P.z);
		sphere(0.1f*SCALE);
		popMatrix();
		arrow(P, I, 0.1f*SCALE);
		arrow(P, J, 0.1f*SCALE);
		arrow(P, K, 0.1f*SCALE);
	}
	
	private void calculateAndShowFrames() {
		if (recalculateFrames) {
			recalculateFrames = false;
			Interpolation interpolation = interpolations[selectedInterpolation];
			float interpolationStep = 1.0f / (interpolatingFramesCount+1);
//			float interpolationStep = 1.0f / (interpolatingFramesCount);
			float extrapolationStep = extrapolatingFramesLength / extrapolatingFramesCount;
			interpolation.setStartEnd(startFrame, endFrame);
			for (int i=0; i<interpolatingFramesCount; ++i) {
				interpolation.interpolate((i+1)*interpolationStep, interpolatingFrames[i]);
			}
			for (int i=0; i<extrapolatingFramesCount; ++i) {
				interpolation.interpolate((i+1)*extrapolationStep + 1, extrapolatingFrames[i]);
			}
		}
		showFrame(startFrame, green);
		showFrame(endFrame, red);
		for (int i=0; i<interpolatingFramesCount; ++i) {
			showFrame(interpolatingFrames[i], blue);
		}
		for (int i=0; i<extrapolatingFramesCount; ++i) {
			showFrame(extrapolatingFrames[i], blue50);
		}
	}

	public void keyPressed() {
		if (key == '1') {
			selectedInterpolation = 0; recalculateFrames = true;
		} else if (key == '2' && interpolations.length >= 2) {
			selectedInterpolation = 1; recalculateFrames = true;
		} else if (key == '3' && interpolations.length >= 3) {
			selectedInterpolation = 2; recalculateFrames = true;
		} else if (key == '4' && interpolations.length >= 4) {
			selectedInterpolation = 3; recalculateFrames = true;
		} else if (key == '5' && interpolations.length >= 5) {
			selectedInterpolation = 4; recalculateFrames = true; //unused
		}
		
		if (key == '?') {
			scribeText = !scribeText;
		}
		if (key == '~') {
			filming = !filming;
		}
		if (key == ']') {
			showControlPolygon = !showControlPolygon;
		}
		if (key == '|') {
			showNormals = !showNormals;
		}
		if (key == 'G') {
			gouraud = !gouraud;
		}

		// if(key=='.') F=P.Picked(); // snaps focus F to the selected vertex of P (easier to rotate and zoom while keeping it in center)
		if (key == 'c') {
			center = !center; // snaps focus F to the selected vertex of P (easier to rotate and zoom while keeping it in center)
		}
		if (key == 't') {
			tracking = !tracking; // snaps focus F to the selected vertex of P (easier to rotate and zoom while keeping it in center)
		}
		if (key == 'a') {
			animating = !animating; // toggle animation
		}
		if (key == ',') {
			viewpoint = !viewpoint;
		}
		if (key == '#') {
			exit();
		}
		change = true;
	}

	public void mouseWheel(MouseEvent event) {
		dz += event.getAmount();
		change = true;
	}

	public void mousePressed() {
		pushMatrix();
		float fov = PI / 3.0f;
		float cameraZ = (height / 2.0f) / tan(fov / 2.0f);
		camera(0, 0, cameraZ, 0, 0, 0, 0, 1, 0);       // sets a standard perspective
		perspective(fov, 1.0f, 0.1f, 10000);
		translate(0, 0, dz); // puts origin of model at screen center and moves forward/away by dz
		rotateX(rx);
		rotateY(ry); // rotates the model around the new origin (center of screen)
		rotateX(PI / 2); // rotates frame around X to make X and Y basis vectors parallel to the floor
		
		float mx = mouseX;
		float my = mouseY;
//		System.out.println("Mouse X="+mx+" Y="+my);
		float dist = Float.MAX_VALUE;
		pickedPoint = 0;
		for (int i=0; i<4; ++i) {
			Vector3f p = startFrame.getPickablePoints()[i].mult(SCALE);
			float x = screenX(p.x, p.y, p.z);
			float y = screenY(p.x, p.y, p.z);
			float z = screenZ(p.x, p.y, p.z);
//			System.out.println("Point "+(i+1)+" X="+x+" Y="+y+" Z="+z);
			if (Math.abs(x-mx) > PICK_TOLERANCE || Math.abs(y-my) > PICK_TOLERANCE)
				continue;
			if (z < dist) {
				dist = z;
				pickedPoint = i + 1;
			}
		}
		for (int i=0; i<4; ++i) {
			Vector3f p = endFrame.getPickablePoints()[i].mult(SCALE);
			float x = screenX(p.x, p.y, p.z);
			float y = screenY(p.x, p.y, p.z);
			float z = screenZ(p.x, p.y, p.z);
//			System.out.println("Point "+(i+5)+" X="+x+" Y="+y+" Z="+z);
			if (Math.abs(x-mx) > PICK_TOLERANCE || Math.abs(y-my) > PICK_TOLERANCE)
				continue;
			if (z < dist) {
				dist = z;
				pickedPoint = i + 5;
			}
		}
		System.out.println("Picked: "+pickedPoint);
		
		popMatrix();
	}

	@Override
	public void mouseReleased() {
		pickedPoint = 0;
	}

	public void mouseMoved() {
		if (keyPressed && key == ' ') {
			rx -= PI * (mouseY - pmouseY) / height;
			ry += PI * (mouseX - pmouseX) / width;
		};
		if (keyPressed && key == 's') {
			dz += (float) (mouseY - pmouseY); // approach view (same as wheel)
		}
	}

	public void mouseDragged() {
		if (pickedPoint == 0) {
			return;
		}
		Vector3f movement = new Vector3f();
		Pv3D.vec v = ToIJ(V((mouseX-pmouseX),0,0));
		movement.addLocal(v.x, v.y, v.z);
		v = ToK(V(0, (mouseY-pmouseY),0));
		movement.addLocal(v.x, v.y, v.z);
		movement.divideLocal(SCALE);
		if (pickedPoint >= 1 && pickedPoint <= 4) {
			startFrame.movePickablePoint(pickedPoint - 1, movement);
		} else if (pickedPoint >= 5) {
			endFrame.movePickablePoint(pickedPoint - 5, movement);
		}
		recalculateFrames = true;
	}

// **** Header, footer, help text on canvas
	private void displayHeader() { // Displays title and authors face on screen
		scribeHeader(title, 0);
		scribeHeader(interpolations[selectedInterpolation].debugString(), 1);
		scribeHeaderRight(name);
		fill(white);
		image(myFace, width - myFace.width / 2, 25, myFace.width / 2, myFace.height / 2);
	}

	private void displayFooter() { // Displays help text at the bottom
		scribeFooter(guide, 1);
		scribeFooter(menu, 0);
	}

	Pv3D.pt ToScreen(Pv3D.pt P) {
		return new Pv3D.pt(screenX(P.x, P.y, P.z), screenY(P.x, P.y, P.z), 0);
	}  // O+xI+yJ+kZ

	Pv3D.pt ToModel(Pv3D.pt P) {
		return new Pv3D.pt(modelX(P.x, P.y, P.z), modelY(P.x, P.y, P.z), modelZ(P.x, P.y, P.z));
	}  // O+xI+yJ+kZ

// ===== mouse
	Pv3D.pt Mouse() {
		return new Pv3D.pt(mouseX, mouseY, 0);
	}

	;                                          // current mouse location
Pv3D.pt Pmouse() {
		return new Pv3D.pt(pmouseX, pmouseY, 0);
	}

	;
	Pv3D.vec MouseDrag() {
		return new Pv3D.vec(mouseX - pmouseX, mouseY - pmouseY, 0);
	}

	;                     // vector representing recent mouse displacement
Pv3D.pt ScreenCenter() {
		return new Pv3D.pt(width / 2, height / 2);
	}                                                        //  point in center of  canvas

// ===== render
	void normal(Pv3D.vec V) {
		normal(V.x, V.y, V.z);
	}

	;                                          // changes normal for smooth shading
void vertex(Pv3D.pt P) {
		vertex(P.x, P.y, P.z);
	}

	;                                           // vertex for shading or drawing
void v(Pv3D.pt P) {
		vertex(P.x, P.y, P.z);
	}

	;                                           // vertex for shading or drawing
void nv(Pv3D.vec N) {
		normal(N.x, N.y, N.z);
	}

	;                                           // vertex for shading or drawing
void vTextured(Pv3D.pt P, float u, float v) {
		vertex(P.x, P.y, P.z, u, v);
	}

	;                          // vertex with texture coordinates
void show(Pv3D.pt P, Pv3D.pt Q) {
		line(Q.x, Q.y, Q.z, P.x, P.y, P.z);
	}

	;                       // draws edge (P,Q)
void show(Pv3D.pt P, Pv3D.vec V) {
		line(P.x, P.y, P.z, P.x + V.x, P.y + V.y, P.z + V.z);
	}

	;          // shows edge from P to P+V
void show(Pv3D.pt P, float d, Pv3D.vec V) {
		line(P.x, P.y, P.z, P.x + d * V.x, P.y + d * V.y, P.z + d * V.z);
	}

	; // shows edge from P to P+dV
void show(Pv3D.pt A, Pv3D.pt B, Pv3D.pt C) {
		beginShape();
		vertex(A);
		vertex(B);
		vertex(C);
		endShape(CLOSE);
	}

	;                      // volume of tet 
void show(Pv3D.pt A, Pv3D.pt B, Pv3D.pt C, Pv3D.pt D) {
		beginShape();
		vertex(A);
		vertex(B);
		vertex(C);
		vertex(D);
		endShape(CLOSE);
	}

	;                      // volume of tet 
void show(Pv3D.pt P, float r) {
		pushMatrix();
		translate(P.x, P.y, P.z);
		sphere(r);
		popMatrix();
	}

	; // render sphere of radius r and center P
void show(Pv3D.pt P, float s, Pv3D.vec I, Pv3D.vec J, Pv3D.vec K) {
		noStroke();
		fill(yellow);
		show(P, 5);
		stroke(red);
		show(P, s, I);
		stroke(green);
		show(P, s, J);
		stroke(blue);
		show(P, s, K);
	}

	; // render sphere of radius r and center P
void show(Pv3D.pt P, String s) {
		text(s, P.x, P.y, P.z);
	}

	; // prints string s in 3D at P
void show(Pv3D.pt P, String s, Pv3D.vec D) {
		text(s, P.x + D.x, P.y + D.y, P.z + D.z);
	}

	; // prints string s in 3D at P+D
void showShadow(Pv3D.pt P, float r) {
		pushMatrix();
		translate(P.x, P.y, 0);
		scale(1, 1, 0.01f);
		sphere(r);
		popMatrix();
	}

	String toText(Pv3D.vec V) {
		return "(" + nf(V.x, 1, 5) + "," + nf(V.y, 1, 5) + "," + nf(V.z, 1, 5) + ")";
	}
// ==== curve

	void bezier(Pv3D.pt A, Pv3D.pt B, Pv3D.pt C, Pv3D.pt D) {
		bezier(A.x, A.y, A.z, B.x, B.y, B.z, C.x, C.y, C.z, D.x, D.y, D.z);
	} // draws a cubic Bezier curve with control points A, B, C, D

	void bezier(Pv3D.pt[] C) {
		bezier(C[0], C[1], C[2], C[3]);
	} // draws a cubic Bezier curve with control points A, B, C, D

	Pv3D.pt bezierPoint(Pv3D.pt[] C, float t) {
		return P(bezierPoint(C[0].x, C[1].x, C[2].x, C[3].x, t), bezierPoint(C[0].y, C[1].y, C[2].y, C[3].y, t), bezierPoint(C[0].z, C[1].z, C[2].z, C[3].z, t));
	}

	Pv3D.vec bezierTangent(Pv3D.pt[] C, float t) {
		return V(bezierTangent(C[0].x, C[1].x, C[2].x, C[3].x, t), bezierTangent(C[0].y, C[1].y, C[2].y, C[3].y, t), bezierTangent(C[0].z, C[1].z, C[2].z, C[3].z, t));
	}

	void PT(Pv3D.pt P0, Pv3D.vec T0, Pv3D.pt P1, Pv3D.vec T1) {
		float d = d(P0, P1) / 3;
		bezier(P0, P(P0, -d, U(T0)), P(P1, -d, U(T1)), P1);
	} // draws cubic Bezier interpolating  (P0,T0) and  (P1,T1) 

	void PTtoBezier(Pv3D.pt P0, Pv3D.vec T0, Pv3D.pt P1, Pv3D.vec T1, Pv3D.pt[] C) {
		float d = d(P0, P1) / 3;
		C[0].set(P0);
		C[1].set(P(P0, -d, U(T0)));
		C[2].set(P(P1, -d, U(T1)));
		C[3].set(P1);
	} // draws cubic Bezier interpolating  (P0,T0) and  (P1,T1) 

	Pv3D.vec vecToCubic(Pv3D.pt A, Pv3D.pt B, Pv3D.pt C, Pv3D.pt D, Pv3D.pt E) {
		return V((-A.x + 4 * B.x - 6 * C.x + 4 * D.x - E.x) / 6, (-A.y + 4 * B.y - 6 * C.y + 4 * D.y - E.y) / 6, (-A.z + 4 * B.z - 6 * C.z + 4 * D.z - E.z) / 6);
	}

	Pv3D.vec vecToProp(Pv3D.pt B, Pv3D.pt C, Pv3D.pt D) {
		float cb = d(C, B);
		float cd = d(C, D);
		return V(C, P(B, cb / (cb + cd), D));
	}

	;  

// ==== perspective
Pv3D.pt Pers(Pv3D.pt P, float d) {
		return P(d * P.x / (d + P.z), d * P.y / (d + P.z), d * P.z / (d + P.z));
	}

	;
	Pv3D.pt InverserPers(Pv3D.pt P, float d) {
		return P(d * P.x / (d - P.z), d * P.y / (d - P.z), d * P.z / (d - P.z));
	}

	;

// ==== intersection
boolean intersect(Pv3D.pt P, Pv3D.pt Q, Pv3D.pt A, Pv3D.pt B, Pv3D.pt C, Pv3D.pt X) {
		return intersect(P, V(P, Q), A, B, C, X);
	} // if (P,Q) intersects (A,B,C), return true and set X to the intersection point

	boolean intersect(Pv3D.pt E, Pv3D.vec T, Pv3D.pt A, Pv3D.pt B, Pv3D.pt C, Pv3D.pt X) { // if ray from E along T intersects triangle (A,B,C), return true and set X to the intersection point
		Pv3D.vec EA = V(E, A), EB = V(E, B), EC = V(E, C), AB = V(A, B), AC = V(A, C);
		boolean s = cw(EA, EB, EC), sA = cw(T, EB, EC), sB = cw(EA, T, EC), sC = cw(EA, EB, T);
		if ((s == sA) && (s == sB) && (s == sC)) {
			return false;
		}
		float t = m(EA, AC, AB) / m(T, AC, AB);
		X.set(P(E, t, T));
		return true;
	}

	boolean rayIntersectsTriangle(Pv3D.pt E, Pv3D.vec T, Pv3D.pt A, Pv3D.pt B, Pv3D.pt C) { // true if ray from E with direction T hits triangle (A,B,C)
		Pv3D.vec EA = V(E, A), EB = V(E, B), EC = V(E, C);
		boolean s = cw(EA, EB, EC), sA = cw(T, EB, EC), sB = cw(EA, T, EC), sC = cw(EA, EB, T);
		return (s == sA) && (s == sB) && (s == sC);
	}

	;
	boolean edgeIntersectsTriangle(Pv3D.pt P, Pv3D.pt Q, Pv3D.pt A, Pv3D.pt B, Pv3D.pt C) {
		Pv3D.vec PA = V(P, A), PQ = V(P, Q), PB = V(P, B), PC = V(P, C), QA = V(Q, A), QB = V(Q, B), QC = V(Q, C);
		boolean p = cw(PA, PB, PC), q = cw(QA, QB, QC), a = cw(PQ, PB, PC), b = cw(PA, PQ, PC), c = cw(PQ, PB, PQ);
		return (p != q) && (p == a) && (p == b) && (p == c);
	}

	float rayParameterToIntersection(Pv3D.pt E, Pv3D.vec T, Pv3D.pt A, Pv3D.pt B, Pv3D.pt C) {
		Pv3D.vec AE = V(A, E), AB = V(A, B), AC = V(A, C);
		return -m(AE, AC, AB) / m(T, AC, AB);
	}

	float angleDraggedAround(Pv3D.pt G) {  // returns angle in 2D dragged by the mouse around the screen projection of G
		Pv3D.pt S = P(screenX(G.x, G.y, G.z), screenY(G.x, G.y, G.z), 0);
		Pv3D.vec T = V(S, Pmouse());
		Pv3D.vec U = V(S, Mouse());
		return atan2(d(R(U), T), d(U, T));
	}

	float scaleDraggedFrom(Pv3D.pt G) {
		Pv3D.pt S = P(screenX(G.x, G.y, G.z), screenY(G.x, G.y, G.z), 0);
		return d(S, Mouse()) / d(S, Pmouse());
	}

// FANS, CONES, AND ARROWS
	void disk(Pv3D.pt P, Pv3D.vec V, float r) {
		Pv3D.vec I = U(Normal(V));
		Pv3D.vec J = U(N(I, V));
		disk(P, I, J, r);
	}

	void disk(Pv3D.pt P, Pv3D.vec I, Pv3D.vec J, float r) {
		float da = TWO_PI / 36;
		beginShape(TRIANGLE_FAN);
		v(P);
		for (float a = 0; a <= TWO_PI + da; a += da) {
			v(P(P, r * cos(a), I, r * sin(a), J));
		}
		endShape();
	}

	void fan(Pv3D.pt P, Pv3D.vec V, float r) {
		Pv3D.vec I = U(Normal(V));
		Pv3D.vec J = U(N(I, V));
		fan(P, V, I, J, r);
	}

	void fan(Pv3D.pt P, Pv3D.vec V, Pv3D.vec I, Pv3D.vec J, float r) {
		float da = TWO_PI / 36;
		beginShape(TRIANGLE_FAN);
		v(P(P, V));
		for (float a = 0; a <= TWO_PI + da; a += da) {
			v(P(P, r * cos(a), I, r * sin(a), J));
		}
		endShape();
	}

	void collar(Pv3D.pt P, Pv3D.vec V, float r, float rd) {
		Pv3D.vec I = U(Normal(V));
		Pv3D.vec J = U(N(I, V));
		collar(P, V, I, J, r, rd);
	}

	void collar(Pv3D.pt P, Pv3D.vec V, Pv3D.vec I, Pv3D.vec J, float r, float rd) {
		float da = TWO_PI / 36;
		beginShape(QUAD_STRIP);
		for (float a = 0; a <= TWO_PI + da; a += da) {
			v(P(P, r * cos(a), I, r * sin(a), J, 0, V));
			v(P(P, rd * cos(a), I, rd * sin(a), J, 1, V));
		}
		endShape();
	}

	void cone(Pv3D.pt P, Pv3D.vec V, float r) {
		fan(P, V, r);
		disk(P, V, r);
	}

	void stub(Pv3D.pt P, Pv3D.vec V, float r, float rd) {
		collar(P, V, r, rd);
		disk(P, V, r);
		disk(P(P, V), V, rd);
	}

	void arrow(Pv3D.pt P, Pv3D.vec V, float r) {
		stub(P, V(.8f, V), r * 2 / 3, r / 3);
		cone(P(P, V(.8f, V)), V(.2f, V), r);
	}

// **************************** PRIMITIVE
	void showFrame(float d) {
		noStroke();
		fill(metal);
		sphere(d / 10);
		fill(blue);
		showArrow(d, d / 10);
		fill(red);
		pushMatrix();
		rotateY(PI / 2);
		showArrow(d, d / 10);
		popMatrix();
		fill(green);
		pushMatrix();
		rotateX(-PI / 2);
		showArrow(d, d / 10);
		popMatrix();
	}

	void showFan(float d, float r) {
		float da = TWO_PI / 36;
		beginShape(TRIANGLE_FAN);
		vertex(0, 0, d);
		for (float a = 0; a <= TWO_PI + da; a += da) {
			vertex(r * cos(a), r * sin(a), 0);
		}
		endShape();
	}

	void showCollar(float d, float r, float rd) {
		float da = TWO_PI / 36;
		beginShape(QUAD_STRIP);
		for (float a = 0; a <= TWO_PI + da; a += da) {
			vertex(r * cos(a), r * sin(a), 0);
			vertex(rd * cos(a), rd * sin(a), d);
		}
		endShape();
	}

	void showCone(float d, float r) {
		showFan(d, r);
		showFan(0, r);
	}

	void showStub(float d, float r, float rd) {
		showCollar(d, r, rd);
		showFan(0, r);
		pushMatrix();
		translate(0, 0, d);
		showFan(0, rd);
		popMatrix();
	}

	void showArrow() {
		showArrow(1, 0.08f);
	}

	void showArrow(float d, float r) {
		float dd = d / 5;
		showStub(d - dd, r * 2 / 3, r / 3);
		pushMatrix();
		translate(0, 0, d - dd);
		showCone(dd, r);
		popMatrix();
	}

	void showBlock(float w, float d, float h, float x, float y, float z, float a) {
		pushMatrix();
		translate(x, y, h / 2);
		rotateZ(TWO_PI * a);
		box(w, d, h);
		popMatrix();
	}

//*********** PICK
	Pv3D.vec I = V(1, 0, 0), J = V(0, 1, 0), K = V(0, 0, 1); // screen projetions of global model frame

	void computeProjectedVectors() {
		Pv3D.pt O = ToScreen(P(0, 0, 0));
		Pv3D.pt A = ToScreen(P(1, 0, 0));
		Pv3D.pt B = ToScreen(P(0, 1, 0));
		Pv3D.pt C = ToScreen(P(0, 0, 1));
		I = V(O, A);
		J = V(O, B);
		K = V(O, C);
	}

	Pv3D.vec ToIJ(Pv3D.vec V) {
		float x = det2(V, J) / det2(I, J);
		float y = det2(V, I) / det2(J, I);
		return V(x, y, 0);
	}

	Pv3D.vec ToK(Pv3D.vec V) {
		float z = dot(V, K) / dot(K, K);
		return V(0, 0, z);
	}

// ******************************************COLORS 
	public static int black = 0xff000000, white = 0xffFFFFFF, // set more colors using Menu >  Tools > Color Selector
			red = 0xffFF0000, green = 0xff00FF01, blue = 0xff0300FF, yellow = 0xffFEFF00, cyan = 0xff00FDFF, magenta = 0xffFF00FB,
			grey = 0xff818181, orange = 0xffFFA600, brown = 0xffB46005, metal = 0xffB5CCDE, dgreen = 0xff157901,
			blue50 = 0x800300FF;

	void pen(int c, float w) {
		stroke(c);
		strokeWeight(w);
	}

// ******************************** TEXT , TITLE, and USER's GUIDE
	Boolean scribeText = true; // toggle for displaying of help text

	void scribe(String S, float x, float y) {
		fill(0);
		text(S, x, y);
		noFill();
	} // writes on screen at (x,y) with current fill color

	void scribeHeader(String S, int i) {
		fill(0);
		text(S, 10, 20 + i * 20);
		noFill();
	} // writes black at line i

	void scribeHeaderRight(String S) {
		fill(0);
		text(S, width - 7.5f * S.length(), 20);
		noFill();
	} // writes black on screen top, right-aligned

	void scribeFooter(String S, int i) {
		fill(0);
		text(S, 10, height - 10 - i * 20);
		noFill();
	} // writes black on screen at line i from bottom

	void scribeAtMouse(String S) {
		fill(0);
		text(S, mouseX, mouseY);
		noFill();
	} // writes on screen near mouse

	void scribeMouseCoordinates() {
		fill(black);
		text("(" + mouseX + "," + mouseY + ")", mouseX + 7, mouseY + 25);
		noFill();
	}
}
