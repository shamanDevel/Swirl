/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package swirl;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import junit.framework.AssertionFailedError;

import static org.junit.Assert.*;

/**
 *
 * @author Sebastian Weiss
 */
public class InterpolationTestHelper {
	private static final double EPSILON = 0.001;
	private static final Random RAND = new Random();
	
	public static void assertRealVector(String vectorName, Vector3f vector) {
		String msg = "Vector "+vectorName+" "+vector+" is not real";
		assertFalse(msg, Float.isInfinite(vector.x));
		assertFalse(msg, Float.isNaN(vector.x));
		assertFalse(msg, Float.isInfinite(vector.y));
		assertFalse(msg, Float.isNaN(vector.y));
		assertFalse(msg, Float.isInfinite(vector.z));
		assertFalse(msg, Float.isNaN(vector.z));
	}
	
	public static void assertFrameValid(Frame frame) {
		assertRealVector("P", frame.P);
		assertRealVector("I", frame.I);
		assertRealVector("J", frame.J);
		assertRealVector("K", frame.K);
		
		float l1 = frame.I.lengthSquared();
		float l2 = frame.J.lengthSquared();
		float l3 = frame.K.lengthSquared();
		assertEquals("I,J,K are not of the same length", l1, l2, EPSILON);
		assertEquals("I,J,K are not of the same length", l1, l3, EPSILON);
		assertTrue("I,J,K are zero", l1 >= EPSILON);
		
		assertEquals("I,J are not orthogonal", 0, frame.I.dot(frame.J), EPSILON);
		assertEquals("I,K are not orthogonal", 0, frame.I.dot(frame.K), EPSILON);
		assertEquals("J,K are not orthogonal", 0, frame.J.dot(frame.K), EPSILON);
		Vector3f K2 = frame.I.cross(frame.J);
		assertEquals("Frame is not a right-handed system", K2.normalize(), frame.K.normalize());
	}
	
	public static void assertFrameEquals(Frame a, Frame b) {
		assertEquals("P", a.P, b.P);
		assertEquals("I", a.I, b.I);
		assertEquals("J", a.J, b.J);
		assertEquals("K", a.K, b.K);
	}
	
	public static void assertInterpolating(Interpolation i, Frame start, Frame end) {
		i.setStartEnd(start, end);
		Frame f = new Frame();
		i.interpolate(0, f);
		assertFrameEquals(start, f);
		i.interpolate(1, f);
		assertFrameEquals(end, f);
	}
	
	public static void assertValidFrames(Interpolation inter, Frame start, Frame end) {
		inter.setStartEnd(start, end);
		List<Float> points = new ArrayList<>(); //sample points
		points.add(0f);
		points.add(1f);
		points.add(2f);
		points.add(-1f);
		points.add(0.5f);
		for (int i=0; i<20; ++i) {
			points.add(RAND.nextFloat()); //interpolate
		}
		for (int i=0; i<20; ++i) {
			points.add(RAND.nextFloat() * 9 + 1); //Extrapolate to 10 in the future
		}
		for (int i=0; i<20; ++i) {
			points.add(RAND.nextFloat() * -10); //Extrapolate to 10 in the past
		}
		//check points
		Frame f = new Frame();
		for (float t : points) {
			inter.interpolate(t, f);
			try {
				assertFrameValid(f);
			} catch (AssertionError e) {
				throw new AssertionError("Illegal frame at time "+t+": "+f+"\n"+e.getMessage(), e);
			}
		}
	}
	
	public static void assertSteadyMotion(Interpolation iter, Frame start, Frame end) {
		float s = 0.2f;
		float e = 0.8f;
		Frame sf = new Frame();
		Frame ef = new Frame();
		float points[] = new float[]{0.3f, 0.4f, 0.5f, 0.6f, 0.7f};
		//when I interpolate from 0 to 1, reaching these points, it must be the same
		//frames as when I interpolate from 0.2 to 0.8, reaching the same points 
		//(in the different scale, of course).
		iter.setStartEnd(start, end);
		Frame[] frames = new Frame[points.length];
		for (int i=0; i<points.length; ++i) {
			frames[i] = new Frame();
			iter.interpolate(points[i], frames[i]);
		}
		iter.interpolate(s, sf);
		iter.interpolate(e, ef);
		//check the frames against the lower scale
		iter.setStartEnd(sf, ef);
		Frame f = new Frame();
		for (int i=0; i<points.length; ++i) {
			float t = (points[i]-s) / (e-s);
			iter.interpolate(t, f);
			try {
				assertFrameEquals(frames[i], f);
			} catch (AssertionError ex) {
				throw new AssertionError("Illegal frame at time "+t+"("+points[i]+")\n"+ex.getMessage(), ex);
			}
		}
	}
}
