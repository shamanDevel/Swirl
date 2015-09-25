/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package swirl;

import org.junit.Test;

import static org.junit.Assert.*;
import static swirl.InterpolationTestHelper.*;
import static swirl.Swirl.DEG_TO_RAD;

/**
 *
 * @author Sebastian Weiss
 */
public class SwirlInterpolationTest {
	
	public SwirlInterpolationTest() {
	}
	
	@Test
	public void testSteadiness() {
		Frame startFrame = new Frame();
		Frame endFrame = new Frame();
		
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
		
		SwirlInterpolation i = new SwirlInterpolation();
		assertValidFrames(i, startFrame, endFrame);
		assertSteadyMotion(i, startFrame, endFrame);
	}
}
