/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package swirl;


public class LinearInterpolation implements Interpolation {
	private Frame start;
	private Frame end;

	@Override
	public void setStartEnd(Frame start, Frame end) {
		this.start = start;
		this.end = end;
	}

	@Override
	public void interpolate(float t, Frame toSet) {
		toSet.P.interpolateLocal(start.P, end.P, t);
		toSet.I.interpolateLocal(start.I, end.I, t);
		toSet.J.interpolateLocal(start.J, end.J, t);
		toSet.K.interpolateLocal(start.K, end.K, t);
	}

	@Override
	public void debugDraw(Swirl swirl) {
	}

	@Override
	public String debugString() {
		return "Linear Interpolation";
	}
	
}
