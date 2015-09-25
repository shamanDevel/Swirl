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
public interface Interpolation {
	
	void setStartEnd(Frame start, Frame end);
	
	void interpolate(float t, Frame toSet);
	
	void debugDraw(Swirl swirl);
}
