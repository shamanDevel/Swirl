/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package swirl;

import java.util.Objects;

/**
 *
 * @author Sebastian Weiss
 */
public class Frame {
	public final Vector3f P;
	public final Vector3f I,J,K;

	public Frame() {
		this.P = new Vector3f();
		this.I = new Vector3f();
		this.J = new Vector3f();
		this.K = new Vector3f();
	}

	public Frame(Vector3f P, Vector3f I, Vector3f J, Vector3f K) {
		this.P = new Vector3f(P);
		this.I = new Vector3f(I);
		this.J = new Vector3f(J);
		this.K = new Vector3f(K);
	}

	@Override
	protected Frame clone() {
		return new Frame(P, I, J, K);
	}
	
	public Vector3f[] getPickablePoints() {
		return new Vector3f[]{P, I.add(P), J.add(P), K.add(P)};
	}
	public void movePickablePoint(int index, Vector3f movement) {
		switch (index) {
			case 0:
				P.addLocal(movement);
				break;
			case 1:
				I.addLocal(movement);
				I.cross(K.normalize(), J);
				J.negateLocal();
				I.cross(J.normalize(), K);
				break;
			case 2:
				J.addLocal(movement);
				J.cross(I.normalize(), K);
				K.negateLocal();
				J.cross(K.normalize(), I);
				break;
			case 3:
				K.addLocal(movement);
				K.cross(J.normalize(), I);
				I.negateLocal();
				K.cross(I.normalize(), J);
				break;
		}
	}

	@Override
	public int hashCode() {
		int hash = 7;
		hash = 29 * hash + Objects.hashCode(this.P);
		hash = 29 * hash + Objects.hashCode(this.I);
		hash = 29 * hash + Objects.hashCode(this.J);
		hash = 29 * hash + Objects.hashCode(this.K);
		return hash;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}
		final Frame other = (Frame) obj;
		if (!Objects.equals(this.P, other.P)) {
			return false;
		}
		if (!Objects.equals(this.I, other.I)) {
			return false;
		}
		if (!Objects.equals(this.J, other.J)) {
			return false;
		}
		if (!Objects.equals(this.K, other.K)) {
			return false;
		}
		return true;
	}

	@Override
	public String toString() {
		return "Frame{" + "P=" + P + ", I=" + I + ", J=" + J + ", K=" + K + '}';
	}
}
