package calc;

import generic.Value;
import generic.Vector;
import generic.Scalar;

public interface NFunction extends Value {

	public NFunction derivative(int varIndex);
	
	public Scalar value(Vector<Scalar> x);

	public NFunctionVector gradient(int maxIndex);

	public NFunctionMatrix hessian(int maxIndex);

}
