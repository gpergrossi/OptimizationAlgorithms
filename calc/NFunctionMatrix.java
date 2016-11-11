package calc;

import generic.Matrix;
import generic.Scalar;
import generic.Vector;

public class NFunctionMatrix extends Matrix<NFunction> {

    // Default constructor, copy array for security
    public NFunctionMatrix(NFunction[][] values) {
        super(values, false);
    }

    @SuppressWarnings("unchecked")
    public NFunctionMatrix(NFunction[][] values, boolean transpose) {
        super(values, transpose);
    }

	public Matrix<Scalar> value(Vector<Scalar> x) {
		return apply(Scalar.class, (nfunc) -> nfunc.value(x));
	}

}
