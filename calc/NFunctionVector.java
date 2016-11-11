package calc;

import generic.Scalar;
import generic.Vector;

public class NFunctionVector extends Vector<NFunction> {

    public NFunctionVector(NFunction... values) {
        super(values, false);
    }

    public NFunctionVector(NFunction[] values, boolean transpose) {
        super(values, transpose);
    }

    public Vector<Scalar> value(Vector<Scalar> x) {
        return apply(Scalar.class, (nfunc) -> nfunc.value(x));
    }

}
