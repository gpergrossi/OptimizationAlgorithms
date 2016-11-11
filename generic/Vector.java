package generic;

import java.util.function.Function;

/**
 * Created by Gregary on 11/9/2016.
 */
public class Vector<T extends Value> extends Matrix<T> {

    public Vector(T[] values) {
        super(values.length, 1, values, false);
    }

    public Vector(T[] values, boolean transpose) {
        super(values.length, 1, values, transpose);
    }

    public Matrix<T> asMatrix() {
        return this;
    }

    public Vector<T> transpose() {
        return new Vector<>(this.values, !this.transpose);
    }

    public Vector<Scalar> toScalars() {
        return this.apply(Scalar.class, (val) -> {
            if (val instanceof Scalar) return (Scalar) val;
            else throw new RuntimeException("Value is not scalar compatible: \""+val.toString()+"\"");
        });
    }

    public Vector<Value> toValues() {
        return this.apply(Value.class, (val) -> (Value) val);
    }

    /**
     * only works on Vector<Scalar>
     */
    public Value norm2() {
        Vector<Scalar> scalars = this.toScalars();
        double d = 0;
        for (int i = 0; i < scalars.getSize(); i++) {
            double v = scalars.get(i).value();
            d += v*v;
        }
        return new Scalar(Math.sqrt(d));
    }

    public boolean isColumn() {
        return getCols() == 1;
    }

    public boolean isRow() {
        return getRows() == 1;
    }

    public <R extends Value> Vector<R> apply(Class<R> clazz, Function<T, R> function) {
        return super.apply(clazz, function).asVector();
    }

    public Vector<Value> negate() {
        return super.negate().asVector();
    }

    public Vector<Value> add(Vector<Value> other) {
        return super.add(other).asVector();
    }

    public Vector<Value> multiply(Scalar other) {
        return super.multiply(other).asVector();
    }

    public Value dot(Vector<Value> other) {
        if (this.getSize() != other.getSize()) throw new RuntimeException("Vector dot product bad dimensions");
        Value sum = null;
        for (int i = 0; i < getSize(); i++) {
            Value product = this.get(i).multiply(other.get(i));
            if (sum == null) sum = product;
            else sum = sum.add(product);
        }
        return sum;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        if (isColumn()) sb.append('<');
        else sb.append('[');

        boolean separate = false;
        for (int i = 0; i < getSize(); i++) {
            if (separate) sb.append(", ");
            sb.append(get(i,0));
            separate = true;
        }

        if (isColumn()) sb.append('>');
        else sb.append(']');

        return sb.toString();
    }

}
