package generic;

/**
 * Created by Gregary on 11/9/2016.
 */
public class Scalar implements Value {

    private final double value;

    public static Scalar[] array(double... vals) {
        Scalar[] out = new Scalar[vals.length];
        for (int i = 0; i < vals.length; i++) {
            out[i] = new Scalar(vals[i]);
        }
        return out;
    }

    public Scalar(double val) {
        this.value = val;
    }

    @Override
    public Value negate() {
        return new Scalar(-this.value);
    }

    @Override
    public Value reciprocal() {
        return new Scalar(1.0/this.value);
    }

    @Override
    public Value add(Value other) {
        if (other instanceof Scalar) return this.add((Scalar) other);
        try {
            return other.add(this);
        } catch (StackOverflowError e) {
            throw new RuntimeException("Addition of "+this.getClass().getName()
                    +" and "+other.getClass().getName()+" is not implemented.");
        }
    }

    @Override
    public Value multiply(Value other) {
        if (other instanceof Scalar) return this.multiply((Scalar) other);
        try {
            return other.multiply(this);
        } catch (StackOverflowError e) {
            throw new RuntimeException("Multiplication of "+this.getClass().getName()
                    +" and "+other.getClass().getName()+" is not implemented.");
        }
    }

    public Scalar add(Scalar other) {
        return new Scalar(this.value + other.value);
    }

    public Scalar multiply(Scalar other) {
        return new Scalar(this.value * other.value);
    }

    public double value() {
        return value;
    }

    public String toString() {
        return Double.toString(value);
    }

}
