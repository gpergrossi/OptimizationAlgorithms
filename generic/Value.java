package generic;

/**
 * Created by Gregary on 11/9/2016.
 */
public interface Value {

    Value negate();

    Value reciprocal();

    Value add(Value other);

    Value multiply(Value other);

}
