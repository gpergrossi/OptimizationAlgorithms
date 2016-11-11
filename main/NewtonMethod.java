package main;

import calc.NFunction;
import calc.NFunctionMatrix;
import calc.NFunctionVector;
import generic.Matrix;
import generic.Scalar;
import generic.Value;
import generic.Vector;

public class NewtonMethod {

    int iteration;
    NFunction func;
    Vector<Scalar> guess;
    NFunctionVector gradient;
    NFunctionMatrix hessian;
    boolean done;

    public NewtonMethod(NFunction func, Vector initial) {
        this.func = func;
        this.guess = initial;
        this.gradient = func.gradient(initial.getSize());
        this.hessian = func.hessian(initial.getSize());
        this.iteration = 0;
        this.done = false;
    }

    private final double EPSILON = 0.0000001;     // Ending Epsilon

    public boolean isDone() {
        return done;
    }

    public int getIteration() {
        return iteration;
    }

    public Vector<Scalar> iterate() {
        if (done) {
            System.out.println("=== Done ===");
            return guess;
        }

        System.out.println("=== Begin Iteration "+iteration+" ===");
        Vector<Scalar> xk = guess;
        Scalar value = func.value(xk);
        if (iteration == 0) {
            System.out.println("x"+iteration+" = "+guess);
            System.out.println("f(x"+iteration+") = "+value);
            iteration++;
            return guess;
        }

        // Calculate the Search Direction
        Vector<Scalar> gk = gradient.value(xk);
        Matrix<Scalar> hk = hessian.value(xk);

        Vector<Value> pk;   Scalar ak;
        try {
            Matrix<Value> hki = hk.reciprocal();
            pk = hki.multiply(gk.negate()).getCol(0);
            ak = new Scalar(1);
        } catch (RuntimeException e) {
            // Steepest descent if non-invertible
            pk = gk.negate().toValues().asVector();
            ak = new Scalar(0.01);
            System.out.println("Non-invertible, using steepest descent");
        }
        System.out.println("p"+iteration+" = "+pk);
        System.out.println("a"+iteration+" = "+ak);

        // Update the xk
        guess = guess.add(pk.multiply(ak)).toScalars();
        System.out.println("x"+iteration+" = "+guess);
        System.out.println("f(x"+iteration+") = "+func.value(guess));

        // Done?
        gk = gradient.value(guess);
        double norm = ((Scalar) gk.norm2()).value();
        double fx = func.value(guess).value();
        double end = norm / (1 + Math.abs(fx));
        if (end < EPSILON) {
            System.out.println("Epsilon condition!");
            done = true;
        }

        // Limit of double precision
        Value delta = xk.negate().add(guess.toValues()).norm2();
        if (((Scalar) delta).value() == 0) {
            System.out.println("Max precision of double arithmetic");
            done = true;
        }

        iteration++;
        return guess;
    }

}
