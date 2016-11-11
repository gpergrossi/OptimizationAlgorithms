package main;

import calc.NFunction;
import calc.NFunctionVector;
import generic.Value;
import generic.Vector;
import generic.Scalar;

public class SteepestDescent {

    int iteration;
    NFunction func;
    Vector<Scalar> guess;
    NFunctionVector gradient;
    boolean done;

    public SteepestDescent(NFunction func, Vector initial) {
        this.func = func;
        this.guess = initial;
        this.gradient = func.gradient(initial.getSize());
        this.iteration = 0;
        this.done = false;
    }

    private final double SMALL   = Double.MIN_VALUE;
    private final double TAU     = 0.5;             // Reduction in step size for each attempt
    private final double BETA    = 0.0001;          // Sufficient reduction in f(x)
    private final double EPSILON = 0.0000001;       // Ending Epsilon

    public boolean isDone() {
        return done;
    }

    public int getIteration() {
        return iteration;
    }

    public static boolean PRINT = false;

    public Vector<Scalar> iterate() {
        if (done) {
            System.out.println("=== Done ===");
            return guess;
        }

        if (PRINT) System.out.println("=== Begin Iteration "+iteration+" ===");
        Vector<Scalar> xk = guess;
        Scalar value = func.value(xk);
        if (iteration == 0) {
            System.out.println("x"+iteration+" = "+guess);
            System.out.println("f(x"+iteration+") = "+value);
            iteration++;
            return guess;
        }

        // Calculate the Search Direction = -gradFunc(f)
        Vector<Scalar> gk = gradient.value(xk);
        Vector<Value> pk = gk.negate();
        if (PRINT) System.out.println("p"+iteration+" = "+pk);

        // We are using the Armijo condition along with a backtracking search
        Scalar gkT_pk = gk.transpose().multiply(pk).toScalars().get(0, 0);
        Scalar armijoCoef = gkT_pk.multiply(new Scalar(BETA));


        // Calculate step length
        double tryStep;
        for (tryStep = 1.0; tryStep >= SMALL; tryStep *= TAU) {
            Scalar ak = new Scalar(tryStep);
            Vector<Scalar> tryX = xk.add(pk.multiply(ak)).toScalars();   // xk + ak*pk
            double tryValue = func.value(tryX).value();               // f(xk + ak*pk)
            double armijo = value.add(armijoCoef.multiply(ak)).value();  // f(xk) + ak*BETA*transpose(gk)*pk
            if (tryValue <= armijo) break; // Armijo condition: f(xk+ak*pk) <= f(xk) + ak*BETA*transpose(gk)*pk
        }
        if (PRINT) System.out.println("a"+iteration+" = "+tryStep);
        Scalar ak = new Scalar(tryStep);

        // Update the xk
        guess = guess.add(pk.multiply(ak)).toScalars();
        if (PRINT) System.out.println("x"+iteration+" = "+guess);
        if (PRINT) System.out.println("f(x"+iteration+") = "+func.value(guess));

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
