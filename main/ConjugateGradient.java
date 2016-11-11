package main;

import calc.NFunction;
import calc.NFunctionVector;
import generic.Scalar;
import generic.Value;
import generic.Vector;

public class ConjugateGradient {

    int iteration;
    NFunction func;
    Vector<Scalar> lastDx;
    Vector<Scalar> lastSn;
    Vector<Scalar> xk;
    NFunctionVector gradFunc;
    boolean done;
    boolean wasReset;

    public ConjugateGradient(NFunction func, Vector initial) {
        this.func = func;
        this.gradFunc = func.gradient(initial.getSize());
        this.lastDx = this.lastSn = initial.multiply(new Scalar(0));
        this.xk = initial;
        this.iteration = 0;
        this.done = false;
        this.wasReset = true;
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

    public Vector<Scalar> iterate() {
        if (done) {
            System.out.println("=== Done ===");
            return xk;
        }

        System.out.println("=== Begin Iteration "+iteration+" ===");
        Scalar value = func.value(xk);
        if (iteration == 0) {
            System.out.println("x"+iteration+" = "+ this.xk);
            System.out.println("f(x"+iteration+") = "+value);
            iteration++;
            return this.xk;
        }

        // Calculate the Search Direction = -gradFunc(f)
        Vector<Scalar> gk = gradFunc.value(xk);
        Vector<Value> dx = gk.negate();
        System.out.println("dx"+iteration+" = "+dx);

        Vector<Value> sn;
        if (!wasReset) {
            // Calculate Beta n, using Fletcher-Reeves formula
            Scalar pkT_pk = dx.transpose().multiply(dx).toScalars().get(0, 0);
            Scalar last_pkT_pk = lastDx.transpose().multiply(lastDx.toValues()).toScalars().get(0, 0);
            Scalar Bn = (Scalar) pkT_pk.multiply(last_pkT_pk.reciprocal());

            sn = dx.add(lastSn.toValues().multiply(Bn));
            System.out.println("Bn" + iteration + " = " + Bn);
        } else {
            System.out.println("Reset using gradient search");
            sn = dx;
        }
        System.out.println("sn"+iteration+" = "+sn);

        lastDx = dx.toScalars();
        lastSn = sn.toScalars();

        // We are using the Armijo condition along with a backtracking search
        Scalar gkT_sn = gk.transpose().multiply(sn).toScalars().get(0, 0);
        Scalar armijoCoef = gkT_sn.multiply(new Scalar(BETA));

        // Calculate step length
        Scalar ak = null;
        for (double tryStep = 1.0; tryStep >= SMALL; tryStep *= TAU) {
            ak = new Scalar(tryStep);
            Vector<Scalar> tryX = xk.add(sn.multiply(ak)).toScalars();   // xk + ak*pk
            double tryValue = func.value(tryX).value();               // f(xk + ak*pk)
            double armijo = value.add(armijoCoef.multiply(ak)).value();  // f(xk) + ak*BETA*transpose(gk)*pk
            if (tryValue <= armijo) break; // Armijo condition: f(xk+ak*pk) <= f(xk) + ak*BETA*transpose(gk)*pk
        }
        System.out.println("a"+iteration+" = "+ak);

        // Update the xk
        Vector<Scalar> xk1 = this.xk.add(sn.multiply(ak)).toScalars();
        System.out.println("x"+iteration+" = "+ xk1);
        System.out.println("f(x"+iteration+") = "+func.value(xk1));

        // Done?
        gk = gradFunc.value(xk1);
        double norm = ((Scalar) gk.norm2()).value();
        double fx = func.value(xk1).value();
        double end = norm / (1 + Math.abs(fx));
        if (end < EPSILON) {
            System.out.println("Epsilon condition!");
            done = true;
        }

        // Limit of double precision
        Value delta = xk.negate().add(xk1.toValues()).norm2();
        if (((Scalar) delta).value() == 0) {
            if (wasReset) {
                System.out.println("Max precision of double arithmetic");
                done = true;
            } else {
                System.out.println("No progress made, retrying with gradient search");
                wasReset = true;
            }
        } else {
            wasReset = false;
        }

        this.xk = xk1;
        iteration++;
        return this.xk;
    }

}
