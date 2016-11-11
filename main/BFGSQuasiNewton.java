package main;

import calc.NFunction;
import calc.NFunctionMatrix;
import calc.NFunctionVector;
import generic.Matrix;
import generic.Scalar;
import generic.Value;
import generic.Vector;

public class BFGSQuasiNewton {

    int iteration;
    NFunction func;

    Vector<Scalar> xk;
    Matrix<Scalar> Bk;

    NFunctionVector gradFunc;

    boolean done;

    public BFGSQuasiNewton(NFunction func, Vector initial) {
        this.func = func;
        this.xk = initial;
        this.Bk = Matrix.identity(initial.getSize());
        this.gradFunc = func.gradient(initial.getSize());
        this.iteration = 0;
        this.done = false;
    }

    private final double SMALL    = Double.MIN_VALUE;
    private final double TAU      = 0.5;        // Reduction in step size for each attempt
    private final double EPSILON  = 0.0000001;  // Ending Epsilon
    private final double WOLFE_C1 = 0.0001;     // Sufficient reduction in f(x)
    private final double WOLFE_C2 = 0.9;        // c1 <= c2 <= 1. Sufficient reduction in g(x)
    private final double BETA    = 0.0001;      // Backup line search Beta

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
            System.out.println("x"+iteration+" = "+xk);
            System.out.println("f(x"+iteration+") = "+value);
            System.out.println("B"+iteration+" = \n"+Bk);
            iteration++;
            return xk;
        }

        // Calculate the Search Direction = -(Bk^-1)*g(xk)
        Vector<Scalar> gk = gradFunc.value(xk);
        Vector<Scalar> pk = Bk.reciprocal().multiply(gk.negate()).getCol(0).toScalars();
        System.out.println("p"+iteration+" = "+pk);

        // Wolfe condition 1 (Armijo)
        Scalar pkT_gk = pk.transpose().multiply(gk.toValues()).toScalars().get(0, 0); // (pk_T)*gk
        Scalar armijoCoef = new Scalar(WOLFE_C1).multiply(pkT_gk); // c1*(pk_T)*gk

        // Wolfe condition 2
        double c2_pkT_gk = new Scalar(WOLFE_C2).multiply(pkT_gk).value();

        // Calculate step length
        double tryStep;
        System.out.print("failed condition: ");
        for (tryStep = 1.0; tryStep >= SMALL; tryStep *= TAU) {
            Scalar ak = new Scalar(tryStep);
            Vector<Scalar> xt = xk.add(pk.multiply(ak)).toScalars();   // xt = xk + ak*pk

            // Wolfe condition 1
            double ft = func.value(xt).value();
            double armijo = value.add(armijoCoef.multiply(ak)).value();  // f(xk) + c1*ak*(pk_T)*gk
            if (ft-0.000001 > armijo) {
                System.out.print("1");
                continue; // Armijo condition: f(xk+ak*pk) <= f(xk) + c1*ak*pk*(gk_T)
            }

            // Wolfe condition 2
            double pk_gt = pk.multiply(gradFunc.value(xt).transpose().toValues()).toScalars().get(0, 0).value();
            if (pk_gt+0.000001 < c2_pkT_gk) {
                System.out.print("2");
                continue; // Wolfe condition 2, sufficient improvement in slope
            }

            break; // Both conditions met
        }
        System.out.println();
        if (tryStep <= SMALL) {
            // Calculate the Search Direction = -gradFunc(f)
            pk = gk.negate().toScalars();
            System.out.println("Using line search");
            System.out.println("p"+iteration+" = "+pk);

            // We are using the Armijo condition along with a backtracking search
            Scalar gkT_pk = gk.transpose().multiply(pk.toValues()).toScalars().get(0, 0);
            armijoCoef = gkT_pk.multiply(new Scalar(BETA));

            // Calculate step length
            for (tryStep = 1.0; tryStep >= SMALL; tryStep *= TAU) {
                Scalar ak = new Scalar(tryStep);
                Vector<Scalar> tryX = xk.add(pk.multiply(ak)).toScalars();   // xk + ak*pk
                double tryValue = func.value(tryX).value();               // f(xk + ak*pk)
                double armijo = value.add(armijoCoef.multiply(ak)).value();  // f(xk) + ak*BETA*transpose(gk)*pk
                if (tryValue <= armijo) break; // Armijo condition: f(xk+ak*pk) <= f(xk) + ak*BETA*transpose(gk)*pk
            }
        }

        Scalar ak = new Scalar(tryStep);
        System.out.println("a"+iteration+" = "+ak);

        // Update the pos
        Vector<Scalar> xk1 = xk.add(pk.multiply(ak)).toScalars();
        System.out.println("x"+iteration+" = "+xk1);
        System.out.println("f(x"+iteration+") = "+func.value(xk1));

        // Update Bk
        Vector<Scalar> gk1 = gradFunc.value(xk1);
        Vector<Scalar> sk = xk.negate().add(xk1.toValues()).toScalars();
        Vector<Scalar> yk = gk.negate().add(gk1.toValues()).toScalars();

        // Update term 1
        Vector<Scalar> Bk_sk = Bk.multiply(sk.toValues()).toScalars().asVector();
        Matrix<Scalar> Bk_sk_skT = Bk_sk.multiply(sk.transpose().toValues()).toScalars();
        Scalar skT_Bk_sk = sk.transpose().toValues().multiply(Bk_sk.toValues()).toScalars().get(0, 0);
        Matrix<Scalar> term1 = Bk_sk_skT.multiply(Bk.toValues()).multiply(skT_Bk_sk.reciprocal()).negate().toScalars();

        // Update term 2
        Matrix<Scalar> yk_ykT = yk.multiply(yk.transpose().toValues()).toScalars();
        Scalar ykT_sk = yk.transpose().multiply(sk.toValues()).toScalars().get(0, 0);
        Matrix<Scalar> term2 = yk_ykT.multiply(ykT_sk.reciprocal()).toScalars();

        // Do update
        Bk = Bk.add(term1.toValues()).add(term2.toValues()).toScalars();
        System.out.println("B"+iteration+" = \n"+Bk);

        // Done?
        double norm = ((Scalar) gk1.norm2()).value();
        double f_xk1 = func.value(xk1).value();
        double end = norm / (1 + Math.abs(f_xk1));
        if (end < EPSILON) {
            System.out.println("Epsilon condition!");
            done = true;
        }

        // Limit of double precision
        Value delta = sk.norm2();
        if (((Scalar) delta).value() < SMALL) {
            System.out.println("Max precision of double arithmetic");
            done = true;
        }

        iteration++;
        xk = xk1;
        return xk;
    }

}
