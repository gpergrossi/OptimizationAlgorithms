package main;

import calc.NPolynomial;
import generic.Scalar;
import generic.Vector;

public class Main {
	
	public static void main(String[] args) {

        //P10 test F1
//        NPolynomial func = NPolynomial.fromString("x[1]^2 + x[2]^2 + x[3]^2");
//        Vector<Scalar> x = new Vector<>(Scalar.array(1, 1, 1));

        //P10 test F2
		NPolynomial func = NPolynomial.fromString("x[1]^2 + 2*x[2]^2 - 2*x[1]*x[2] - 2*x[2]");
        Vector<Scalar> x = new Vector<>(Scalar.array(1, 0.5));

        //P10 test F3
//        NPolynomial func = NPolynomial.fromString("100*x[2]^2 - 200*x[1]^2*x[2] + 100*x[1]^4 + 1 - 2*x[1] + x[1]^2");
//        Vector<Scalar> x = new Vector<>(Scalar.array(-1.2, 1));

        //P10 test F4
//        NPolynomial func = NPolynomial.fromString("x[1]^4 + 4*x[1]^3*x[2] + 6*x[1]^2*x[2]^2 + 4*x[1]*x[2]^3 + x[2]^4 + x[2]^2");
//        Vector<Scalar> x = new Vector<>(Scalar.array(2, -2));

        //P10 test F4
//        Scalar c = new Scalar(10000);
//        NPolynomial func_1 = NPolynomial.fromString("x[1]^2 - 2*x[1] + x[2]^2 - 2*x[2] + 2");
//        NPolynomial func_2 = NPolynomial.fromString("x[1]^4 + 2*x[1]^2*x[2]^2 - 0.5*x[1]^2 + x[2]^4 - 0.5*x[2]^2 + 0.0625");
//        NPolynomial func = func_2.multiply(c).add(func_1);
//        Vector<Scalar> x = new Vector<>(Scalar.array(1, -1));

        System.out.println(func);
        System.out.println(func.gradient(x.getSize()));
        System.out.println(func.hessian(x.getSize()));

        SteepestDescent solve = new SteepestDescent(func, x);
        //NewtonMethod solve = new NewtonMethod(func, x);
        //BFGSQuasiNewton solve = new BFGSQuasiNewton(func, x);
        //ConjugateGradient solve = new ConjugateGradient(func, x);

        long nanos = System.nanoTime();
        while (!solve.isDone()) {
            Vector<Scalar> xk = solve.iterate();
            System.out.println(xk);
        }
		long time = System.nanoTime() - nanos;

        double s = time / 1000000000.0;
        System.out.println("Finished in "+s+" seconds");


        // step = dt = -H(xt)^-1*g(xt)
        // H(xt)*dt = -g(xt)
	}
	
}
