package calc;

import generic.Scalar;
import generic.Value;
import generic.Vector;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class NPolynomial implements NFunction {

	final List<NPolyTerm> parts;
	
	/**
	 * Format is "coefficient * x[0]^p1 * x[2]^p2 * x[3]^p3 + ..."
	 * only valid variables are x[0], ..., x[n]; n is an integer
	 * @param form
	 * @return 
	 */
	public static NPolynomial fromString(String form) throws RuntimeException {
		List<NPolyTerm> build = new ArrayList<NPolyTerm>();
		
		// Remove all whitespace
		form = form.replaceAll("\\s", "");
		
		// Divide the string just in front of each + or - using Regex look-ahead
		String[] split = form.split("(?=[+-])");
		
		for (int i = 0; i < split.length; i++) {
			if (split[i].startsWith("+")) {
				if (i == 0) {
					throw new RuntimeException("Bad form, first element has '+'");
				} else {
					split[i] = split[i].substring(1);
				}
			}
			NPolyTerm prod = NPolyTerm.fromString(split[i]);
			if (!prod.isZero()) build.add(prod);
		}
		
		return new NPolynomial(build);
	}
	
	private static List<NPolyTerm> combineTerms(List<NPolyTerm> parts) {
		List<NPolyTerm> combined = new ArrayList<NPolyTerm>();
		List<Integer> skipIndex = new ArrayList<Integer>();
		
		for (int i = 0; i < parts.size(); i++) {
			if (skipIndex.contains(i)) continue;
			NPolyTerm term = parts.get(i);
			double coef = term.coefficient;
			for (int j = i+1; j < parts.size(); j++) {
				NPolyTerm other = parts.get(j);
				if (term.canCombine(other)) {
					skipIndex.add(j);
					coef += other.coefficient;
				}
			}
			combined.add(new NPolyTerm(coef, term.powersBeginIndex, term.powers));
		}
		
		return combined;
	}
	
	public NPolynomial(List<NPolyTerm> parts) {
		// add parts to a private list
		List<NPolyTerm> build = combineTerms(parts);
		build.sort(NPolyTerm.SORT_ORDER);
		
		// only remaining reference to the private list is unmodifiable
		this.parts = Collections.unmodifiableList(build);
	}
	
	public NPolynomial derivative(int varIndex) {
		List<NPolyTerm> build = new ArrayList<>();
		for (NPolyTerm prod : parts) {
			NPolyTerm deriv = prod.derivative(varIndex);
			if (!deriv.isZero()) build.add(deriv);
		}
		if (build.size() == 0) build.add(NPolyTerm.ZERO);
		return new NPolynomial(build);		
	}
	
	public Scalar value(Vector<Scalar> x) {
		double sum = 0;
		for (NPolyTerm part : parts) {
			sum += part.value(x).value();
		}
		return new Scalar(sum);
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		int num = 0;
		for (NPolyTerm prod : parts) {
			sb.append(prod);
			num++;
			if (num < parts.size()) sb.append(" + ");
		}
		return sb.toString(); 
	}

	public NPolynomial negate() {
        List<NPolyTerm> build = new ArrayList<>();
        for (NPolyTerm prod : parts) {
            NPolyTerm negate = prod.negate();
            build.add(negate);
        }
        return new NPolynomial(build);
    }

    public NPolynomial reciprocal() {
        List<NPolyTerm> build = new ArrayList<>();
        for (NPolyTerm prod : parts) {
            NPolyTerm recip = prod.reciprocal();
            build.add(recip);
        }
        return new NPolynomial(build);
    }
	
	public Value add(Value other) {
        if (other instanceof NPolynomial) return this.add((NPolynomial) other);
        if (other instanceof NPolyTerm) return this.add((NPolyTerm) other);
        if (other instanceof Scalar) return this.add((Scalar) other);
        try {
		    return other.add(this);
        } catch (StackOverflowError e) {
            throw new RuntimeException("Addition of "+this.getClass().getName()
                    +" and "+other.getClass().getName()+" is not implemented.");
        }
	}

	public NPolynomial add(NPolynomial other) {
		List<NPolyTerm> newParts = new ArrayList<NPolyTerm>();
		for (NPolyTerm term : parts) newParts.add(term);
		for (NPolyTerm term : other.parts) newParts.add(term);
		return new NPolynomial(newParts);
	}
	
	public NPolynomial add(NPolyTerm other) {
		List<NPolyTerm> newParts = new ArrayList<NPolyTerm>();
		for (NPolyTerm term : parts) newParts.add(term);
		newParts.add(other);
		return new NPolynomial(newParts);
	}

    public NPolynomial add(Scalar other) {
        return this.add(new NPolyTerm(other.value()));
    }
	
	public Value multiply(Value other) {
        if (other instanceof NPolynomial) return this.multiply((NPolynomial) other);
        if (other instanceof NPolyTerm) return this.multiply((NPolyTerm) other);
        if (other instanceof Scalar) return this.multiply((Scalar) other);
        try {
            return other.multiply(this);
        } catch (StackOverflowError e) {
            throw new RuntimeException("Multiplication of "+this.getClass().getName()
                    +" and "+other.getClass().getName()+" is not implemented.");
        }
	}

	public NPolynomial multiply(NPolynomial other) {
		NPolynomial sumPoly = null;
		for (NPolyTerm term : other.parts) {
			NPolynomial newPoly = new NPolynomial(this.parts);
			newPoly = newPoly.multiply(term);
			
			if (sumPoly == null) sumPoly = newPoly;
			else sumPoly = sumPoly.add(newPoly);
		}
		
		return sumPoly;
	}
	
	public NPolynomial multiply(NPolyTerm other) {
		List<NPolyTerm> newParts = new ArrayList<NPolyTerm>();
		
		if (!other.isZero()) {
			for (NPolyTerm term : parts) {
				NPolyTerm newTerm = term.multiply(other);
				if (!newTerm.isZero()) newParts.add(newTerm);
			}
		}
		if (newParts.size() == 0) newParts.add(NPolyTerm.ZERO);
	
		return new NPolynomial(newParts);
	}

    public NPolynomial multiply(Scalar other) {
        return this.multiply(new NPolyTerm(other.value()));
    }

	public NFunctionVector gradient(int maxIndex) {
		NPolynomial[] values = new NPolynomial[maxIndex];
		for (int i = 1; i <= maxIndex; i++) {
			values[i-1] = this.derivative(i);
		}
		return new NFunctionVector(values);
	}

	public NFunctionMatrix hessian(int maxIndex) {
		NFunction[][] hessian = new NFunction[maxIndex][maxIndex];
		for (int r = 0; r < maxIndex; r++) {
			for (int c = 0; c < maxIndex; c++) {
				hessian[r][c] = this.derivative(r+1).derivative(c+1);
			}
		}
		return new NFunctionMatrix(hessian);
	}

	
}
