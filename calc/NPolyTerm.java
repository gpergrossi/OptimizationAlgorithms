package calc;

import generic.Value;
import generic.Scalar;
import generic.Vector;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class NPolyTerm implements NFunction {

	public static final NPolyTerm ZERO = new NPolyTerm(0);
	
	int powersBeginIndex;
	final double coefficient;
	final double[] powers;
	
	public static final Comparator<NPolyTerm> SORT_ORDER = (o1, o2) -> {
        if (o1.powersBeginIndex < o2.powersBeginIndex) return -1;
        if (o1.powersBeginIndex > o2.powersBeginIndex) return 1;

        for (int i = 0; i < o1.powers.length; i++) {
            if (i >= o2.powers.length) return 1;
            if (o1.powers[i] > o2.powers[i]) return 1;
        }

        return 0;
    };

	/**
	 * Format is "coefficient * x1^p1 * x2^p2 * x3^p3 + ..."
	 * only valid variables are x1, x2, ..., x34, x35, ..., x[num]
	 * @param form
	 * @return 
	 */
	public static NPolyTerm fromString(String form) throws RuntimeException {
		double coefficient = 1;
		Map<Integer, Double> exponents = new HashMap<>();
		int largestSeenIndex = -1;
		int smallestSeenIndex = Integer.MAX_VALUE;
		
		// Remove all whitespace
		form = form.replaceAll("\\s", "");
		
		// Divide the string just in front of each * or / using Regex look-ahead
		String[] split = form.split("(?=[*\\/])");
		
		for (int i = 0; i < split.length; i++) {
			
			// Check operator, all parts after the first must have a * or /
			boolean negativeExponent = false;
			if (i == 0) {
				if (split[i].startsWith("[*\\/]")) {
					throw new RuntimeException("Bad form, first element has operator '"+split[i].charAt(0)+"'");
				}
			} else {
				if (split[i].startsWith("/")) negativeExponent = true;
				split[i] = split[i].substring(1);
			}
				
			
			// Get parts before and after ^
			String[] parts = split[i].split("\\^");
			if (parts.length > 2) { 
				throw new RuntimeException("Bad form, multiple '^' between * or /");
			}
						
			// Check for straight coefficients
			if (!parts[0].startsWith("x")) {
				try {
					double coeff = Double.parseDouble(parts[0]);
					coefficient *= coeff;
					continue;
				} catch(Exception e) {
					throw new RuntimeException("Bad form, variable or coefficient expected, got '"+parts[0]+"'");
				}
			}
						
			// Get variable index: 'x[this number]'
			if (!parts[0].startsWith("x[") || !parts[0].endsWith("]")) {
				throw new RuntimeException("Bad form, variable name must start with 'x[' and end with ']', got '"+parts[0]+"'");
			}
			parts[0] = parts[0].substring(2, parts[0].length()-1);
			int varIndex = 0;
			try {
				varIndex = Integer.parseInt(parts[0]);
			} catch(Exception e) {
				throw new RuntimeException("Bad form, variable expected integer value, received '"+parts[0]+"'");
			}
			if (varIndex < 1) {
				throw new RuntimeException("Bad form, negative or zero variable index");
			}
			if (varIndex > largestSeenIndex) largestSeenIndex = varIndex;
			if (varIndex < smallestSeenIndex) smallestSeenIndex = varIndex;
			
			// Get exponent value '^ [this number]'
			double exponent = 1;
			if (parts.length == 2) {
				try {
					exponent = Double.parseDouble(parts[1]);
				} catch(Exception e) {
					throw new RuntimeException("Bad form, exponent expected floating point value, received '"+parts[1]+"'");
				}
			}
			if (negativeExponent) exponent = -exponent;
			
			// Create an entry for the exponent of this variable index, later entries add together
			double currentExp = exponents.getOrDefault(varIndex, 0.0);
			exponents.put(varIndex, currentExp + exponent);
		}
		
		if (largestSeenIndex-smallestSeenIndex < 0) {
			return new NPolyTerm(coefficient);
		}
		
		// Create actual powers array
		double[] powers = new double[largestSeenIndex-smallestSeenIndex+1];
		for (Integer index : exponents.keySet()) {
			powers[index-smallestSeenIndex] = exponents.get(index);
		}
		
		return new NPolyTerm(coefficient, smallestSeenIndex, powers);
	}
	
	NPolyTerm(double coefficient, int beginIndex, double... powers) {
		this.coefficient = coefficient;
		this.powers = new double[powers.length];
		System.arraycopy(powers, 0, this.powers, 0, powers.length);
		this.powersBeginIndex = beginIndex;
	}
	
	public NPolyTerm(double coefficient, double... powers) {
		this(coefficient, 0, powers);
	}
	
	public boolean isZero() {
		return this.coefficient == 0;
	}
	
	public double getPower(int varIndex) {
		varIndex -= powersBeginIndex;
		if (varIndex < 0) return 0;
		if (varIndex >= powers.length) return 0;
		return powers[varIndex];
	}
	
	public boolean canCombine(NPolyTerm other) {
        if (this.powers.length == 0 && other.powers.length == 0) return true;
		if (this.powersBeginIndex != other.powersBeginIndex) return false;
		if (this.powers.length != other.powers.length) return false;
		for (int i = 0; i < powers.length; i++) {
			if (powers[i] != other.powers[i]) return false;
		}
		return true;
	}
	
	public NPolyTerm derivative(int varIndex) {
		double power = getPower(varIndex);
		double newCoef = coefficient*power;
		if (newCoef == 0) return ZERO;
		
		int index = varIndex - powersBeginIndex;
		double[] newPowers = new double[powers.length];
		System.arraycopy(powers, 0, newPowers, 0, powers.length);
		double newPower = (newPowers[index] -= 1);

		if (newPower == 0) {
			int minIndex = getMinIndex(newPowers);
			int maxIndex = getMaxIndex(newPowers);
			int len = maxIndex - minIndex + 1;
			if (len < 0) len = 0;
			double[] shrink = new double[len];
			System.arraycopy(newPowers, minIndex, shrink, 0, len);
			return new NPolyTerm(newCoef, powersBeginIndex+minIndex, shrink);
		}

		return new NPolyTerm(newCoef, powersBeginIndex, newPowers);
	}
	
	private int getMinIndex(double[] arr) {
		int min = arr.length;
		for (int i = arr.length-1; i >= 0; i--) {
			if (arr[i] != 0) {
				min = i;
				break;
			}
		}
		return min;
	}

	private int getMaxIndex(double[] arr) {
		int max = -1;
		for (int i = 0; i < arr.length; i++) {
			if (arr[i] != 0) {
				max = i;
				break;
			}
		}
		return max;
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		
		boolean appendTimes = false;
		if (this.coefficient != 1.0) {
			appendTimes = true;
			if (coefficient == Math.floor(coefficient)) {
				sb.append((int) coefficient);
			} else {
				sb.append(coefficient);
			}
		}
		
		for (int i = 0; i < powers.length; i++) {
			if (powers[i] == 0) continue;
			
			if (appendTimes) sb.append('*');
			appendTimes = true;
			
			sb.append("x[").append(i+powersBeginIndex).append(']');
			if (powers[i] == 1) continue;
			
			sb.append('^');
			if (powers[i] == Math.floor(powers[i])) {
				sb.append((int) powers[i]);
			} else {
				sb.append(powers[i]);
			}
		}
		return sb.toString(); 
	}

	public Scalar value(Vector<Scalar> x) {
		double product = coefficient;
		for (int i = 0; i < x.getSize(); i++) {
			double value = x.get(i).value();
			int ind = i+1-powersBeginIndex;
			if (ind < 0 || ind >= powers.length) continue;
			double power = powers[ind];
			product *= Math.pow(value, power);
		}
		return new Scalar(product);
	}

	public NPolyTerm negate() {
        return new NPolyTerm(-coefficient, powersBeginIndex, powers);
    }

    public NPolyTerm reciprocal() {
        double[] recip = new double[powers.length];
        for (int i = 0; i < powers.length; i++) {
            recip[i] = -powers[i];
        }
        return new NPolyTerm(1.0/coefficient, powersBeginIndex, powers);
    }

	public Value add(Value other) {
        if (other instanceof NPolyTerm) return this.add((NPolyTerm) other);
        if (other instanceof Scalar) return this.add((Scalar) other);
        try {
            return other.add(this);
        } catch (StackOverflowError e) {
            throw new RuntimeException("Addition of "+this.getClass().getName()
                    +" and "+other.getClass().getName()+" is not implemented.");
        }
	}
	
	public NFunction add(NPolyTerm other) {
		if (this.canCombine(other)) {
			return new NPolyTerm(this.coefficient + other.coefficient, this.powersBeginIndex, this.powers);
		}
		
		List<NPolyTerm> parts = new ArrayList<>();
		parts.add(this);
		parts.add(other);
		return new NPolynomial(parts);
	}

	public NFunction add(Scalar other) {
        return this.add(new NPolyTerm(other.value()));
    }
	
	public Value multiply(Value other) {
        if (other instanceof NPolyTerm) return this.multiply((NPolyTerm) other);
        if (other instanceof Scalar) return this.multiply((Scalar) other);
        try {
            return other.multiply(this);
        } catch (StackOverflowError e) {
            throw new RuntimeException("Multiplication of "+this.getClass().getName()
                    +" and "+other.getClass().getName()+" is not implemented.");
        }
	}
	
	public NPolyTerm multiply(NPolyTerm other) {
		int minA = this.powersBeginIndex;
		int minB = other.powersBeginIndex;
		int maxA = this.powersBeginIndex+this.powers.length;
		int maxB = other.powersBeginIndex+other.powers.length;
		int min = Math.min(minA, minB);
		int max = Math.max(maxA, maxB);
		
		double[] output = new double[max-min+1];
		int newMin = Integer.MAX_VALUE, newMax = -1;

		double coef = this.coefficient * other.coefficient;
		if (coef == 0) return ZERO;
		
		for (int i = min; i < max; i++) {
			int indexA = i-this.powersBeginIndex;
			int indexB = i-other.powersBeginIndex;
			output[i-min] = 0;
			if (indexA >= 0 && indexA < this.powers.length)  output[i-min] += this.powers[indexA];
			if (indexB >= 0 && indexB < other.powers.length) output[i-min] += other.powers[indexB];
			if (output[i-min] != 0) {
				if (i < newMin) newMin = i;
				if (i > newMax) newMax = i;
			}
		}

		if (newMax-newMin < 0) return new NPolyTerm(coef);
		double[] trunc = new double[newMax-newMin+1];
		System.arraycopy(output, newMin-min, trunc, 0, newMax-newMin+1);
		
		return new NPolyTerm(coef, newMin, trunc);
	}

    public NPolyTerm multiply(Scalar other) {
        return new NPolyTerm(coefficient*other.value(), powersBeginIndex, powers);
    }

	public NFunctionVector gradient(int maxIndex) {
		NPolyTerm[] values = new NPolyTerm[maxIndex];
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
