package generic;

import java.lang.reflect.Array;
import java.util.function.Function;

/**
 * Created by Gregary on 11/9/2016.
 */
public class Matrix<T extends Value> {

    final Class<T> clazz;
    final int rows, cols;
    final T[] values;
    final boolean transpose;

    public static Matrix<Scalar> identity(int size) {
        Scalar[] ident = new Scalar[size*size];
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (i == j) ident[i * size + j] = new Scalar(1);
                else ident[i * size + j] = new Scalar(0);
            }
        }
        return new Matrix<>(size, size, ident, false);
    }

    @SuppressWarnings("unchecked")
    Matrix(int rows, int cols, T[] values, boolean transpose) {
        if (values.length < rows*cols) throw new RuntimeException("Too few values");
        if (values.length > rows*cols) throw new RuntimeException("Too many values");
        this.rows = rows;
        this.cols = cols;
        this.values = values;
        this.transpose = transpose;
        this.clazz = (Class<T>) values[0].getClass();
    }

    // Default constructor, copy array for security
    public Matrix(T[][] values) {
        this(values, false);
    }

    @SuppressWarnings("unchecked")
    public Matrix(T[][] values, boolean transpose) {
        this.rows = values.length;
        if (rows == 0) throw new RuntimeException("Zero dimension");
        this.cols = values[0].length;
        if (cols == 0) throw new RuntimeException("Zero dimension");

       this.clazz = (Class<T>) values[0][0].getClass();
        this.values = (T[]) Array.newInstance(clazz, rows*cols);

        for (int i = 0; i < rows; i++) {
            if (values[i].length != cols) throw new RuntimeException("Non-square matrix");
            System.arraycopy(values[i], 0, this.values, i*cols, cols);
        }
        this.transpose = transpose;
    }

    public Vector<T> asVector() {
        if (this instanceof Vector) return (Vector<T>) this;
        if (this.getCols() == 1) {
            return new Vector<>(this.values, false);
        } else if (this.getRows() == 1) {
            return new Vector<>(this.values, true);
        } else {
            throw new RuntimeException("Matrix is not a vector");
        }
    }

    public Matrix<Scalar> toScalars() {
        return this.apply(Scalar.class, (val) -> {
            if (val instanceof Scalar) return (Scalar) val;
            else throw new RuntimeException("Value is not scalar compatible: \""+val.toString()+"\"");
        });
    }

    public Matrix<Value> toValues() {
        return this.apply(Value.class, (val) -> (Value) val);
    }

    public Matrix<T> transpose() {
        return new Matrix<>(rows, cols, values, !transpose);
    }

    public int getSize() {
        return getRows() * getCols();
    }

    public int getRows() {
        if (transpose) return cols;
        else return rows;
    }

    public int getCols() {
        if (transpose) return rows;
        else return cols;
    }

    public T get(int row, int col) {
        if (transpose) { int swap = row; row = col; col = swap; }
        return values[row*cols+col];
    }

    public T get(int index) {
        int row = index / getCols();
        int col = index - row*getCols();
        return get(row, col);
    }

    @SuppressWarnings("unchecked")
    public Vector<T> getRow(int row) {
        T[] vector = (T[]) Array.newInstance(clazz, getCols());
        for (int i = 0; i < getCols(); i++) {
            vector[i] = get(row, i);
        }
        return new Vector<>(vector, true); // row vector
    }

    @SuppressWarnings("unchecked")
    public Vector<T> getCol(int col) {
        T[] vector = (T[]) Array.newInstance(clazz, getRows());
        for (int i = 0; i < getRows(); i++) {
            vector[i] = get(i, col);
        }
        return new Vector<>(vector, false); // column vector
    }

    @SuppressWarnings("unchecked")
    public Matrix<T> reduce(int row, int col) {
        int rows = getRows()-1;
        int cols = getCols()-1;
        if (rows <= 0 || cols <= 0) throw new RuntimeException("Can't reduce matrix, too small");

        T[] reduced = (T[]) Array.newInstance(clazz, rows*cols);
        int reducedR = 0;
        for (int r = 0; r < getRows(); r++) {
            if (r == row) continue;
            int reducedC = 0;
            for (int c = 0; c < getCols(); c++) {
                if (c == col) continue;
                reduced[reducedR*cols+reducedC] = values[r*getCols()+c];
                reducedC++;
            }
            reducedR++;
        }

        return new Matrix<T>(rows, cols, reduced, false);
    }

    public Value cofactor(int row, int col) {
        double sign = 1;
        if ((row+col) % 2 == 1) sign = -1;
        return reduce(row, col).determinant().multiply(new Scalar(sign));
    }

    public Matrix<Value> adjoint() {
        if (getRows() != getCols()) throw new RuntimeException("Can't computer adjoint of non-square matrix");
        Value[] adjoint = new Value[rows*cols];
        for (int r = 0; r < getRows(); r++) {
            for (int c = 0; c < getCols(); c++) {
                adjoint[r*getCols()+c] = cofactor(r, c);
            }
        }
        return new Matrix<>(getRows(), getCols(), adjoint, true);
    }

    public Value determinant() {
        if (getRows() != getCols()) throw new RuntimeException("Can't computer determinant of non-square matrix");
        if (getRows() == 1) return get(0,0);
        if (getRows() == 2) return get(0,0).multiply(get(1,1)).add(get(1,0).multiply(get(0,1)).negate());
        Value det = new Scalar(0);
        for (int i = 0; i < getRows(); i++) {
            Value cofactor = cofactor(0, i);
            det = det.add(get(0,i).multiply(cofactor));
        }
        return det;
    }

    @SuppressWarnings("unchecked")
    public Matrix<Value> inverse() {
        if (getRows() != getCols()) throw new RuntimeException("Can't computer inverse of non-square matrix");
        Value det = determinant();
        if (((Scalar) det).value() == 0) return null; // Non-invertible
        return adjoint().multiply(det.reciprocal());
    }

    @SuppressWarnings("unchecked")
    public <R extends Value> Matrix<R> apply(Class<R> clazz, Function<T, R> function) {
        R[] output = (R[]) Array.newInstance(clazz, rows*cols);
        for (int i = 0; i < rows*cols; i++) {
            output[i] = function.apply(values[i]);
        }
        return new Matrix<>(rows, cols, output, transpose);
    }

    public Matrix<Value> negate() {
        return this.multiply(new Scalar(-1));
    }

    public Matrix<Value> reciprocal() {
        Matrix<Value> inverse = inverse();
        if (inverse == null) throw new RuntimeException("Non-invertible matrix:\n"+this);
        return inverse;
    }

    public Matrix<Value> add(Matrix<Value> other) {
        if (other.getRows() != getRows() || other.getCols() != getCols()) {
            throw new RuntimeException("Cannot add matrices of different sizes");
        }
        Value[] values = new Value[getSize()];
        for (int r = 0; r < getRows(); r++) {
            for (int c = 0; c < getCols(); c++) {
                values[r*getCols()+c] = this.get(r,c).add(other.get(r,c));
            }
        }
        return new Matrix<>(getRows(), getCols(), values, false);
    }

    public Matrix<Value> multiply(Matrix<Value> other) {
        if (this.getCols() != other.getRows()) throw new RuntimeException("Matrix multiplication bad dimensions");
        Value[] vals = new Value[this.getRows() * other.getCols()];
        for (int r = 0; r < this.getRows(); r++) {
            for (int c = 0; c < other.getCols(); c++) {
                Vector<T> aRow = getRow(r);
                Vector<Value> bCol = other.getCol(c);
                vals[r * other.getCols() + c] = aRow.dot(bCol);
            }
        }
        return new Matrix<>(this.getRows(), other.getCols(), vals, false);
    }

    public Matrix<Value> multiply(Scalar other) {
        return this.apply(Value.class, (elem)->elem.multiply(other));
    }

    public Matrix<Value> multiply(Value other) {
        return this.apply(Value.class, (elem)->elem.multiply(other));
    }

    private static final int PADDING = 2;

    public String toString() {
        String[][] valStr = new String[getRows()][getCols()];
        for (int y = 0; y < getRows(); y++) {
            for (int x = 0; x < getCols(); x++) {
                valStr[y][x] = String.valueOf(get(y, x));
            }
        }

        int totalWidth = 0;
        for (int x = 0; x < getCols(); x++) {
            int width = 0;
            for (int y = 0; y < getCols(); y++) {
                int len = valStr[y][x].length();
                if (len > width) width = len;
            }
            width += PADDING;
            for (int y = 0; y < getCols(); y++) {
                int len = valStr[y][x].length();
                valStr[y][x] = valStr[y][x] + spaces(width-len);
            }
            totalWidth += width;
        }

        String pad = spaces(PADDING);
        int rowLength = totalWidth+PADDING+3;

        StringBuilder sb = new StringBuilder(rowLength * values.length);
        for (int y = 0; y < getRows(); y++) {
            sb.append('[').append(pad);
            for (int x = 0; x < getCols(); x++) {
                sb.append(valStr[y][x]);
            }
            sb.append(']').append('\n');
        }

        return sb.toString();
    }

    private String spaces(int spaces) {
        StringBuilder sb = new StringBuilder(spaces);
        for (int i = 0; i < spaces; i++) sb.append(' ');
        return sb.toString();
    }

}
