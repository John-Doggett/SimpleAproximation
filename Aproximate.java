import java.util.ArrayList;

public class Aproximate {
	private static final complexNumber negative = new complexNumber(-1.0, 0);

	public static double[][] aproximate(Function func, double lowerBound, double upperBound, int numPoints, int polynomialPower) {
		double increment = (upperBound - lowerBound) / ((double) (numPoints - 1));
		complexNumber[][] A = new complexNumber[numPoints][polynomialPower];
		for (int a = 0; a < A.length; a++) {
			A[a][0] = new complexNumber(1, 0);
		}
		for (int a = 0; a < A.length; a++) {
			A[a][1] = new complexNumber(lowerBound + a * increment, 0);
		}
		for (int a = 0; a < A.length; a++) {
			for (int b = 2; b < A[0].length; b++) {
				A[a][b] = new complexNumber(A[a][1].pow(b));
			}
		}
		// A setup

		complexNumber[][] B = new complexNumber[A.length][1];
		for (int a = 0; a < B.length; a++) {
			double[] temp = func.evaulate(a * increment + lowerBound, 0);
			B[a][0] = new complexNumber(temp[0], temp[1]);
		}
		// B setup

		complexNumber[][] AT = transpose(A);
		complexNumber[][] C = matrixMultiply(AT, A);
		complexNumber[][] D = matrixMultiply(AT, B);
		complexNumber[][] E = concat(C, D);
		complexNumber[][] FINAL = RREF(E);
		double[][] output = new double[FINAL.length][2];
		for(int a = 0; a < FINAL.length; a++) {
			complexNumber temp = FINAL[a][FINAL[a].length-1];
			output[a][0] = temp.realConstant;
			output[a][1] = temp.imagineryConstant;
		}
		return output;
	}

	private static complexNumber[][] RREF(complexNumber[][] inputMatrix) {
		complexNumber[][] input = new complexNumber[inputMatrix.length][inputMatrix[0].length];
		for (int a = 0; a < inputMatrix.length; a++) {
			for (int b = 0; b < inputMatrix[0].length; b++) {
				input[a][b] = new complexNumber(inputMatrix[a][b]);
			}
		}

		ArrayList<Integer> pivotPositions = new ArrayList<>();
		int length = input[0].length;
		for (int a = 0; a < input.length; a++) {
			if (input[a].length != length) {
				throw new IllegalArgumentException("Matrix is not mxn!");
			}
		}

		// REF Forward Operation
		int rowAdjustor = 0;
		for (int a = 0; a < input[0].length && a - rowAdjustor < input.length; a++) {

			int largest = a - rowAdjustor;
			for (int b = a - rowAdjustor; b < input.length; b++) {
				if (input[b][a].abs() > input[largest][a].abs()) {
					largest = b;
				}
			}
			swapRow(input, a - rowAdjustor, largest);
			if (input[a - rowAdjustor][a].abs() == 0.0) {
				rowAdjustor++;
			} else {
				pivotPositions.add(a);

				for (int b = a + 1 - rowAdjustor; b < input.length; b++) {
					multiplyAddRow(input, a - rowAdjustor, b,
							negative.multiply(input[b][a].divide(input[a - rowAdjustor][a])));
				}

			}
		}

		// BACKWARD SUBSTITUTION
		// We know where the pivot positions are, which we can use to ignore critical
		// double error.
		// We are in REF.
		for (int a = pivotPositions.size(); a < input.length; a++) {
			for (int b = 0; b < input[a].length; b++) {
				if (a != pivotPositions.size() - 1 && b != input[0].length - 1
						&& pivotPositions.get(pivotPositions.size() - 1) == input[0].length - 1) {
					input[a][b] = new complexNumber(0, 0);
				}
			}
		}

		for (int a = pivotPositions.size() - 1; a >= 0; a--) {
			for (int b = pivotPositions.get(a) - 1; b >= 0; b--) { // Set left of column = 0
				input[a][b] = new complexNumber(0, 0);
			}

			// Divide row by pivot position
			divideRow(input, a, input[a][pivotPositions.get(a)]);
			// multiplyAddRow to each row above
			for (int b = a - 1; b >= 0; b--) {
				multiplyAddRow(input, a, b,
						input[b][pivotPositions.get(a)].divide(input[a][pivotPositions.get(a)]).multiply(negative));
			}
		}

		return input;
	}

	private static void swapRow(complexNumber[][] input, int rowA, int rowB) {
		complexNumber[] temp = input[rowA];
		input[rowA] = input[rowB];
		input[rowB] = temp;
	}

	private static void multiplyAddRow(complexNumber[][] input, int giver, int receiver, complexNumber multible) {
		complexNumber[] temp = new complexNumber[input[receiver].length];
		for (int a = 0; a < temp.length; a++) {
			temp[a] = input[receiver][a].add(input[giver][a].multiply(multible));
		}
		input[receiver] = temp;
	}

	private static void divideRow(complexNumber[][] input, int row, complexNumber divder) {
		for (int a = 0; a < input[row].length; a++) {
			input[row][a] = input[row][a].divide(divder);
		}
	}

	private static complexNumber[][] transpose(complexNumber[][] input) {
		complexNumber[][] output = new complexNumber[input[0].length][input.length];
		for (int a = 0; a < input.length; a++) {
			for (int b = 0; b < input[a].length; b++) {
				output[b][a] = new complexNumber(input[a][b]);
			}
		}
		return output;
	}

	private static complexNumber[][] concat(complexNumber[][] A, complexNumber[][] B) {
		if (A.length != B.length) {
			throw new IllegalArgumentException("Matricies do not have some amount of rows!");
		}
		complexNumber[][] output = new complexNumber[A.length][A[0].length + B[0].length];
		for (int a = 0; a < output.length; a++) {
			for (int b = 0; b < A[0].length; b++) {
				output[a][b] = new complexNumber(A[a][b]);
			}
		}
		for (int a = 0; a < output.length; a++) {
			for (int b = A[0].length; b < output[0].length; b++) {
				output[a][b] = new complexNumber(B[a][b - A[0].length]);
			}
		}
		return output;
	}

	private static complexNumber[][] matrixMultiply(complexNumber[][] matrixOne, complexNumber[][] matrixTwo) {
		complexNumber[][] output = new complexNumber[matrixOne.length][matrixTwo[0].length];
		for (int a = 0; a < output.length; a++) {
			for (int b = 0; b < output[a].length; b++) {
				output[a][b] = new complexNumber(0, 0);
			}
		}
		for (int a = 0; a < matrixTwo[0].length; a++) {
			for (int b = 0; b < matrixOne.length; b++) {
				for (int c = 0; c < matrixTwo.length; c++) {
					output[b][a] = output[b][a].add(matrixOne[b][c].multiply(matrixTwo[c][a]));
				}
			}
		}
		return output;
	}

	private static class complexNumber {
		private double realConstant;
		private double imagineryConstant;

		private complexNumber(double realConstant, double imagineryConstant) {
			super();
			this.realConstant = realConstant;
			this.imagineryConstant = imagineryConstant;
		}

		private complexNumber(complexNumber other) {
			super();
			this.realConstant = other.realConstant;
			this.imagineryConstant = other.imagineryConstant;
		}

		private complexNumber add(complexNumber other) {
			return new complexNumber(realConstant + other.realConstant, imagineryConstant + other.imagineryConstant);
		}

		private complexNumber subtract(complexNumber other) {
			return new complexNumber(realConstant - other.realConstant, imagineryConstant - other.imagineryConstant);
		}

		private complexNumber multiply(complexNumber other) {
			return new complexNumber(realConstant * other.realConstant - imagineryConstant * other.imagineryConstant,
					realConstant * other.imagineryConstant + imagineryConstant * other.realConstant);
		}

		private complexNumber conjugate() {
			return new complexNumber(realConstant, -imagineryConstant);
		}

		private complexNumber divide(complexNumber other) {
			complexNumber conjugate = other.conjugate();
			complexNumber nominator = this.multiply(conjugate);
			double denominator = other.multiply(conjugate).realConstant;
			return new complexNumber(nominator.realConstant / denominator, nominator.imagineryConstant / denominator);
		}

		private complexNumber pow(int power) {
			if (power == 0) {
				return new complexNumber(1, 0);
			} else {
				complexNumber temp = this;
				while (power - 1 > 0) {
					temp = temp.multiply(this);
					power--;
				}
				return temp;
			}
		}

		private double abs() {
			return Math.hypot(realConstant, imagineryConstant);
		}

		@Override
		public String toString() {
			return realConstant + "";
		}
	}
}
