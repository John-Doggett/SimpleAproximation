
public class Driver implements Function{

	@Override
	public double[] evaulate(double realX, double ImagX) {
		return new double[] {Math.pow(realX, 3) - Math.pow(2, realX),0};
	}
	
	public static void main(String[] args) {
		Function fun = new Driver();
		double[][] polynomial = Aproximate.aproximate(fun,-5,5,100, 20);
		for(int a = 0; a < polynomial.length; a++) {
			System.out.println(polynomial[a][0] + " " + polynomial[a][1] + "I");
		}
		double x=3.0;
		double sum = 0.0;
		for(int a = 0; a< polynomial.length; a++) {
			sum+=polynomial[a][0]*Math.pow(x, a);
		}
		System.out.println(sum);
	}
}
