package main;

public class FisherTest {
	
	private double[] log_factorial = null;

	/**
	 * new FisherTest must be creat with max num possiblely present
	 * @param max_num max num possiblely present when using this test
	 */
	public FisherTest(int max_num) {
        this.log_factorial = new double[max_num+1];
        this.log_factorial[0] = 0.0;
        for (int i = 1; i <= max_num; ++i) {
            this.log_factorial[i] = this.log_factorial[i-1] + Math.log(i);
        }
	}
	
    /** 
     * Calculate a p-value for Fisher's Exact Test.
     * @param mode 1 means less; 2 means greater; 3 means two-tailed
     * @return the p-value
     */
	public double calpValue(int a, int b, int c, int d, int mode) {
		double p_value = 0.0;
		if	(!this.canTest(a, b, c, d)) {
			System.out.println("Cannot test negative numbers in fisher_exact_test");
			System.out.printf("\tWhen a = %d, b = %d, c = %d, d = %d\n", a, b, c, d);
			return -1;
		}
		if (a > d) {
			a=a^d;
			d=a^d;
			a=a^d;
		}
		if (b > c) {
			b=b^c;
			c=b^c;
			b=b^c;
		}
		
		switch (mode) {
		case 1:
	        while (a >= 0) {
				double p = fisherSub(a, b, c, d, this.log_factorial);
	            p_value += p;
	            --a;
	            ++b;
	            ++c;
	            --d;
	        }
	        break;
		case 2:
			while (b >= 0) {
				double p = fisherSub(a, b, c, d, this.log_factorial);
	            p_value += p;
	            ++a;
	            --b;
	            --c;
	            ++d;
			}
			break;
		case 3:
	        if (a * d > b * c) {
	            a=a^b;
	            b=a^b;
	            a=a^b; 
	            c=c^d;
	            d=c^d;
	            c=c^d;
	        }
			int a_org = a;
			double p = fisherSub(a, b, c, d, this.log_factorial);
			p_value += p;
			double p_org = p;
	        while (a > 0) {
	        	--a;
	            ++b;
	            ++c;
	            --d;
				p = fisherSub(a, b, c, d, this.log_factorial);
	            p_value += p;
	        }
	        a = b;
	        b = 0;
	        c = c - a;
	        d = d + a;
	        p = fisherSub(a, b, c, d, this.log_factorial);

	        while (p < p_org) {
	            if (a == a_org) {
	            	break;
	            }
	            p_value += p;
	            --a;
	            ++b;
	            ++c;
	            --d;
	            p = fisherSub(a, b, c, d, this.log_factorial);
	        }
	        break;
	    default:
	    	System.out.println("NO such fisher_exact_test mode");
		}
		return p_value;
	}
	
    private double fisherSub(int a, int b, int c, int d, double[] logFactorial) {
        return Math.exp(logFactorial[a + b] +
                        logFactorial[c + d] +
                        logFactorial[a + c] +
                        logFactorial[b + d] -
                        logFactorial[a + b + c + d] -
                        logFactorial[a] -
                        logFactorial[b] -
                        logFactorial[c] -
                        logFactorial[d]);
    }
	
	private boolean canTest(int a, int b, int c, int d) {
		return a >= 0 && b >= 0 && c >= 0 && d >= 0;
	}
}
