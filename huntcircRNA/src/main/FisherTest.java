package main;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;

public class FisherTest {
	
	public static void main (String[] args) {
		List<Double> l = new ArrayList<>();
		for (String s : args) {
			l.add(Double.parseDouble(s));
		}
		FisherTest ft = new FisherTest(100694394 + 49026762 + 32384);
		System.out.println(ft.calpValue(3, 5, 18821682 - 3, 100694394 -5, 2));
		System.out.println("END");
	}
	
	private double[][] log_factorial = null;
	private double[] small = null;
	
	/**
	 * new FisherTest must be creat with max num possiblely present
	 * @param max_num max num possiblely present when using this test
	 */
	public FisherTest(int max_num) {
        this.log_factorial = new double[1][max_num];
        this.log_factorial[0][0] = 0.0;
        for (int i = 1; i < max_num; ++i) {
            this.log_factorial[0][i] = this.log_factorial[0][i-1] + Math.log(i);
        }
	}
	
	public FisherTest(long ip, long input, int range) {
		this.log_factorial = new double[4][range];
		double fac = 0.0;
		this.log_factorial[0][0] = 0.0;
		for (long i = 1; i <= ip + input; ++i) {
			fac += Math.log(i);
			if (i < range) {
				this.log_factorial[0][(int) i] = fac;
			}
			if (ip - i < range && ip - i >= 0) {
				this.log_factorial[1][(int) (ip - i)] = fac;
			}
			if (input - i < range && input - i >= 0) {
				this.log_factorial[2][(int) (input - i)] = fac;
			}
			if (ip + input - i < range) {
				this.log_factorial[3][(int) (input + ip - i)] = fac;
			}
		}
		small = new double[]{ip, input, ip + input};
	}
	
    /** 
     * Calculate a p-value for Fisher's Exact Test.
     * @param mode 1 means less; 2 means greater; 3 means two-tailed
     * @return the p-value
     */
	public double calpValue(int a, int b, long c, long d, int mode) {
		double p_value = 0.0;
		if	(!this.canTest(a, b, c, d)) {
			System.out.println("Cannot test negative numbers in fisher_exact_test");
			System.out.printf("\tWhen a = %d, b = %d, c = %d, d = %d\n", a, b, c, d);
			return -1;
		}
		
		switch (mode) {
		case 1:
	        while (a >= 0) {
				double p = fisherSub(a, b, c, d);
	            p_value += p;
	            --a;
	            ++b;
	            ++c;
	            --d;
	        }
	        break;
		case 2:
			while (b >= 0) {
				double p = fisherSub(a, b, c, d);
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
			double p = fisherSub(a, b, c, d);
			p_value += p;
			double p_org = p;
	        while (a > 0) {
	        	--a;
	            ++b;
	            ++c;
	            --d;
				p = fisherSub(a, b, c, d);
	            p_value += p;
	        }
	        a = b;
	        b = 0;
	        c = c - a;
	        d = d + a;
	        p = fisherSub(a, b, c, d);

	        while (p < p_org) {
	            if (a == a_org) {
	            	break;
	            }
	            p_value += p;
	            --a;
	            ++b;
	            ++c;
	            --d;
	            p = fisherSub(a, b, c, d);
	        }
	        break;
	    default:
	    	System.out.println("NO such fisher_exact_test mode");
		}
		return effectiveNum(p_value, 6);
	}
	
    private double effectiveNum(double p_value, int n) {
    	if (p_value == 0.0) {
    		return 0.0;
    	}

        double magnitude = Math.pow(10, n - (int) Math.ceil(Math.log10(p_value < 0 ? -p_value: p_value)));
        long shifted = Math.round(p_value * magnitude);
		return shifted / magnitude;
	}

	private double fisherSub(int a, int b, long c, long d) {
        return small == null ? Math.exp(log_factorial[0][a + b] +
        		log_factorial[0][(int) (c + d)] +
        		log_factorial[0][(int) (a + c)] +
        		log_factorial[0][(int) (b + d)] -
        		log_factorial[0][(int) (a + b + c + d)] -
        		log_factorial[0][a] -
        		log_factorial[0][b] -
        		log_factorial[0][(int) c] -
        		log_factorial[0][(int) d])
        		: Math.exp(log_factorial[0][(int) (a + b)] +
        		log_factorial[3][(int) Math.round(small[2] - c - d)] +
        		log_factorial[1][0] +
        		log_factorial[2][0] -
        		log_factorial[3][0] - 
        		log_factorial[0][a] -
        		log_factorial[0][b] -
        		log_factorial[1][(int) Math.round(small[0] - c)] -
        		log_factorial[2][(int) Math.round(small[1] - d)]);
    }
	
	private boolean canTest(long a, long b, long c, long d) {
		return a >= 0 && b >= 0 && c >= 0 && d >= 0;
	}
	
	public static double combinePvalue(Collection<Double> p_values) {
		if (p_values == null || p_values.size() == 0) {
			return Double.NaN;
		}
		double p_value = 0.0;
		switch (InParam.getParams().getCombine_p()) {
		case "ave":
			for (Double d : p_values) {
				p_value += d;
			}
			p_value = p_value / p_values.size();
			break;
			
		case "fm":
			for (Double d : p_values) {
				p_value += d <= 0.0 ? Math.log(Double.MIN_VALUE) : Math.log(d);
			}
			p_value *= -2.0;
			ChiSquaredDistribution csd = new ChiSquaredDistribution(2.0 * p_values.size());
			p_value = 1.0 - csd.cumulativeProbability(p_value);
			break;
			
		default:
			System.err.println("Warning: Unknown method of combining p-value, using average of -log");
		case "log":
			for (Double d : p_values) {
				p_value += d <= 0.0 ? -Math.log(Double.MIN_VALUE) : -Math.log(d);
			}
			p_value = p_value / p_values.size();
			break;
		}
		return p_value;
	}
}
