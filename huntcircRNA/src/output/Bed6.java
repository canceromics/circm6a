package output;

public class Bed6 extends Bed3{

	public final static String HEADER = Bed3.HEADER + "\tname\tscore\tstrand";
	private String name = null;
	private double score = 0.0;
	private char strand = '.';
	
	public Bed6(String chr, int start, int end, String name, double score, char strand) {
		super(chr, start, end);
		this.name = name;
		this.score = score;
		this.strand = strand;
	}

	public String getName() {
		return name;
	}

	public double getScore() {
		return score;
	}
	
	public void setScore(double score) {
		this.score = score;
	}
	
	public char getStrand() {
		return strand;
	}

	public void setName(String name) {
		this.name = name;
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(super.toString());
		sb.append(Bed3.getSep());
		sb.append(name);
		sb.append(Bed3.getSep());
		sb.append(score);
		sb.append(Bed3.getSep());
		sb.append(strand);
		return sb.toString();
	}
	
}
