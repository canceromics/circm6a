package output;

public class PeakOut extends Bed12{

	private Integer confidence;
	private double fdr;
	private double propor;

	public PeakOut(String chr, int start, int end, String name, double score, char strand, int t_start, int t_end,
			int rgb, int count, int[] sizes, int[] starts, Integer confidence, double fdr, double propor) {
		super(chr, start, end, name, score, strand, t_start, t_end, rgb, count, sizes, starts);
		this.confidence = confidence;
		this.fdr = fdr;
		this.propor = propor;
	}
	
	public void setPeakName() {
		StringBuilder sb = new StringBuilder();
		if (getName() != null) {
			sb.append(getName());
			sb.append(':');
		}
		sb.append(getChr());
		for (int i = 0; i < getCount(); i++) {
			sb.append(':');
			sb.append(getStarts()[i]);
			sb.append('-');
			sb.append(getStarts()[i] + getSizes()[i] - 1);
		}
	}
	
	public void setFdr(double fdr) {
		this.fdr = fdr;
	}
	
	public PeakOut confidencePeak(Integer confidence) {
		return new PeakOut(getChr(), getStart(), getEnd(), getName(), getScore(), getStrand(), getT_start(), getT_end(), getRgb(),
				getCount(), getSizes(), getStarts(), confidence, fdr, propor);
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(super.toString());
		if (confidence != null) {
			sb.append(Bed3.getSep());
			if (confidence > 0) {
				sb.append("high");
			}
			else if (confidence == -1) {
				sb.append("low");
			}
			else if (confidence < -1){
				sb.append("superlow");
			}
			else {
				sb.append("moderate");
			}
			sb.append(" confidence");
		}
		sb.append(Bed3.getSep());
		sb.append(fdr);
		sb.append(Bed3.getSep());
		sb.append(propor);
		return sb.toString();
	}
	
	public static String getHeader() {
		return Bed12.HEADER + "\tConfidence\tFDR\tProportion";
	}
	
	public static String getLinearHeader() {
		return Bed12.HEADER + "\tFDR\tProportion";
	}
	
	public int getConfidence() {
		return confidence;
	}
}
