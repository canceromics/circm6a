package output;

import main.InParam;

public class CircOut extends Bed6 {
	
	private int input_reads;
	private int ip_reads;
	private int linear;
	private double circ_ratio;
	
	public CircOut(String chr, int start, int end, String name, double score, char strand, int input_reads,
			int ip_reads, int linear, double circ_ratio) {
		super(chr, start, end, name, score, strand);
		this.input_reads = input_reads;
		this.ip_reads = ip_reads;
		this.linear = linear;
		this.circ_ratio = circ_ratio;
	}
	
	public static String getHeader(boolean ip_flag) {
		String header = "#Chr\tStart\tEnd\tGene Name\tScore\tStrand\tINPUT JunctionReads\tLinearReads\tCircRatio";
		if (InParam.getParams().isOut_detail()) {
			header = header + "\tIDs";
		}
		if (ip_flag) {
			StringBuilder sb = new StringBuilder();
			sb.append(header);
			sb.insert(header.indexOf("INPUT "), "IP JunctionReads\t");
			return sb.toString();
		}
		else {
			return header;
		}
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(super.toString());
		if (InParam.getParams().getIp_file() != null) {
			sb.append(Bed3.getSep());
			sb.append(ip_reads);
		}
		sb.append(Bed3.getSep());
		sb.append(input_reads);
		sb.append(Bed3.getSep());
		sb.append(linear);
		sb.append(Bed3.getSep());
		sb.append(circ_ratio);
		return sb.toString();
	}
}
