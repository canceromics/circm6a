package output;

public class Bed12 extends Bed6 {

	public final static String HEADER = Bed6.HEADER + "\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts";
	
	private int t_start;
	private int t_end;
	private int rgb;
	private int count;
	private int[] sizes;
	private int[] starts;
	
	public Bed12(String chr, int start, int end, String name, double score, char strand, int t_start, int t_end,
			int rgb, int count, int[] sizes, int[] starts) {
		super(chr, start, end, name, score, strand);
		this.t_start = t_start;
		this.t_end = t_end;
		this.rgb = rgb;
		this.count = count;
		this.sizes = sizes;
		this.starts = starts;
	}
	
	public int getT_start() {
		return t_start;
	}

	public int getT_end() {
		return t_end;
	}

	public int getRgb() {
		return rgb;
	}

	public int getCount() {
		return count;
	}

	public int[] getSizes() {
		return sizes;
	}

	public int[] getStarts() {
		return starts;
	}

	public boolean checkBlock() {
		return count > 0 && sizes != null && sizes.length >= count && starts != null && starts.length >= count;
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(super.toString());
		if (checkBlock()) {
			sb.append(Bed3.getSep());
			sb.append(t_start);
			sb.append(Bed3.getSep());
			sb.append(t_end);
			sb.append(Bed3.getSep());
			sb.append(rgb);
			sb.append(Bed3.getSep());
			sb.append(count);
			sb.append(Bed3.getSep());
			sb.append(sizes[0]);
			for (int i = 1; i < count; ++i) {
				sb.append(',');
				sb.append(sizes[i]);
			}
			sb.append(Bed3.getSep());
			sb.append(starts[0]);
			for (int i = 1; i < count; ++i) {
				sb.append(',');
				sb.append(starts[i]);
			}
		}
		return sb.toString();
	}
}
