package mapping;

public class MissingAlignment {

	private Segment seg = null;
	private boolean front;
	private boolean ip_flag;
	private String seq = null;
	private MissingAlignment same_position = null;
	
	public MissingAlignment(boolean front_missing, boolean ip, String seq, Segment seg) {
		super();
		this.front = front_missing;
		this.seq = seq;
		this.ip_flag = ip;
		this.seg = seg;
	}
	
	public boolean isFrontMissing() {
		return front;
	}
	
	public String getSeq() {
		return seq;
	}
	
	public Segment getSeg() {
		return seg;
	}
	
	public void addAnotherMissingAlignment(MissingAlignment mis) {
		if (mis != null) {
			mis.same_position = this.same_position;
			this.same_position = mis;
		}
	}
	
	public MissingAlignment getAnotherMissingAlignment() {
		return same_position;
	}
	
	public boolean isIP() {
		return ip_flag;
	}
}
