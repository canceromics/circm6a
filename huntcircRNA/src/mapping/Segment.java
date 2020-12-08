package mapping;

import genome.IntRegion;
import main.InParam;

public class Segment extends IntRegion{
	
	private boolean circ_flag = false;
	private boolean negative = false;
	private Segment same_position = null;
	private int quali = 0;
	
	
	public Segment(int start, int end, boolean negative, int quali) {
		super(start, end);
		this.negative = negative;
		this.quali = quali;
	}

	public void expandSegment(int start, int end) {
		this.resetStartAndEnd(Math.min(start, getStart()), Math.max(end, getEnd()));
	}
	
	public void addAnotherSeg(Segment seg) {
		if (seg != null) {
			seg.same_position = this.same_position;
			this.same_position = seg;
		}
	}
	
	public Segment getAnotherSeg() {
		return same_position;
	}
	
	public boolean isNegative() {
		return negative;
	}
	
	public boolean isOnceCirc() {
		return circ_flag ? true : !(circ_flag = true);
	}
	
	public int getQuality() {
		return quali;
	}
	
	public boolean isEnoughQuality() {
		return quali >= InParam.getParams().getMapQuality();
	}
}
