package genome;

import java.util.ArrayList;
import java.util.Collections;

public abstract class IntRegion implements Comparable<IntRegion>{
	
	private int start;
	private int end;
	
	public IntRegion(int start, int end){
		super();
		this.start = start > end ? end : start;
		this.end = start ^ end ^ this.start;
	}
	
	public boolean setStart(int start) {
		if (start <= end) {
			this.start = start;
			return true;
		}
		return false;
	}
	
	public int getStart() {
		return start;
	}
	
	public boolean setEnd(int end) {
		if (start <= end) {
			this.end = end;
			return true;
		}
		return false;
	}

	
	public int getEnd() {
		return end;
	}
	
	public int getLength() {
		return end - start + 1;
	}
	
	public void resetStartAndEnd(int start, int end) {
		this.start = start > end ? end : start;
		this.end = start ^ end ^ this.start;
	}
	
	public static <T extends IntRegion> ArrayList<T> mergeIntRegion(ArrayList<T> regions) {
		if (regions == null || regions.size() == 0) {
			return regions;
		}
		ArrayList<T> merge = new ArrayList<>();
		Collections.sort(regions);
		merge.add(regions.get(0));
		int last_end = merge.get(0).getEnd();
		for (int i = 1; i < regions.size(); i++) {
			if (regions.get(i).getStart() <= last_end + 1) {
				if (last_end < regions.get(i).getEnd()) {
					merge.get(merge.size() - 1).setEnd(regions.get(i).getEnd());
				}
			}
			else {
				merge.add(regions.get(i));
			}
		}
		return merge;
	}
	
	@Override
	public int compareTo(IntRegion anotherRegion) {
		return start < anotherRegion.start ? -1 : 
			start > anotherRegion.start ? 1 : 
			end < anotherRegion.end ? -1 : 
			end > anotherRegion.end ? 1 : 0;
	}
	
	@Override
	public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o instanceof IntRegion) {
        	IntRegion r = (IntRegion) o;
            return r.start == this.start && r.end == this.end;
        }
        return false;
	}
}
