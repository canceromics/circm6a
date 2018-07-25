package main;

import java.util.ArrayList;
import java.util.Collections;

public class Transcript {

	private Gene gene=null;
	private String id=null;
	private ArrayList<ExonInfo> exons=null;
	private int start=0;
	private int end=0;
	private int base_sum=0;
	private char strand='+';
	private boolean circ_flag=false;
	
	public Transcript(){
		this.exons = new ArrayList<>();
	}
	
	public Transcript(Gene gene, String id, ArrayList<ExonInfo> exons, int start, int end, char strand,int base_sum, boolean circ_flag) {
		super();
		this.gene = gene;
		this.id = id;
		this.exons = exons;
		this.start = start;
		this.end = end;
		this.base_sum = base_sum;
		this.circ_flag = circ_flag;
	}

	public Gene getGene() {
		return gene;
	}
	
	public void setGene(Gene gene) {
		this.gene = gene;
	}
	
	public String getId() {
		return id;
	}
	
	public void setId(String id) {
		this.id = id;
	}
	
	public ArrayList<ExonInfo> getExons() {
		return exons;
	}
	
	public ExonInfo getExon(int index) {
		return exons.get(index);
	}
	
	public void addExon(ExonInfo exon) {
		this.exons.add(exon);
	}
	
	public void setExons(ArrayList<ExonInfo> exons) {
		this.exons = exons;
	} 
	
	public int getStart() {
		return start;
	}

	public void setStart(int start) {
		this.start = start;
	}

	public int getEnd() {
		return end;
	}

	public void setEnd(int end) {
		this.end = end;
	}

	public char getStrand() {
		return strand;
	}

	public void setStrand(char strand) {
		this.strand = strand;
	}
	
	public int getBase_sum() {
		return base_sum;
	}

	public void setBase_sum(int base_sum) {
		this.base_sum = base_sum;
	}

	public boolean isCirc_flag() {
		return circ_flag;
	}

	public void setCirc_flag(boolean circ_flag) {
		this.circ_flag = circ_flag;
	}

	/**
	 * sort exon regions according to their start position
	 * @param less_front means sort from less to greater while true
	 */
	public void sortExons(boolean less_front) {
		this.quickSortExons(0, this.exons.size()- 1);
		if (!less_front) {
			Collections.reverse(this.exons);
		}
		for (int i = 1; i < this.exons.size(); ++i) {
			this.exons.get(i - 1).setNext(this.exons.get(i));
		}
	}
	
	/**
	 * quick sort funtion
	 * @param start unsorted start
	 * @param end unsorted end
	 */
	private void quickSortExons(int start, int end) {
		if (end - start <= 8) {
			this.insertSortExons(start, end);
			return;
		}
		int left = start;
		int right = end;
		int middle = (left + right) >> 1;
		int key = this.exons.get(middle).getStart_position();
		
		while(left < right) {
			while (this.exons.get(left).getStart_position() <= key) {
				left++;
			}
			while (this.exons.get(right).getStart_position() >= key) {
				right--;
			}
			if (left < right) {
				this.exons.set(right, this.exons.set(left, this.exons.get(right)));
			}
			else if (left < (start+end) >> 1) {
				this.exons.set(left, this.exons.set(middle, this.exons.get(left)));
				right = left;
			}
			else if (right > middle){
				this.exons.set(right, this.exons.set(middle, this.exons.get(right)));
				left = right;
			}
		}
		
		this.quickSortExons(start, left - 1);
		this.quickSortExons(right + 1, end);
	}
	
	/**
	 * insert sort function
	 * @param start unsorted start
	 * @param end unsorted end
	 */
	private void insertSortExons(int start, int end) {
		for (int i=start + 1; i <= end; ++i) {
			int key = this.exons.get(i).getStart_position();
			for (int j=start; j < i; ++j) {
				if (key < this.exons.get(j).getStart_position()) {
					this.exons.add(j, this.exons.get(i));
					this.exons.remove(i + 1);
					break;
				}
			}
		}
	}
	
}
