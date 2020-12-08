package genome;

import java.util.ArrayList;
import java.util.Collections;

import main.Method;

public class Transcript extends IntRegion{

	private Gene gene = null;
	private String script_id = null;
	private ArrayList<Exon> exons = null;
	
	public Transcript(Gene gene, String script_id, int start, int end, ArrayList<Exon> exons) {
		super(start, end);
		this.gene = gene;
		this.script_id = script_id;
		this.exons = exons;
	}

	public String getChr() {
		return gene==null? null : gene.getChr_symbol();
	}
	
	public Gene getGene() {
		return gene;
	}
	
	public String getScript_id() {
		return script_id;
	}

	public ArrayList<Exon> getExons() {
		return exons;
	}
	
	public int fixToExonStart(int target, int dev) {
		int out = -128;
		for (Exon exon : exons) {
			out = Method.nearestOne(target, exon.fixToExonStart(target, dev), out);
			if (out == target) {
				break;
			}
		}
		return out;
	}
	
	public int fixToExonEnd(int target, int dev) {
		int out = -128;
		for (Exon exon : exons) {
			out = Method.nearestOne(target, exon.fixToExonEnd(target, dev), out);
			if (out == target) {
				break;
			}
		}
		return out;
	}
	
	public int indexOfExons(int pos) {
		int out = 0;
		if (exons != null) {
			for (int i = 0; i < exons.size(); i++) {
				if (pos < exons.get(i).getEnd() && pos >= exons.get(i).getStart()) {
					out += pos - exons.get(i).getStart();
					return out;
				}
				out += exons.get(i).getEnd() - exons.get(i).getStart();
			}
		}
		return -1;
	}

	public boolean checkPos() {
		boolean out = true;
		if (exons != null) {

			Collections.sort(exons);
			for (int i = exons.size() - 1; i >= 0; i--) {
				out &= getStart() <= exons.get(i).getStart()
						&& getEnd() >= exons.get(i).getEnd();
				if (i > 0 && exons.get(i).getStart() <= exons.get(i - 1).getEnd()) {
					exons.get(i - 1).resetStartAndEnd(exons.get(i - 1).getStart(), Math.max(exons.get(i).getEnd(), exons.get(i - 1).getEnd()));
					exons.remove(i);
				}
			}
			resetStartAndEnd(getStart() > exons.get(0).getStart() ? exons.get(0).getStart() : getStart(), 
					getEnd() < exons.get(exons.size() - 1).getEnd() ? exons.get(exons.size() - 1).getEnd() : getEnd());
		}
		return out;
	}
}
