package genome;

import java.util.ArrayList;
import java.util.Collections;

import main.Method;

public class Gene extends IntRegion{
	
	private String chr_symbol = null;
	private String gene_id = null;
	private String gene_symbol = null;
	private char strand = '.';
	private ArrayList<Transcript> scripts = null;
	private boolean sorted = false;
	
	public Gene(String chr_symbol, String gene_id, String gene_symbol, char strand, int start, int end, ArrayList<Transcript> scripts) {
		super(start, end);
		this.chr_symbol = chr_symbol;
		this.gene_symbol = gene_symbol;
		this.gene_id = gene_id;
		this.strand = strand;
		this.scripts = scripts;
	}

	public String getChr_symbol() {
		return chr_symbol;
	}
	
	public String getGene_id() {
		return gene_id;
	}

	public void appendGeneID(String id) {
		gene_id = gene_id + "," + id;
	}
	
	public String getGene_symbol() {
		return gene_symbol;
	}
	
	public void appendGeneSymbol(String symbol) {
		gene_symbol = gene_symbol + "," + symbol;
	}

	public char getStrand() {
		return strand;
	}

	public void checkStrand(char strand) {
		this.strand = strand == this.strand ? this.strand : '.';
	}
	
	public void mergeGene(Gene gene) {
		appendGeneID(gene.gene_id);
		appendGeneSymbol(gene.gene_symbol);
		checkStrand(gene.strand);
		scripts.addAll(gene.scripts);
	}
	
	public int fixToExonStart(int target, int dev) {
		int out = -128;
		for (Transcript script : scripts) {
			out = Method.nearestOne(target, script.fixToExonStart(target, dev), out);
			if (out == target) {
				break;
			}
		}
		return out;
	}
	
	public int fixToExonEnd(int target, int dev) {
		int out = -128;
		for (Transcript script : scripts) {
			out = Method.nearestOne(target, script.fixToExonEnd(target, dev), out);
			if (out == target) {
				break;
			}
		}
		return out;
	}
	
	public boolean checkPos() {
		boolean out = true;
		if (!sorted && scripts != null && scripts.size() > 0) {
			for (int i = 0; i < scripts.size(); i++) {
				out &= scripts.get(i).checkPos()
						&& getStart() <= scripts.get(i).getStart()
						&& getEnd() >= scripts.get(i).getEnd();
			}
			Collections.sort(scripts);
			resetStartAndEnd(getStart() > scripts.get(0).getStart() ? scripts.get(0).getStart() : getStart(), 
					getEnd() < scripts.get(scripts.size() - 1).getEnd() ? scripts.get(scripts.size() - 1).getEnd() : getEnd());
		}
		sorted = true;
		return out;
	}

	public ArrayList<Transcript> getScripts() {
		return scripts;
	}
	
}
