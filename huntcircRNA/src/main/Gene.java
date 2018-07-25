package main;

import java.util.ArrayList;
import java.util.Iterator;

import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;

public class Gene {
	private String gene_symbol=null;
	private String gene_id=null;
	private int start=0;
	private int end=0;
	private char strand='.';
	private ArrayList<Transcript> transcripts=null;
	
	public Gene() {
		this.setTranscripts(new ArrayList<>());
	}
	
	public String getGene_symbol() {
		return gene_symbol;
	}
	
	public void setGene_symbol(String gene_symbol) {
		this.gene_symbol = gene_symbol;
	}
	
	public String getGene_id() {
		return gene_id;
	}

	public void setGene_id(String gene_id) {
		this.gene_id = gene_id;
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

	public ArrayList<Transcript> getTranscripts() {
		return transcripts;
	}

	public void setTranscripts(ArrayList<Transcript> transcripts) {
		this.transcripts = transcripts;
	}
	
	/**
	 * judging whether the junction is near the bound of exon in this gene
	 * @param junc the junction;
	 * @param dev deviation of the judging;
	 * @return true if at least one site near the bound of exons
	 */
	public ArrayList<Integer> isIntron(JuncInfo junc, int dev) {
		ArrayList<Integer> out = new ArrayList<>();
		IntervalTree<ExonInfo> exon_tree = this.bulidExonTree();
		Iterator<Node<ExonInfo>> nodes = exon_tree.overlappers(junc.getSP(), junc.getEP());
		while (nodes.hasNext()) {
			Node<ExonInfo> node = nodes.next();
			ExonInfo exon = node.getValue();
			int start = exon.getStart_position();
			int end = exon.getEnd_position();
			if (junc.getSP() - dev <= exon.getStart_position() && junc.getSP() + dev >= exon.getStart_position()) {
				junc.setSP(exon.getStart_position());
				junc.setFix_start_exon(true);
			}
			if (junc.getEP() - dev <= exon.getEnd_position() && junc.getEP() + dev >= exon.getEnd_position()) {
				junc.setEP(exon.getEnd_position());
				junc.setFix_end_exon(true);
			}
			start = Math.max(start, junc.getSP());
			end = Math.min(end, junc.getEP());
			if (out.size() == 0 || out.get(out.size() -1) < exon.getStart_position()) {
				out.add(start);
				out.add(end);
			}
		}
		return out;
	}
	
	public ArrayList<Integer> scriptExons(JuncInfo junc, int dev, int[] intron) {
		ArrayList<Integer> out = new ArrayList<>();
		for (int i = 0; i < this.transcripts.size(); ++i) {
			Transcript script = this.transcripts.get(i);
			int fixed_start = -1;
			int fixed_end = -1;
			ArrayList<Integer> exons = new ArrayList<>();
			int[] temp_in = {-1, Integer.MAX_VALUE};
			for (int j = 0; j < script.getExons().size(); ++j) {
				ExonInfo exon = script.getExon(j);
				int start = exon.getStart_position();
				int end = exon.getEnd_position();
				if (junc.getSP() - dev <= start && junc.getSP() + dev >= start) {
					fixed_start = start;
				}
				if (junc.getEP() - dev <= end && junc.getEP() + dev >= end) {
					fixed_end = end;
				}
				if (junc.getSP() <= end && junc.getEP() >= start) {
					exons.add(start);
					exons.add(end);
				}
				else if (junc.getSP() > end){
					temp_in[1] = Math.min(temp_in[1], junc.getSP() - 1);
				}
				else {
					temp_in[0] = Math.max(temp_in[0], junc.getEP() + 1);
				}
			}
			if (exons.size() > out.size()) {
				boolean fix_flag = false;
				if (fixed_start >= 0) {
					junc.setSP(fixed_start);
					junc.setFix_start_exon(true);
					fix_flag = true;
				}
				if (fixed_end >= 0) {
					junc.setEP(fixed_end);
					junc.setFix_end_exon(true);
					fix_flag = true;
				}
				if (fix_flag || (!junc.isFix_start_exon() && !junc.isFix_end_exon())){
					out = exons;
					intron[0] = temp_in[0];
					intron[1] = temp_in[1];
				}
			}
			else {
				boolean fix_flag = false;
				if (!junc.isFix_start_exon() && fixed_start >= 0) {
					junc.setSP(fixed_start);
					junc.setFix_start_exon(true);
					fix_flag = true;
				}
				if (!junc.isFix_end_exon() && fixed_end >= 0) {
					junc.setEP(fixed_end);
					junc.setFix_end_exon(true);
					fix_flag = true;
				}
				if (fix_flag) {
					out = exons;
					intron[0] = temp_in[0];
					intron[1] = temp_in[1];
				}
			}
		}
		return out;
	}
	
	/**
	 * build a tree of all exons in this gene
	 * @return tree of all exons
	 */
	public IntervalTree<ExonInfo> bulidExonTree(){
		IntervalTree<ExonInfo> exon_tree = new IntervalTree<>();
		for (int i = 0; i < this.transcripts.size(); ++i) {
			Transcript script = this.transcripts.get(i);
			for (int j = 0; j < script.getExons().size(); ++j) {
				ExonInfo exon = script.getExon(j);
				exon_tree.put(exon.getStart_position(), exon.getEnd_position(), exon);
			}
		}
		return exon_tree;
	}
}
