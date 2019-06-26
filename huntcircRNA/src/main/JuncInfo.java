package main;

import java.util.ArrayList;
import java.util.HashSet;

public class JuncInfo {
	private ArrayList<String> readids=new ArrayList<String>();
	private ArrayList<String> inputids=new ArrayList<String>();
	private HashSet<String> single_ip_ids=null;
	private HashSet<String> single_input_ids=null;
	private ArrayList<String> genes = new ArrayList<String>();
	private ArrayList<Integer> exons= new ArrayList<>();
	private int startPos=0;
	private int endPos=0;
	private int single_ip_reads=0;
	private int single_input_reads=0;
	private int totalReads=0;
	private int inputReads=0;
	private int[] intron=null;
	private char strand='.';
	private boolean intron_flag=false;
	private boolean fix_start_exon=false;
	private boolean fix_end_exon=false;
	
	public JuncInfo(){
		this.single_ip_ids = new HashSet<>();
		this.single_input_ids = new HashSet<>();
		this.startPos = 0;
		this.endPos = 0;
		this.totalReads = 0;
	}
	
	public JuncInfo(int SP, int EP) {
		this.startPos = SP;
		this.endPos = EP;
		this.totalReads = 0;
	}
	
	public ArrayList<String> getReadID(){
		return this.readids;
	}
	
	public void setReadID(ArrayList<String> readids){
		this.readids = readids;
	}
	
	public void addReadID(String readid){
		if (!this.readids.contains(readid)) {
			this.readids.add(readid);
		}
	}
	
	public void addReadID(String readid, boolean ip_flag){
		if (ip_flag) {
			if (readid == null || !this.readids.contains(readid)) {
				this.readids.add(readid);
			}
		}
		else {
			if (readid == null || !this.inputids.contains(readid)) {
				this.inputids.add(readid);
			}
		}
	}
	
	public void addReadID(ArrayList<String> readid){
		this.readids.addAll(readid);
	}
	
	public int getSP(){
		return this.startPos;
	}
	
	public void setSP(Integer SP){
		this.startPos = SP;
	}
	
	public int getEP(){
		return this.endPos;
	}
	
	public void setEP(Integer EP){
		this.endPos = EP;
	}
	
	public int getTR() {
		return this.totalReads;
	}
	
	public void setTR(int TR) {
		this.totalReads = TR;
	}
	
	public void incTR() {
		this.totalReads++;
	}
	
	public void addReads(int reads, boolean ip_flag) {
		if (ip_flag) {
			this.totalReads += reads;
		}
		else {
			this.inputReads += reads;
		}
	}
	
	public ArrayList<String> getInputids() {
		return inputids;
	}

	public void setInputids(ArrayList<String> inputids) {
		this.inputids = inputids;
	}
	
	public void addInputID(String readid){
		if (!this.inputids.contains(readid)) {
			this.inputids.add(readid);
		}
	}
	
	public void addInputID(ArrayList<String> readid){
		this.inputids.addAll(readid);
	}
	
	public int getInputReads() {
		return inputReads;
	}

	public void setInputReads(Integer inputreads) {
		this.inputReads = inputreads;
	}

	public ArrayList<String> getGenes(){
		return genes;
	}
	
	public void setGenes(ArrayList<String> genes) {
		this.genes = genes;
	}
	
	public void addGene(String gene) {
		this.genes.add(gene);
	}
	
	public ArrayList<Integer> getExons() {
		return exons;
	}

	public void setExons(ArrayList<Integer> exons) {
		this.exons = exons;
	}
	
	public HashSet<String> getSingle_ip_ids() {
		return single_ip_ids;
	}

	public HashSet<String> getSingle_input_ids() {
		return single_input_ids;
	}

	public int[] getIntron() {
		return intron;
	}

	public void setIntron(int[] intron) {
		this.intron = intron;
	}

	public char getStrand() {
		return strand;
	}

	public void setStrand(char strand) {
		this.strand = strand;
	}

	public int getSingle_ip_reads() {
		return single_ip_reads;
	}

	public void setSingle_ip_reads(int single_ip_reads) {
		this.single_ip_reads = single_ip_reads;
	}

	public int getSingle_input_reads() {
		return single_input_reads;
	}

	public void setSingle_input_reads(int single_input_reads) {
		this.single_input_reads = single_input_reads;
	}

	public boolean isIntron_flag() {
		return intron_flag;
	}

	public void setIntron_flag(boolean intron_flag) {
		this.intron_flag = intron_flag;
	}

	public boolean isFix_start_exon() {
		return fix_start_exon;
	}

	public void setFix_start_exon(boolean fix_start_exon) {
		this.fix_start_exon = fix_start_exon;
	}

	public boolean isFix_end_exon() {
		return fix_end_exon;
	}

	public void setFix_end_exon(boolean fix_end_exon) {
		this.fix_end_exon = fix_end_exon;
	}
	
}
