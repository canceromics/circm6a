package main;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import htsjdk.samtools.util.IntervalTree;

public class Main {
	
	public static void main (String[] args){
		InParam in_args = new InParam(args);
		if (in_args.completeArgs()) {
			FileRW.printNow("Start at");
			
			HashMap<String, HashMap<String,JuncInfo>> juncTable = new HashMap<>();// junctions keyed by their start position and end position split by the chr
			HashMap<String, IntervalTree<ReadInfo>> itree = new HashMap<>();// tree for all reads
			HashMap<String, IntervalTree<ReadInfo>> itree_input = new HashMap<>();
			HashMap<String, Integer> chr_lengths = new HashMap<>();
			chr_lengths.put("chr1", 250000000);
			chr_lengths.put("chr2", 245000000);
			chr_lengths.put("chr3", 200000000);
			chr_lengths.put("chr4", 200000000);
			chr_lengths.put("chr5", 180000000);
			chr_lengths.put("chr6", 175000000);
			chr_lengths.put("chr7", 160000000);
			chr_lengths.put("chr8", 150000000);
			chr_lengths.put("chr9", 145000000);
			chr_lengths.put("chr10", 140000000);
			chr_lengths.put("chr11", 140000000);
			chr_lengths.put("chr12", 135000000);
			chr_lengths.put("chr13", 120000000);
			chr_lengths.put("chr14", 110000000);
			chr_lengths.put("chr15", 105000000);
			chr_lengths.put("chr16", 95000000);
			chr_lengths.put("chr17", 85000000);
			chr_lengths.put("chr18", 80000000);
			chr_lengths.put("chr19", 65000000);
			chr_lengths.put("chr20", 60000000);
			chr_lengths.put("chr21", 55000000);
			chr_lengths.put("chr22", 50000000);
			chr_lengths.put("chrX", 160000000);
			chr_lengths.put("chrY", 60000000);
			chr_lengths.put("chrM", 20000);
			FileRW.halfDev = in_args.getRead_dev();
			FileRW.sup_read = in_args.getSup_reads();
			if (in_args.getCirc_bed() != null) {
				 FileRW.readJuncs(in_args.getCirc_bed(), juncTable);
			}
			// first time read sam
//			if (in_args.isRetain_test()) {
//				juncTable = FileRW.filtTxt(in_args, itree);
//			}
//			else {
				 FileRW.filtBam(juncTable, in_args, itree, itree_input);
//			}
			// report number of back junctions
			int sum = 0;
			for (Entry<String, HashMap<String, JuncInfo>> entry : juncTable.entrySet()){
				sum += entry.getValue().size();
			}
			System.out.println("Total Possible BSJ:" + sum);
			//  filt junctions within gene and fix their site to exon bound if possible
			HashMap<String, IntervalTree<ArrayList<Gene>>> genes = new HashMap<>();
			FileRW.loadGenes(in_args, juncTable, genes);
			// filt paired clipping signal such as GT-AG
			FileRW.filtGTAG(in_args.getGenome_file(), juncTable, chr_lengths);
			FileRW.countReadsExon(juncTable, itree, itree_input);
			// transform junctions to strings for output
			FileRW.printNow("Calculate p-value at");
			ArrayList<String> outList = null;
			outList = FileRW.juncsToBed(juncTable, in_args);
			System.out.println("Number of BSJ: " + (outList.size() - 1));
			FileRW.fileWrite(in_args.getOut_prefix() + "_circ.bed", outList);
			// writing details of these junctions
			if (in_args.isOut_detail()) {
				outList = FileRW.juncsToArray(in_args, juncTable, itree, itree_input);
				FileRW.fileWrite(in_args.getOut_prefix() + "_circ_detail.bed", outList);
			}
			// peakcalling using information that stored
			FileRW.readTrim(in_args, itree_input);
			FileRW.calPeak(in_args, itree, itree_input, juncTable, genes, chr_lengths);
			
			FileRW.printNow("Finish at");
		}
		System.out.println("END");
	}
	
}