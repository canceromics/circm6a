package main;

import java.util.ArrayList;
import java.util.HashMap;

import htsjdk.samtools.util.IntervalTree;

public class Main {
	
	public static void main (String[] args){
		InParam in_args = new InParam(args);
		if (in_args.completeArgs()) {
			FileRW.printNow("Start at");

			ArrayList<HashMap<String,JuncInfo>> juncTable = null;// junctions keyed by their start position and end position split by the chr
			ArrayList<IntervalTree<ReadInfo>> itree = new ArrayList<>();// tree for all reads
			ArrayList<IntervalTree<ReadInfo>> itree_input = new ArrayList<>();
			ArrayList<Integer> chr_lengths = new ArrayList<>();
			chr_lengths.add(250000000);
			chr_lengths.add(245000000);
			chr_lengths.add(200000000);
			chr_lengths.add(200000000);
			chr_lengths.add(180000000);
			chr_lengths.add(175000000);
			chr_lengths.add(160000000);
			chr_lengths.add(150000000);
			chr_lengths.add(145000000);
			chr_lengths.add(140000000);
			chr_lengths.add(140000000);
			chr_lengths.add(135000000);
			chr_lengths.add(120000000);
			chr_lengths.add(110000000);
			chr_lengths.add(105000000);
			chr_lengths.add(95000000);
			chr_lengths.add(85000000);
			chr_lengths.add(80000000);
			chr_lengths.add(65000000);
			chr_lengths.add(60000000);
			chr_lengths.add(55000000);
			chr_lengths.add(50000000);
			chr_lengths.add(160000000);
			chr_lengths.add(60000000);
			chr_lengths.add(20000);
			for (int i = 0; i < 25; ++i) {
				IntervalTree<ReadInfo> temp = new IntervalTree<>();
				temp.setSentinel(null);
				itree.add(temp);
				temp = new IntervalTree<>();
				temp.setSentinel(null);
				itree_input.add(temp);
			}
			FileRW.halfDev = in_args.getRead_dev();
			FileRW.sup_read = in_args.getSup_reads();
			if (in_args.getCirc_bed() != null) {
				juncTable = FileRW.readJuncs(in_args.getCirc_bed());
			}
			// first time read sam
			if (in_args.isRetain_test()) {
				juncTable = FileRW.filtTxt(in_args, itree);
			}
			else {
				juncTable = FileRW.filtBam(in_args, itree, itree_input);
			}
			// report number of back junctions
			int sum = 0;
			for (int i = 0; i < juncTable.size(); ++i){
				sum += juncTable.get(i).size();
			}
			System.out.println("Total Possible BSJ:" + sum);
			//  filt junctions within gene and fix their site to exon bound if possible
			ArrayList<IntervalTree<ArrayList<Gene>>> genes = new ArrayList<>();
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