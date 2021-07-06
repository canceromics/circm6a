package main;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;

import genome.Exon;
import genome.Gene;
import genome.IntRegion;
import genome.Transcript;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import mapping.Alignment;
import mapping.Junction;
import mapping.MappingMethod;
import mapping.MappingStat;
import mapping.MissingAlignment;
import output.Bed3;
import output.CircOut;
import output.PeakOut;

public class Method {
	
	private static Method m = new Method();
	private boolean ip_flag = false;
//	private boolean debug_flag = false;
//	private boolean search_flag = true;
	private RemoverRNA r = null;
	public final String[] AG = {"AG","AC","AG"};
	public final String[] GT = {"GT","AT","GC"};
	public final String[] AG_neg = {"AC","AT","GC"};
	public final String[] GT_neg = {"CT","GT","CT"};
	private HashMap<String, String> chr_seqs = new HashMap<>();
	private HashMap<String, Integer> chr_lengths = null;
	private HashMap<String, HashMap<Integer, Integer>> exon_boundary_map = new HashMap<>();
	
	private Method() {
	}

	public static void run() {
		
		m.r = new RemoverRNA(InParam.getParams().getRrna_bed()); // initial rRNA removing function
		
		HashMap<String, IntervalTree<Junction>> junc_table = new HashMap<>();
		HashMap<String, IntervalTree<MappingStat>> mapping_tree = new HashMap<>();
		HashMap<String, IntervalTree<MissingAlignment>> missing_tree = new HashMap<>();

		m.chr_seqs = m.loadGenome(InParam.getParams().getGenome_file());
		m.chr_lengths = m.getDefaultChrLength(null);
		HashMap<String, IntervalTree<Gene>> genes = m.loadGenes(InParam.getParams().getGtf_file());
		
		if (InParam.getParams().getCirc_bed() == null) {
			if (!m.handleBam(junc_table, mapping_tree, missing_tree, true)) {
				return;
			}
			m.handleMis(junc_table, mapping_tree, missing_tree);
			Bed3.setChrOrder(m.chr_lengths.keySet());
	
			m.fixJunctionBounder(genes, junc_table);
			
			fileWrite(InParam.getParams().getOut_prefix() + "_circRNAs.txt", m.juncsToTXT(junc_table, mapping_tree));
			m.handleTrim(junc_table, mapping_tree);
		}
		else {
			junc_table = m.loadCircRNA(InParam.getParams().getCirc_bed());
//			m.search_flag = false;
			if (!m.handleBam(junc_table, mapping_tree, missing_tree, false)) {
				return;
			}
			m.handleMis(junc_table, mapping_tree, missing_tree);
			Bed3.setChrOrder(m.chr_lengths.keySet());
		}
		m.chr_seqs.clear();
		m.exon_boundary_map.clear();
		m.calPeak(mapping_tree, junc_table, genes);
		
	}
	
	private HashMap<String, String> loadGenome(String genome_file) {
		HashMap<String, String> chr_seqs = new HashMap<>();
		if (genome_file == null) {
			System.out.println("No Genome File For GT-AG Signal");
			return chr_seqs;
		}
		try (BufferedReader reader = new BufferedReader(new FileReader(new File(genome_file)))){
			String tempString = null;
			StringBuffer bases = new StringBuffer();
			String chr = null;
			while ((tempString = reader.readLine()) != null){
				if (tempString.charAt(0) == '>'){
					if (bases.length() > 0) {
						chr_seqs.put(chr, bases.toString());
					}
					String[] temp_chr = tempString.split("\\s+");
					chr = temp_chr[0].substring(1);
					bases.setLength(0);
				}
				else {
					bases.append(tempString);
				}
			}
			if (bases.length() > 0) {
				chr_seqs.put(chr, bases.toString());
			}
			reader.close();
		} catch(IOException e){
			e.printStackTrace();
		}
		return chr_seqs;
	}

	public HashMap<String, Integer> getDefaultChrLength(SAMFileHeader header){
		HashMap<String, Integer> chr_lengths = new HashMap<>();
		if (header == null) {
			if (chr_seqs == null) {
//				chr_lengths.put("chr1", 250000000);
//				chr_lengths.put("chr2", 245000000);
//				chr_lengths.put("chr3", 200000000);
//				chr_lengths.put("chr4", 200000000);
//				chr_lengths.put("chr5", 180000000);
//				chr_lengths.put("chr6", 175000000);
//				chr_lengths.put("chr7", 160000000);
//				chr_lengths.put("chr8", 150000000);
//				chr_lengths.put("chr9", 145000000);
//				chr_lengths.put("chr10", 140000000);
//				chr_lengths.put("chr11", 140000000);
//				chr_lengths.put("chr12", 135000000);
//				chr_lengths.put("chr13", 120000000);
//				chr_lengths.put("chr14", 110000000);
//				chr_lengths.put("chr15", 105000000);
//				chr_lengths.put("chr16", 95000000);
//				chr_lengths.put("chr17", 85000000);
//				chr_lengths.put("chr18", 80000000);
//				chr_lengths.put("chr19", 65000000);
//				chr_lengths.put("chr20", 60000000);
//				chr_lengths.put("chr21", 55000000);
//				chr_lengths.put("chr22", 50000000);
//				chr_lengths.put("chrX", 160000000);
//				chr_lengths.put("chrY", 60000000);
//				chr_lengths.put("chrM", 20000);
			}
			else {
				for (Entry<String, String> entry : chr_seqs.entrySet()) {
					chr_lengths.put(entry.getKey(), entry.getValue().length());
				}
			}
		}
		else {
			for (SAMSequenceRecord srecord : header.getSequenceDictionary().getSequences()) {
				String chr = srecord.getSequenceName();
				chr_lengths.put(chr, getChrLength(chr) == 0 ? srecord.getSequenceLength() : 
					Math.min(srecord.getSequenceLength(), getChrLength(chr)));
			}
		}
		return chr_lengths;
	}

	private HashMap<String, IntervalTree<Gene>> loadGenes(String file_name){
		HashMap<String, IntervalTree<Gene>> out = new HashMap<>();
		if (file_name == null) {
			return out;
		}
		try (BufferedReader reader = new BufferedReader(new FileReader(new File(file_name)))){
			String line = null;
			Gene gene = null;
			Transcript script = null;
			String gene_id = "gene_id";
			String gene_symbol = "gene_name";
			String script_id = "script_id";
			char sep = '"';
			ArrayList<Gene> unsort_genes = new ArrayList<>();
			while ((line = reader.readLine()) != null) {
				String[] cols = line.split("\t");
				if (cols.length < 7) {
					continue;
				}
				String key = cols[2].toLowerCase();
				String chr = cols[0];
				int start = Integer.parseInt(cols[3]);
				int end = Integer.parseInt(cols[4]);
				char strand = cols[6].charAt(0);
				if ("gene".equals(key)) {
					if (!out.containsKey(chr)) {
						out.put(chr, new IntervalTree<>());
					}
					gene = out.get(chr).put(start, end, null);
					if (gene == null) {
						gene = new Gene(chr, getIDSym(line, gene_id, sep), getIDSym(line, gene_symbol, sep), strand, start, end, new ArrayList<>());
					}
					else {
						gene.appendGeneID(getIDSym(line, gene_id, sep));
						gene.appendGeneSymbol(getIDSym(line, gene_symbol, sep));
						gene.checkStrand(strand);
					}
					out.get(chr).put(start, end, gene);
				}
				else if ("transcript".equals(key)) {
					if (gene == null || gene.getGene_id().indexOf(getIDSym(line, gene_id, sep)) == -1){
						gene = new Gene(chr, getIDSym(line, gene_id, sep), getIDSym(line, gene_symbol, sep), strand, start, end, new ArrayList<>());
						unsort_genes.add(gene);
					}
					script = new Transcript(gene, getIDSym(line, script_id, sep), start, end, new ArrayList<>());
					gene.getScripts().add(script);
				}
				else if ("exon".equals(key)) {
					if (gene == null || gene.getGene_id().indexOf(getIDSym(line, gene_id, sep)) == -1){
						gene = new Gene(chr, getIDSym(line, gene_id, sep), getIDSym(line, gene_symbol, sep), strand, start, end, new ArrayList<>());
						unsort_genes.add(gene);
						script = new Transcript(gene, getIDSym(line, script_id, sep), start, end, new ArrayList<>());
						gene.getScripts().add(script);
					}
					else if (script == null || script.getScript_id().indexOf(getIDSym(line, script_id, sep)) == -1) {
						script = new Transcript(gene, getIDSym(line, script_id, sep), start, end, new ArrayList<>());
						gene.getScripts().add(script);
					}
					Exon exon = new Exon(script, start, end);
					if (!exon_boundary_map.containsKey(chr)) {
						exon_boundary_map.put(chr, new HashMap<>());
					}
					exon_boundary_map.get(chr).put(start, exon_boundary_map.get(chr).containsKey(start) ? exon_boundary_map.get(chr).get(start) & 1 : 1);
					exon_boundary_map.get(chr).put(end, exon_boundary_map.get(chr).containsKey(end) ? exon_boundary_map.get(chr).get(end) & 2 : 2);
					script.getExons().add(exon);
				}
			}
			for (int i = 0; i < unsort_genes.size(); i++) {
				gene = unsort_genes.get(i);
				gene.checkPos();
				if (!out.containsKey(gene.getChr_symbol())) {
					out.put(gene.getChr_symbol(), new IntervalTree<>());
				}
				Gene old = out.get(gene.getChr_symbol()).put(gene.getStart(), gene.getEnd(), gene);
				if (old != null) {
					gene.mergeGene(old);
				}
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return out;
	}

	private String getIDSym(String line, String key, char spl_sym) {
		int index = line.indexOf(key);
		if (index != -1) {
			index += key.length() + 2;
			return line.substring(index, line.indexOf(spl_sym, index));
		}
		return "";
	}

	private boolean handleBam(HashMap<String, IntervalTree<Junction>> junc_table, HashMap<String, IntervalTree<MappingStat>> mapping_tree,
			HashMap<String, IntervalTree<MissingAlignment>> missing_tree, boolean search) {
		String ip_file = InParam.getParams().getIp_file();
		String input_file = InParam.getParams().getInput_file();
		SAMFileHeader header = null;
		if (input_file != null) {
			ip_flag = false;
			Method.printNow("Scaning " + input_file + " at");
			try (SamReader reader = SamReaderFactory.makeDefault().open(new File(input_file))){
				header = reader.getFileHeader();
				if (!handleSAMHeader(header)) {
					System.err.println("Warning:input bam have different chrs or chr lengths with genome file");
					System.err.println("genome file:");
					chr_lengths.forEach((chr, len) -> {
						System.err.println(chr + "\t" + len);
					});
					chr_lengths = getDefaultChrLength(reader.getFileHeader());
					System.err.println("BAM file:");
					chr_lengths.forEach((chr, len) -> {
						System.err.println(chr + "\t" + len);
					});
				}
				if (search) {
					handlRecords(reader, junc_table, mapping_tree, missing_tree);
				}
				else {
					handleRecordsWithoutSearch(reader, junc_table, mapping_tree, missing_tree);
				}
				System.out.println(input_file + " total alignments: " + MappingStat.getInputReads());
				System.out.println(input_file + " total mapping alignments: " + MappingStat.getMappingInputReads());
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		if (ip_file != null) {
			ip_flag = true;
			Method.printNow("Scaning " + ip_file + " at");
			try (SamReader reader = SamReaderFactory.makeDefault().open(new File(ip_file))){
				if (getDefaultChrLength(header).equals(getDefaultChrLength(reader.getFileHeader()))) {
					if (search) {
						handlRecords(reader, junc_table, mapping_tree, missing_tree);
					}
					else {
						handleRecordsWithoutSearch(reader, junc_table, mapping_tree, missing_tree);
					}
					System.out.println(ip_file + " total alignments: " + MappingStat.getIpReads());
					System.out.println(ip_file + " total mapping alignments: " + MappingStat.getMappingIpReads());
				}
				else {
					System.err.println("Error:IP bam have different chrs or chr lengths with input bam");
					return false;
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		return true;
	}

	private boolean handleSAMHeader(SAMFileHeader header) {
		if (chr_lengths == null) {
			chr_lengths = getDefaultChrLength(header);
		}
		else {
			return chr_lengths.equals(getDefaultChrLength(header));
		}
		return true;
	}

	private int handlRecords(SamReader reader, HashMap<String, IntervalTree<Junction>> junc_table, 
			HashMap<String, IntervalTree<MappingStat>> mapping_tree, HashMap<String, IntervalTree<MissingAlignment>> missing_tree) {
		SAMRecordIterator samit = reader.iterator();
		String id = null;
		SAMRecord record = null;
		Alignment read = null;
		String unmap_id = null;
		int reads = 0;
		while (samit.hasNext()) {
			record = samit.next();
			if (record.getReadUnmappedFlag() || record.getReadName().equals(unmap_id)) {
				unmap_id = record.getReadName();
				if (read != null && !unmap_id.equals(read.getID())) {
					read.handleReads(junc_table, mapping_tree, missing_tree, true);
					read = null;
				}
				continue;
			}
			if (++reads % 1000000 == 0) {
				printNow("Runing " + reads + " alignments at");
			}
			if (record.getReadName().equals(id)) {
				putRecordToRead(record, read, mapping_tree);
			}
			else {
				if (read != null) {
					read.handleReads(junc_table, mapping_tree, missing_tree, true);
				}
				id = record.getReadName();
				read = new Alignment(0, 0, id);
				putRecordToRead(record, read, mapping_tree);
			}
		}
		if (read != null) {
			read.handleReads(junc_table, mapping_tree, missing_tree, true);
		}
		return reads;
	}

	private void putRecordToRead(SAMRecord record, Alignment read, HashMap<String, IntervalTree<MappingStat>> mapping_tree) {
		if (record.getReadPairedFlag() && record.getSecondOfPairFlag()) {
			if (read.getPairRead() == null) {
				read.setPairRead(new Alignment(0, 0, read.getID()));
				read.getPairRead().setPairRead(read);
			}
			if (!read.getPairRead().addSegs(record, InParam.getParams().getMapQuality(), InParam.getParams().isUniq_mode())) {
				putQualityMapping(record, mapping_tree);
			}
		}
		else {
			if (!read.addSegs(record, InParam.getParams().getMapQuality(), InParam.getParams().isUniq_mode())) {
				putQualityMapping(record, mapping_tree);
			}
		}
	}

	private void putQualityMapping(SAMRecord record, HashMap<String, IntervalTree<MappingStat>> mapping_tree) {
		if (!mapping_tree.containsKey(record.getReferenceName())) {
			mapping_tree.put(record.getReferenceName(), new IntervalTree<>());
		}
		MappingMethod.putMappingStat(mapping_tree.get(record.getReferenceName()), record.getAlignmentStart(), record.getAlignmentEnd(), ip_flag, record.getMappingQuality());
	}

	private void handleMis(HashMap<String, IntervalTree<Junction>> junc_table, HashMap<String, IntervalTree<MappingStat>> mapping_tree, HashMap<String, IntervalTree<MissingAlignment>> mis_table) {
			junc_table.forEach((chr, junc_tree) -> {
				if (mis_table.containsKey(chr)) {
					for (Node<Junction> junc_node : junc_tree) {
						Iterator<Node<MissingAlignment>> nodes = mis_table.get(chr).overlappers(junc_node.getStart(), junc_node.getStart());
						while (nodes.hasNext()) {
							Node<MissingAlignment> node = nodes.next();
							Junction junc = junc_node.getValue();
							MissingAlignment mis = node.getValue();
							while (mis != null) {
								if (mis.isFrontMissing()) {
									junc.incLoose(mis.isIP());
									if (isCircClip(getChrSeq(chr, Math.max(junc_node.getStart() - 1, junc_node.getEnd() - 20000), junc_node.getEnd()),
										getChrSeq(chr, Math.max(0, junc_node.getStart() - 50000), junc_node.getStart() - 1), mis.getSeq())) {
										junc.incJunc(mis.isIP());
										if (!mis.getSeg().isOnceCirc()) {
											MappingMethod.putCircMappingStat(mapping_tree.get(chr), mis.getSeg().getStart(), mis.getSeg().getEnd(), mis.isIP());
										}
									}
								}
								mis = mis.getAnotherMissingAlignment();
							}
						}
						nodes = mis_table.get(chr).overlappers(junc_node.getEnd(), junc_node.getEnd());
						while (nodes.hasNext()) {
							Node<MissingAlignment> node = nodes.next();
							Junction junc = junc_node.getValue();
							MissingAlignment mis = node.getValue();
							while (mis != null) {
								if (!mis.isFrontMissing()) {
									junc.incLoose(mis.isIP());
									if (isCircClip(getChrSeq(chr, junc_node.getStart() - 1, Math.min(junc_node.getStart() + 20000, junc_node.getEnd())),
										getChrSeq(chr, junc_node.getEnd(), Math.min(junc_node.getEnd() + 50000, getChrLength(chr))), mis.getSeq())) {
										junc.incJunc(mis.isIP());
										if (!mis.getSeg().isOnceCirc()) {
											MappingMethod.putCircMappingStat(mapping_tree.get(chr), mis.getSeg().getStart(), mis.getSeg().getEnd(), mis.isIP());
										}
									}
								}
								mis = mis.getAnotherMissingAlignment();
							}
						}
					}
				}
			});
		}

	public boolean isCircClip(String circ_seq, String linear_seq, String clip_seg) {
		if (linear_seq == null) {
			return true;
		}
		if (circ_seq == null) {
			return false;
		}
		if (circ_seq.indexOf(clip_seg) != -1) {
			if (linear_seq.indexOf(clip_seg) != -1) {
				return false;
			}
			else {
				return true;
			}
		}
		else if (linear_seq.indexOf(clip_seg) != -1) {
			return false;
		}
		int[] base_units = {9, 7, 5, 4, 3};
		Boolean circ_flag = null;
		for (int step : base_units) {
			String[] seeds = new String[(clip_seg.length() - 1) / step];
			for (int i = 0; i < seeds.length; i++) {
				if (clip_seg.length() >= (i + 2) * step) {
					seeds[i] = clip_seg.substring(i * step, (i + 2) * step);
				}
				else {
					seeds[i] = clip_seg.substring(Math.max(0, clip_seg.length() - 2 * step), clip_seg.length());
				}
			}
			circ_flag = isCircClip(circ_seq, linear_seq, seeds);
			if (circ_flag != null) {
				return circ_flag;
			}
		}
		return false;
	}

	private Boolean isCircClip(String circ_seq, String linear_seq, String[] seeds) {
		int match_count = 0;
		int match_dis = 0;
		int last_index = -1;
		for (String seed : seeds) {
			int index = circ_seq.indexOf(seed, last_index + 1);
			if (index != -1) {
				match_count++;
				match_dis += last_index >= 0 ? Math.abs(index - last_index) * 20 : 0;
				last_index = index;
			}
		}
		boolean enough_seeds = match_count + 3 >= Math.max(4, seeds.length);
		boolean full_match = match_count == seeds.length;
		last_index = -1;
		for (String seed : seeds) {
			int index = linear_seq.indexOf(seed, last_index + 1);
			if (index != -1) {
				match_count--;
				match_dis -= last_index >= 0 ? Math.abs(index - last_index) : 0;
				last_index = index;
			}
		}
		if (match_count > 0 && enough_seeds) {
			return true;
		}
		if (match_count < 0) {
			return false;
		}
		if (match_count == 0 && enough_seeds && full_match) {
			if (match_dis > 0) {
				return false;
			}
			else {
				return true;
			}
		}
		return null;
	}

	private void fixJunctionBounder(HashMap<String, IntervalTree<Gene>> genes, HashMap<String, IntervalTree<Junction>> junctions) {
		for (Entry<String, IntervalTree<Junction>> entry : junctions.entrySet()) {
			if (genes.containsKey(entry.getKey())) {
				Iterator<Node<Junction>> nodes = entry.getValue().iterator();
				while (nodes.hasNext()) {
					Node<Junction> node = nodes.next();
					if (node.getValue().checkBounderAndGene(false)) {
						continue;
					}
					Iterator<Node<Gene>> gene_nodes = genes.get(entry.getKey()).overlappers(node.getStart() - getDev(), node.getEnd() + getDev());
					while (gene_nodes.hasNext()) {
						Node<Gene> gene_node = gene_nodes.next();
						gene_node.getValue().checkPos();
						if (gene_node.getEnd() + getDev() >= node.getEnd() && gene_node.getStart() - getDev() <= node.getStart()) {
							node.getValue().resetGene(gene_node.getValue());
						}
					}
				}
			}
		}
	}

	private List<String> juncsToTXT(HashMap<String, IntervalTree<Junction>> junc_table, HashMap<String, IntervalTree<MappingStat>> mapping_stat){
		List<String> out = new ArrayList<String>();
		out.add(CircOut.getHeader(InParam.getParams().getIp_file() != null));
		for (String chr : Bed3.getChrsInOrder()) {
			IntervalTree<MappingStat> mapping_tree = mapping_stat.getOrDefault(chr, new IntervalTree<>());
			Iterator<Node<Junction>> nodes = junc_table.getOrDefault(chr, new IntervalTree<>()).overlappers(0, Integer.MAX_VALUE);
			while(nodes.hasNext()) {
				Node<Junction> node = nodes.next();
				if (node.getValue().reachOutputStandard(0)) {
					out.add(InParam.getParams().isOut_detail() ? 
						getOutPut(node, mapping_tree, chr).toString() + "\t" + node.getValue().getIDString()
							: getOutPut(node, mapping_tree, chr).toString());
				}
			}
		}
		return out;
	}

	private CircOut getOutPut(Node<Junction> node, IntervalTree<MappingStat> mapping_tree, String chr) {
		int start = 0;
		int end = 0;
		int start_ip = 0;
		int end_ip = 0;
		int total = 0;
		int total_ip = 0;
		Iterator<Node<MappingStat>> map_nodes = mapping_tree.overlappers(node.getStart() - getDev(), node.getStart() + getDev());
		while (map_nodes.hasNext()) {
			Node<MappingStat> map_node = map_nodes.next();
			if (Math.abs(map_node.getStart() - node.getStart()) <= getDev()) {
				start += map_node.getValue().getCircReads();
				start_ip += map_node.getValue().getCircReadsIP();
			}
			total += map_node.getValue().getReads();
			total_ip += map_node.getValue().getReadsIP();
		}
		map_nodes = mapping_tree.overlappers(node.getEnd() - getDev(), node.getEnd() + getDev());
		while (map_nodes.hasNext()) {
			Node<MappingStat> map_node = map_nodes.next();
			if (map_node.getStart() > node.getStart() + getDev()) {
				total += map_node.getValue().getReads();
				total_ip += map_node.getValue().getReadsIP();
				if (Math.abs(map_node.getEnd() - node.getEnd()) <= getDev()) {
					end += map_node.getValue().getCircReads();
					end_ip += map_node.getValue().getCircReadsIP();
				}
			}
		}
		return new CircOut(chr, node.getStart() - 1, node.getEnd(), node.getValue().getGene() == null ? "None" : node.getValue().getGene().getGene_symbol(), total + total_ip,
				node.getValue().getGene() == null ? '.' : node.getValue().getGene().getStrand(), start + end - node.getValue().getIDs().size(), 
						start_ip + end_ip - node.getValue().getIPIDs().size(), total - start - end + total_ip - start_ip - end_ip, (double) (start + end) / (double) total);
	}

	private void handleTrim(HashMap<String, IntervalTree<Junction>> junc_table, HashMap<String, IntervalTree<MappingStat>> mapping_stat) {
		if (InParam.getParams().getTrim_file() == null) {
			return;
		}
		ip_flag = false;
		MappingMethod.resetMappingStat(mapping_stat, false);
		Method.printNow("Scaning " + InParam.getParams().getTrim_file() + " at");
		try (SamReader reader = SamReaderFactory.makeDefault().open(new File(InParam.getParams().getTrim_file()))){
			handleTrimRecords(reader, junc_table, mapping_stat);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void handleTrimRecords(SamReader reader, HashMap<String, IntervalTree<Junction>> junc_table,
			HashMap<String, IntervalTree<MappingStat>> mapping_tree) {
		SAMRecordIterator samit = reader.iterator();
		String id = null;
		String unmap_id = null;
		SAMRecord record = null;
		SAMRecord[] records = new SAMRecord[32];
		int len = 0;
		int reads = 0;
		while (samit.hasNext()) {
			record = samit.next();
			if (record.getReadUnmappedFlag() || record.getReadName().equals(unmap_id)) {
				unmap_id = record.getReadName();
				if (len != 0 && !records[0].getReadName().equals(unmap_id)) {
					handleTrimRecords(junc_table, mapping_tree, records, len);
					len = 0;
				}
				continue;
			}
			if (++reads % 1000000 == 0) {
				printNow("Runing " + reads + " alignments at");
			}
			if (record.getReadName().equals(id)) {
				records[len] = record;
				++len;
			}
			else {
				handleTrimRecords(junc_table, mapping_tree, records, len);
				id = record.getReadName();
				records[0] = record;
				len = 1;
			}
		}
		if (len != 0) {
			handleTrimRecords(junc_table, mapping_tree, records, len);
		}
	}

	private void handleTrimRecords(HashMap<String, IntervalTree<Junction>> junc_table,
			HashMap<String, IntervalTree<MappingStat>> mapping_tree, SAMRecord[] records, int len) {
		String seq = "";
		String mate_seq = "";
		for (int i = 0; i < len; i++) {
			if (records[i].getReadPairedFlag() && records[i].getSecondOfPairFlag()) {
				mate_seq = mate_seq.length() < records[i].getReadString().length() ? records[i].getReadString() : mate_seq;
			}
			else {
				seq = seq.length() < records[i].getReadString().length() ? records[i].getReadString() : seq;
			}
		}
		for (int i = 0; i < len; i++) {
			int[] cigar_value = Alignment.cigarToRegion(records[i].getCigarString());
			boolean mate_flag = records[i].getReadPairedFlag() && records[i].getSecondOfPairFlag();
			String chr = records[i].getReferenceName();
			if (cigar_value[0] > getDev()) {
				Iterator<Node<Junction>> junc_nodes = junc_table.getOrDefault(chr, new IntervalTree<>()).overlappers(
						records[i].getAlignmentStart() - getDev(), records[i].getAlignmentStart() + getDev());
				while (junc_nodes.hasNext()) {
					Node<Junction> junc_node = junc_nodes.next();
					if (Math.abs(records[i].getAlignmentEnd() - junc_node.getEnd()) < getDev() && Math.abs(records[i].getAlignmentStart() - junc_node.getStart()) < getDev() &&
							isCircClip(getChrSeq(chr, junc_node.getStart() - 1, Math.min(junc_node.getStart() + 20000, junc_node.getEnd())),
							getChrSeq(chr, Math.max(0, junc_node.getStart() - 50000), junc_node.getStart() - 1),
							(mate_flag ? mate_seq : seq).substring(0, cigar_value[0]))) {
						MappingMethod.putCircMappingStat(mapping_tree.get(chr), records[i].getAlignmentStart(), records[i].getAlignmentEnd(), isIP());
						break;
					}
				}
			}
			else if (cigar_value[1] < (mate_flag ? mate_seq : seq).length() - getDev()){
				Iterator<Node<Junction>> junc_nodes = junc_table.getOrDefault(chr, new IntervalTree<>()).overlappers(
						records[i].getAlignmentStart() - getDev(), records[i].getAlignmentStart() + getDev());
				while (junc_nodes.hasNext()) {
					Node<Junction> junc_node = junc_nodes.next();
					if (Math.abs(records[i].getAlignmentStart() - junc_node.getStart()) < getDev() && Math.abs(records[i].getAlignmentStart() - junc_node.getStart()) < getDev() &&
							isCircClip(getChrSeq(chr, junc_node.getStart() - 1, Math.min(junc_node.getStart() + 20000, junc_node.getEnd())),
							getChrSeq(chr, junc_node.getEnd(), Math.min(junc_node.getEnd() + 50000, getChrLength(chr))),
							(mate_flag ? mate_seq : seq).substring(cigar_value[1], (mate_flag ? mate_seq : seq).length()))) {
						MappingMethod.putCircMappingStat(mapping_tree.get(chr), records[i].getAlignmentStart(), records[i].getAlignmentEnd(), isIP());
						break;
					}
	
				}
			}
			MappingMethod.putMappingStat(mapping_tree.get(chr), records[i].getAlignmentStart(), records[i].getAlignmentEnd(), ip_flag, records[i].getMappingQuality());
		}
	}

	private HashMap<String, IntervalTree<Junction>> loadCircRNA(String circ_file){
		HashMap<String, IntervalTree<Junction>> junc_table = new HashMap<>();
		try (BufferedReader reader = new BufferedReader(new FileReader(new File(circ_file)))){
			String tempString = null;
			while ((tempString = reader.readLine()) != null){
				if (tempString.charAt(0) == '#'){
					continue;
				}
				String[] cols = tempString.split("\t");
				if (cols.length < 3) {
					continue;
				}
				if (!junc_table.containsKey(cols[0])) {
					junc_table.put(cols[0], new IntervalTree<>());
				}
				int start = Integer.parseInt(cols[1]);
				int end = Integer.parseInt(cols[2]);
				junc_table.get(cols[0]).put(start + 1, end, new Junction(start + 1, end));
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return junc_table;
	}
	
	private int handleRecordsWithoutSearch(SamReader reader, HashMap<String, IntervalTree<Junction>> junc_table,
			HashMap<String, IntervalTree<MappingStat>> mapping_tree, HashMap<String, IntervalTree<MissingAlignment>> missing_tree) {
		SAMRecordIterator samit = reader.iterator();
		String id = null;
		SAMRecord record = null;
		Alignment read = null;
		String unmap_id = null;
		int reads = 0;
		boolean search = false;
		while (samit.hasNext()) {
			record = samit.next();
			if (record.getReadUnmappedFlag() || record.getReadName().equals(unmap_id)) {
				unmap_id = record.getReadName();
				if (read != null && !unmap_id.equals(read.getID())) {
					read.handleReads(junc_table, mapping_tree, missing_tree, search);
					read = null;
				}
				continue;
			}
			if (++reads % 1000000 == 0) {
				printNow("Runing " + reads + " alignments at");
			}
			if (record.getReadName().equals(id)) {
				putRecordToRead(record, read, mapping_tree);
				search |= isCircRegion(junc_table, record);
			}
			else {
				if (read != null) {
					read.handleReads(junc_table, mapping_tree, missing_tree, search);
				}
				id = record.getReadName();
				read = new Alignment(0, 0, id);
				putRecordToRead(record, read, mapping_tree);
				search |= isCircRegion(junc_table, record);
			}
		}
		if (read != null) {
			read.handleReads(junc_table, mapping_tree, missing_tree, search);
		}
		return reads;
	}
	
	private boolean isCircRegion(HashMap<String, IntervalTree<Junction>> junc_table, SAMRecord record) {
		IntervalTree<Junction> junc_tree = junc_table.get(record.getReadName());
		if (junc_tree == null) {
			return false;
		}
		Iterator<Node<Junction>> junc_nodes = junc_tree.overlappers(record.getAlignmentStart(), record.getAlignmentEnd());
		while (junc_nodes.hasNext()) {
			Node<Junction> junc_node = junc_nodes.next();
			if (Math.abs(junc_node.getStart() - record.getAlignmentStart()) <= getDev() 
					|| Math.abs(junc_node.getEnd() - record.getAlignmentEnd()) <= getDev()) {
				return true;
			}
		}
		return false;
	}
	
	private void calPeak(HashMap<String, IntervalTree<MappingStat>> mapping_table, HashMap<String, IntervalTree<Junction>> junc_table, HashMap<String, IntervalTree<Gene>> gene_table){
		if (MappingStat.getInputReads() == 0 || MappingStat.getIpReads() == 0) {
			return;
		}
		
		FisherTest fisher_test = MappingStat.getInputReads() + MappingStat.getIpReads() + 32384 > Integer.MAX_VALUE ? 
				new FisherTest(MappingStat.getIpReads(), MappingStat.getInputReads(), 32384) :
				new FisherTest((int) (MappingStat.getInputReads() + MappingStat.getIpReads() + 16192));
		
		int window_size = InParam.getParams().getWindow_size();
		
		String peak_file = InParam.getParams().getOut_prefix() + "_linear_peak.bed";
		String circ_peak_file = InParam.getParams().getOut_prefix() + "_circ_peak.bed";
		if (chr_lengths == null) {
			chr_lengths = getDefaultChrLength(null);
		}

//		HashMap<String, ArrayList<PeakOut>> circ_map = new HashMap<>();
//		circ_map.put("bsj", new ArrayList<>());
//		double propor = InParam.getParams().getPropor();
//		for (double i = 0.0; i <= propor; i += 0.01) {
//			circ_map.put(i + "_2", new ArrayList<>());
//			circ_map.put(i + "_0", new ArrayList<>());
//		}
		ArrayList<PeakOut> circ_peak = new ArrayList<>();
		ArrayList<PeakOut> linear_peak = new ArrayList<>();
		for (String chr : Bed3.getChrsInOrder()) {
			if (junc_table.containsKey(chr)) {
				Iterator<Node<Junction>> nodes = junc_table.getOrDefault(chr, new IntervalTree<>()).overlappers(0, Integer.MAX_VALUE);
				while (nodes.hasNext()) {
					Node<Junction> node = nodes.next();
//					this.peak_read = 'c';
					PeakOut peak = p_valueCircPeak(fisher_test, node, mapping_table.get(chr), chr);
//					if (peak != null && (InParam.getParams().isOut_detail() || peak.getConfidence() > 0)) {
//						circ_map.get("bsj").add(peak);
//					}
//					this.peak_read = 't';
//					peak = p_valueCircPeak(fisher_test, node, mapping_table.get(chr), chr);
					if (peak != null && (InParam.getParams().isOut_detail() || peak.getConfidence() > 0)) {
//						for (double i = 0.0; i <= propor; i += 0.01) {
//							InParam.getParams().setPropor(i);
//							circ_map.get(i + "_2").add(peak.confidencePeak(node.getValue().getConfidence(2, 2)));
//							circ_map.get(i + "_0").add(peak.confidencePeak(node.getValue().getConfidence(0, 0)));
//						}
						circ_peak.add(peak);
					}
				}
			}
		}
		
//		m.search_flag = true;
		fisher_test = MappingStat.getInputReads() + MappingStat.getIpReads() + 32384 > Integer.MAX_VALUE ? 
				new FisherTest(MappingStat.getMappingIpReads(), MappingStat.getMappingInputReads(), 32384) : fisher_test;
		HashMap<String, IntervalTree<Double>> p_value_tree = calPvalueTree(fisher_test, window_size, mapping_table);
		for (String chr : Bed3.getChrsInOrder()) {
			List<Double> p_value_list = new ArrayList<>();
			int peak_start = 1;
			int peak_end = 1;
			Iterator<Node<Double>> p_nodes = p_value_tree.getOrDefault(chr, new IntervalTree<>()).overlappers(Integer.MIN_VALUE, Integer.MAX_VALUE);
			while (p_nodes.hasNext()) {
				Node<Double> p_node = p_nodes.next();
				if (p_node.getStart() <= peak_end + 1) {
					p_value_list.add(p_node.getValue());
				}
				else {
					if (peak_end - peak_start + 1 >= InParam.getParams().getPeak_length()) {
						List<IntRegion> list = new ArrayList<>();
						list.add(new Exon(null, peak_start - 1, peak_end));
						PeakOut peak = new PeakOut(chr, peak_start - 1, peak_end, findGeneSymbol(peak_start, peak_end, gene_table.get(chr)), FisherTest.combinePvalue(p_value_list), '.',
								peak_start - 1, peak_end, 0, 1, new int[] {peak_end - peak_start + 1}, new int[] {peak_start - 1}, null, 0.0, getPropor(list, mapping_table.get(chr)));
						linear_peak.add(peak);
					}
					peak_start = p_node.getStart();
					p_value_list = new ArrayList<>();
					p_value_list.add(p_node.getValue());
				}
				peak_end = p_node.getEnd();
			}
			if (peak_end - peak_start + 1 >= InParam.getParams().getPeak_length()) {
				List<IntRegion> list = new ArrayList<>();
				list.add(new Exon(null, peak_start - 1, peak_end));
				PeakOut peak = new PeakOut(chr, peak_start, peak_end, findGeneSymbol(peak_start, peak_end, gene_table.get(chr)), FisherTest.combinePvalue(p_value_list), '.',
						peak_start, peak_end, 0, 1, new int[] {peak_end - peak_start + 1}, new int[] {peak_start}, null, 0.0, getPropor(list, mapping_table.get(chr)));
				linear_peak.add(peak);
			}
		}
		
		if (InParam.getParams().getAdjsut_p() != null) {
//			circ_map.forEach((file, circ_peak) -> {
//				adjustP_Value(circ_peak, InParam.getParams().getAdjsut_p());
//				fileWriteWithHeader(InParam.getParams().getOut_prefix() + "_" + file + "_circ_peak.bed", filtPeak(circ_peak), PeakOut.getHeader());
//			});
			adjustP_Value(circ_peak, InParam.getParams().getAdjsut_p());
			fileWriteWithHeader(circ_peak_file, filtPeak(circ_peak), PeakOut.getLinearHeader());
			adjustP_Value(linear_peak, InParam.getParams().getAdjsut_p());
			fileWriteWithHeader(peak_file, filtPeak(linear_peak), PeakOut.getLinearHeader());
		}
	}
	
	private PeakOut p_valueCircPeak(FisherTest fisher_test, Node<Junction> junc_node, IntervalTree<MappingStat> mapping_tree, String chr) {
		List<Double> p_value_list = calCircPvalue(fisher_test, chr, junc_node.getStart(), junc_node.getEnd(), InParam.getParams().getWindow_size(), mapping_tree);
		int window_size = InParam.getParams().getWindow_size();
		int start_len = 0;
		int end_len = 0;
		boolean end_flag = false;
		for (int i = 0; i < p_value_list.size(); i++) {
			if (p_value_list.get(i) <= InParam.getParams().getP_value()) {
				if (end_flag) {
					end_len += window_size;
				}
				else {
					start_len += window_size;
				}
			}
			else {
				end_flag = true;
				p_value_list.remove(i);
				--i;
			}
		}
		if (start_len >= junc_node.getLength()) {
			String name = junc_node.getValue().getGene() != null ? junc_node.getValue().getGene().getGene_symbol() : "None";
			char strand = junc_node.getValue().getGene() != null ? junc_node.getValue().getGene().getStrand() : '.';
			List<IntRegion> list = new ArrayList<>(2);
			list.add(new Exon(null, junc_node.getStart() - 1, junc_node.getStart() - 1 + window_size));
			list.add(new Exon(null, junc_node.getEnd() - window_size, junc_node.getEnd()));
			return new PeakOut(chr, junc_node.getStart() - 1, junc_node.getEnd(), name, FisherTest.combinePvalue(p_value_list), strand, junc_node.getStart() - 1, junc_node.getEnd(), 
					junc_node.getValue().getIPIDs().size(), 1, new int[]{junc_node.getLength()}, new int[] {junc_node.getStart() - 1}, 
					junc_node.getValue().getConfidence(3 * getIPSup(mapping_tree, junc_node.getStart(), junc_node.getEnd()),
					3 * getIPSup(mapping_tree, junc_node.getStart(), junc_node.getEnd())), 0.0, getPropor(list, mapping_tree));
		}
		if (start_len > 0 && end_len > 0 && start_len + end_len >= InParam.getParams().getPeak_length()) {
			String name = junc_node.getValue().getGene() != null ? junc_node.getValue().getGene().getGene_symbol() : "None";
			char strand = junc_node.getValue().getGene() != null ? junc_node.getValue().getGene().getStrand() : '.';
			List<IntRegion> list = new ArrayList<>(2);
			list.add(new Exon(null, junc_node.getStart() - 1, junc_node.getStart() - 1 + window_size));
			list.add(new Exon(null, junc_node.getEnd() - window_size, junc_node.getEnd()));
			return new PeakOut(chr, junc_node.getStart() - 1, junc_node.getEnd(), name, FisherTest.combinePvalue(p_value_list), strand, junc_node.getStart() - 1, junc_node.getEnd(),
					junc_node.getValue().getIPIDs().size(), 2, new int[]{start_len, end_len}, new int[] {junc_node.getStart() - 1, junc_node.getEnd() - end_len}, 
					junc_node.getValue().getConfidence(getIPSup(mapping_tree, junc_node.getStart(), junc_node.getEnd()),
					getIPSup(mapping_tree, junc_node.getStart(), junc_node.getEnd())), 0.0, getPropor(list, mapping_tree));
		}
		return null;
	}
	
	private List<Double> calCircPvalue(FisherTest fisher_test, String chr, int start, int end, int window_size, IntervalTree<MappingStat> mapping_tree){
		ArrayList<Double> p_value_list = new ArrayList<>();
		if (mapping_tree == null) {
			return p_value_list;
		}
		for (int i = start; i <= end; i += window_size) {
			p_value_list.add(calP_ValueBack(fisher_test, chr, i, i + window_size - 1, mapping_tree, InParam.getParams().isPeakBSJ() ? 'c' : 't', false));
			if (p_value_list.get(p_value_list.size() - 1) > InParam.getParams().getP_value()) {
				break;
			}
		}
		int start_len = p_value_list.size();
		if (start_len < 1 || p_value_list.get(start_len - 1) <= InParam.getParams().getP_value()) {
			return p_value_list;
		}
		for (int i = end; i >= start; i -= window_size) {
			double p_value = calP_ValueBack(fisher_test, chr, i - window_size + 1, i, mapping_tree, InParam.getParams().isPeakBSJ() ? 'c' : 't', false);
			if (p_value > InParam.getParams().getP_value()) {
				break;
			}
			else {
				p_value_list.add(start_len, p_value);
			}
		}
		return p_value_list;
	}
	
	private double getPropor(Collection<? extends IntRegion> list, IntervalTree<MappingStat> mapping_tree) {
		int ip_total = 1; //void be divided by 0
		int input_total = 1;
		for (IntRegion region : list) {
			for (Iterator<Node<MappingStat>> map_nodes = mapping_tree.overlappers(region.getStart() + 1, region.getEnd()); map_nodes.hasNext();) {
				Node<MappingStat> map_node = map_nodes.next();
				ip_total += map_node.getValue().getReadsIP();
				input_total += map_node.getValue().getReads();
			}
		}
		return (double) ip_total / (double) MappingStat.getIpReads() / (double) input_total * (double) MappingStat.getInputReads();
	}
	
	private HashMap<String, IntervalTree<Double>> calPvalueTree(FisherTest fisher_test, int window_size, HashMap<String, IntervalTree<MappingStat>> mapping_table){
		HashMap<String, IntervalTree<Double>> p_value_tree = new HashMap<>();
		if (mapping_table == null) {
			return p_value_tree;
		}
		List<String> chrs = new ArrayList<>(mapping_table.keySet());
		List<Double> p_value_list = new ArrayList<>();
		for (String chr : chrs) {
			p_value_tree.put(chr, new IntervalTree<>());
			IntervalTree<MappingStat> mapping_tree = mapping_table.get(chr);
			if (mapping_tree == null) {
				continue;
			}
			for (int i = 1; i <= chr_lengths.getOrDefault(chr, 0); i += window_size) {
				p_value_list.add(calP_ValueBack(fisher_test, chr, i, i + window_size - 1, mapping_tree, 
						InParam.getParams().isPeakBSJ() ? 'l' : 't', InParam.getParams().isPeakQuality()));
			}
		}
		adjustPvalue(p_value_list, InParam.getParams().getAdjsut_p());
		int used_len = 0;
		int chr_index = 0;
		for (int i = 0; i < p_value_list.size() && chr_index < chrs.size();) {
			String chr = chrs.get(chr_index);
			if (i * window_size - used_len >= chr_lengths.getOrDefault(chr, 0)) {
				++chr_index;
				used_len = i * window_size;
				continue;
			}
			if (p_value_list.get(i) <= InParam.getParams().getP_value()) {
				p_value_tree.get(chr).put(i * window_size - used_len + 1, i * window_size + window_size - used_len, p_value_list.get(i));
			}
			++i;
		}
		return p_value_tree;
	}

	private void adjustPvalue(List<Double> p_value, String method){
		if (p_value == null || p_value.size() <= 1) {
			return;
		}
		
		if ("bon".equals(method)) {
			double n = (double) p_value.size();
			for (int i = 0; i < p_value.size() ; ++i) {
				p_value.set(i, p_value.get(i) * n);
			}
		}
		else if ("bh".equals(method)) {
			List<double[]> out = new ArrayList<>();
			for (Double d : p_value) {
				out.add(new double[] {d, 0.0});
			}
			List<double[]> sort = new ArrayList<>(out);
			sort.sort((p1, p2) -> {
				return Double.compare(p1[0], p2[0]);
			});
			double value = Double.MAX_VALUE;
			for (int i = sort.size(); i > 0 ; --i) {
				sort.get(i - 1)[1] = Math.min(value, sort.get(i - 1)[0] * (double) sort.size() / (double) i);
			}
			for (int i = 0; i < out.size(); i++) {
				p_value.set(i, out.get(i)[1]);
			}
		}
		else {
			System.out.println("Warning: unkown method of adjusting p-value, adjusting is disabled");
		}
	}
	
	private static double calP_ValueBack(FisherTest fisher_test, String chr, int start, int end, IntervalTree<MappingStat> mapping_tree, char method, boolean mapping){
		double p_value = 1.0;
		int ip_circ = 0;
		int ip_total = 0;
		int input_circ = 0;
		int input_total = 0;
		long ip_reads = mapping ? MappingStat.getMappingIpReads() : MappingStat.getIpReads();
		long input_reads = mapping ? MappingStat.getMappingInputReads() : MappingStat.getInputReads();
//		if (m.debug_flag) {
//			System.out.printf("Total Input Reads:\t%d\tMapping Reads:\t%d\nTotal IP Reads:\t%d\tMapping Reads:\t%d\nWindow:\t%d\t%d\n", 
//					MappingStat.getInputReads(), MappingStat.getMappingInputReads(), MappingStat.getIpReads(), MappingStat.getMappingIpReads(), start, end);
//		}
		Iterator<Node<MappingStat>> nodes = mapping_tree.overlappers(start, end);
		while (nodes.hasNext()){
			Node<MappingStat> node = nodes.next();
//			if (m.debug_flag) {
//				System.out.printf("Node:\t%d\t%d\tCN: %d\tCI: %d\tTN: %d\tMN: %d\tTI: %d\tMI: %d\n", node.getStart(), node.getEnd(), node.getValue().getCircReads(), 
//						node.getValue().getCircReadsIP(), node.getValue().getReads(), node.getValue().getMapReads(), node.getValue().getReadsIP(), node.getValue().getMapReadsIP());
//			}
			ip_circ += node.getValue().getCircReadsIP();
			input_circ += node.getValue().getCircReads();
			ip_total += mapping ? node.getValue().getMapReadsIP() : node.getValue().getReadsIP();
			input_total += mapping ? node.getValue().getCircReadsIP() : node.getValue().getReads();
		}
		switch (method) {
		case 't':
			p_value = fisher_test.calpValue(ip_total, input_total, ip_reads - ip_total,input_reads - input_total, 2);
			break;
			
		case 'c':
			p_value = fisher_test.calpValue(ip_circ, input_circ, ip_reads - ip_circ, input_reads - input_circ, 2);
			break;
			
		case 'l':
			p_value = fisher_test.calpValue(ip_total - ip_circ, input_total - input_circ, ip_reads - ip_total + ip_circ,
					input_reads - input_total + input_circ, 2);
			break;
			
		default:
			System.err.println("Unknown method for peak-calling: " + method);
			break;
		}
//		if (!m.search_flag) {
//			try (BufferedWriter writer = new BufferedWriter(new FileWriter(new File (InParam.getParams().getOut_prefix() + "_windows.bed"), true))){
//				writer.write(chr);
//				writer.write('\t');
//				writer.write(String.valueOf(start));
//				writer.write('\t');
//				writer.write(String.valueOf(end));
//				writer.write('\t');
//				writer.write(String.valueOf(ip_total));
//				writer.write('\t');
//				writer.write(String.valueOf(input_total));
//				writer.write('\t');
//				writer.write(String.valueOf(ip_reads));
//				writer.write('\t');
//				writer.write(String.valueOf(input_reads));
//				writer.write('\t');
//				writer.write(String.valueOf(p_value));
//				writer.newLine();
//				writer.flush();
//				writer.close();
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
//		}
		return p_value;
	}
	
	private String findGeneSymbol(int start, int end, IntervalTree<Gene> gene_tree) {
		if (gene_tree == null || start > end) {
			return "None";
		}
		Gene gene = null;
		int overlap = -1;
		int gene_len = -1;
		for (Iterator<Node<Gene>> gene_nodes = gene_tree.overlappers(start, end); gene_nodes.hasNext();) {
			Node<Gene> gene_node = gene_nodes.next();
			int overlap_len = Math.min(end, gene_node.getEnd()) - Math.max(start, gene_node.getStart());
			if (overlap_len > overlap) {
				gene = gene_node.getValue();
				gene_len = gene.getLength();
				overlap = overlap_len;
			}
			else if (overlap_len == overlap && gene_len < gene_node.getValue().getLength()) {
				gene = gene_node.getValue();
				gene_len = gene.getLength();
			}
		}
		return gene == null ? "None" : gene.getGene_symbol();
	}
	
	private void adjustP_Value(List<PeakOut> p_value, String method){
		if ("bon".equals(method)) {
			double n = (double) p_value.size();
			for (int i = 0; i < p_value.size() ; ++i) {
				p_value.get(i).setFdr(p_value.get(i).getScore() * n / (i + 1));
			}
		}
		else if ("bh".equals(method)) {
			List<PeakOut> sort = new ArrayList<>(p_value);
			sort.sort((p1, p2) -> {
				Double d = p1.getScore();
				return d.compareTo(p2.getScore());
			});
			double value = Double.MAX_VALUE;
			for (int i = sort.size(); i > 0 ; --i) {
				value = Math.min(value, sort.get(i - 1).getScore() * (double) sort.size() / (double) i);
				sort.get(i - 1).setFdr(value);
			}
		}
		else {
			System.out.println("Warning: unkown method of adjusting p-value, adjusting is disabled");
		}
	}

	private List<PeakOut> filtPeak(List<PeakOut> in) {
		List<PeakOut> out = new ArrayList<>();
		for (PeakOut peak : in) {
			if (isReachPvalue(peak.getScore())) {
				out.add(peak);
			}
		}
		return out;
	}
	
	private boolean isReachPvalue(double score) {
		switch (InParam.getParams().getCombine_p()) {
		case "ave":
		case "fm":
			return score <= InParam.getParams().getP_value();

		case "log":
		default:
			return score >= -Math.log(InParam.getParams().getP_value());
		}
	}
	
	/**
	 * overwrite a string list into a file
	 * @param fileName the file name/path to output a list of string;
	 * @param out the content to be output;
	 */
	public static <T> void fileWrite(String fileName, Collection<T> out){
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(new File (fileName)));){
			writeList(writer, out);
		}
		catch(IOException e) {
			e.printStackTrace();
		}
	}

	public static <T> void fileWriteWithHeader(String fileName, Collection<T> out, String header){
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(new File (fileName)));){
			writer.write(header);
			writer.newLine();
			writeList(writer, out);
		}
		catch(IOException e) {
			e.printStackTrace();
		}
	}
	
	public static <T> void writeList(BufferedWriter writer, Collection<T> out) throws IOException {
		for (T line : out) {
			writer.write(line.toString());
			writer.newLine();
		}
		writer.flush();
	}

	public static <T> void fileWrite(String fileName, HashMap<String, IntervalTree<T>> out){
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(new File (fileName)));){
			Bed3.writeBedInOrder(writer, out);
		}
		catch(IOException e) {
			e.printStackTrace();
		}
	}
	
	public static <T> void writeTreeMap(BufferedWriter writer, HashMap<String, IntervalTree<T>> out) throws IOException {
		for (Entry<String, IntervalTree<T>> entry : out.entrySet()) {
			Iterator<Node<T>> nodes = entry.getValue().overlappers(Integer.MIN_VALUE, Integer.MAX_VALUE);
			while (nodes.hasNext()) {
				Node<T> node = nodes.next();
				writer.write(entry.getKey());
				writer.write(Bed3.getSep());
				writer.write(node.getStart() - 1);
				writer.write(Bed3.getSep());
				writer.write(node.getEnd());
				writer.write(Bed3.getSep());
				writer.write(node.getValue().toString());
				writer.newLine();
			}
			writer.flush();
		}
	}
	
	/**
	 * append a string list into a file
	 * @param fileName the file name to output;
	 * @param out the content to be put out;
	 */
	public static <T> void fileAppend(String fileName, Collection<T> out){
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(new File(fileName), true))){
			synchronized (writer) {
				writeList(writer, out);
			}
		}
		catch(IOException e) {
			e.printStackTrace();
		}
	}
	
	
	/**
	 * give the index of the first integer which is not less than the target and return -1 while the list is null, return 0 while the list is empty
	 * @param target the target number, we want to search the first number not less than it
	 * @param inc_seq an increasing sequence include integers
	 * @return the index of the integer
	 */
	public static int searchMinNoLess(int target, ArrayList<Integer> inc_seq) {
		int out = -1;
		if (inc_seq == null) {
			return -1;
		}
		int l = 0;
		int r = inc_seq.size() - 1;
		if (inc_seq.size() == 0||target <= inc_seq.get(0)) {
			out = 0;
		}
		else if(target <= inc_seq.get(r)){
			int m = 0;
			while (l < r) {
				m = (l + r) >> 1;
				if (l == m) {
					out = r;
					break;
				}
				if (target < inc_seq.get(m)) {
					r = m;
				}
				else if (target > inc_seq.get(m)) {
					l = m;
				}
				else {
					while (target == inc_seq.get(m)) {
						out = m;
						m--;
					}
					break;
				}
			}
		}
		else {
			out = r + 1;
		}
		return out;
	}
	
	public static <T> String listToString(List<T> list, String sep) {
		if (list == null) {
			return null;
		}
		if (list.size() == 0) {
			return "";
		}
		StringBuffer sb = new StringBuffer();
		Iterator<T> it = list.iterator();
		sb.append(it.next());
		while (it.hasNext()) {
			sb.append(sep);
			sb.append(it.next());
		}
		return sb.toString();
	}
	
	public static void toSamFile(String file_name) {
		try (SamReader reader = SamReaderFactory.makeDefault().open(new File(file_name));
				BufferedWriter writer = new BufferedWriter(new FileWriter(file_name.replace(".bam", ".sam")))){
			writer.write(reader.getFileHeader().getSAMString());
			for (SAMRecord record : reader) {
				writer.write(record.getSAMString());
			}
			writer.flush();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public int[] adjustToClip(String chr, int start, int end, String AG_str, String GT_str) {
		if (!chr_seqs.containsKey(chr) || isExonBoundary(chr, start, end)) {
			return new int[]{start, end};
		}
		for (int i = 0; i < GT.length; i++) {
			String clip_GT = (GT_str.length() > getDev() ? GT_str.substring(GT_str.length() - getDev(), GT_str.length()) : GT_str) + GT[i];
			String clip_AG = AG[i] + (AG_str.length() > getDev() ? AG_str.substring(0, getDev()) : AG_str);
			int max_start_adjust = Math.min(start - 3, getDev());
			int adjust_start = getChrSeq(chr, start - 3 - max_start_adjust, start + getDev() - 3 + clip_AG.length()).
					toUpperCase().indexOf(clip_AG.toUpperCase()) - max_start_adjust;
			int adjust_end = getChrSeq(chr, end - getDev() + 2 - clip_GT.length(), 
					Math.min(getChrLength(chr), end + 2 + getDev())).toUpperCase().indexOf(clip_GT.toUpperCase()) - getDev();
			if (adjust_start + max_start_adjust >= 0 && adjust_end + getDev() >= 0) {
				return new int[] {start + adjust_start, end + adjust_end};
			}
			clip_GT = clip_GT.substring(0, clip_GT.length() - 2) + GT_neg[i];
			clip_AG = AG_neg[i] + clip_AG.substring(2);
			adjust_start = getChrSeq(chr, start - 3 - max_start_adjust, start + getDev() - 3 + clip_AG.length()).
					toUpperCase().indexOf(clip_AG.toUpperCase()) - max_start_adjust;
			adjust_end = getChrSeq(chr, end - getDev() + 2 - clip_GT.length(), 
					Math.min(getChrLength(chr), end + 2 + getDev())).toUpperCase().indexOf(clip_GT.toUpperCase()) - getDev();
			if (adjust_start + max_start_adjust >= 0 && adjust_end + getDev() >= 0) {
				return new int[] {start + adjust_start, end + adjust_end};
			}
		}
		return null;
	}
	
	public boolean isExonBoundary(String chr, int start, int end) {
		return exon_boundary_map.containsKey(chr) && exon_boundary_map.get(chr).containsKey(start) && (exon_boundary_map.get(chr).get(start) & 1) != 0
				&& exon_boundary_map.get(chr).containsKey(end) && (exon_boundary_map.get(chr).get(end) & 2) != 0;
	}
	
	public static String reverseFliq(String bases) {
		if (bases == null || bases.length() == 0) {
			return bases;
		}
		StringBuilder sb = new StringBuilder();
		for (int i = bases.length() - 1; i >= 0; i--) {
			switch (bases.charAt(i)) {
			case 'A':
				sb.append('T');
				break;
			case 'T':
				sb.append('A');
				break;	
			case 'G':
				sb.append('C');
				break;
			case 'C':
				sb.append('G');
				break;	
			default:
				sb.append(bases.charAt(i));
				break;
			}
		}
		return sb.toString();
	}
	
	public boolean checkAllGene(HashMap<String, IntervalTree<ArrayList<Gene>>> gene_info) {
		boolean out = true;
		for (Entry<String, IntervalTree<ArrayList<Gene>>> gene_lists : gene_info.entrySet()) {
			Iterator<Node<ArrayList<Gene>>> nodes = gene_lists.getValue().iterator();
			while (nodes.hasNext()) {
				Node<ArrayList<Gene>> node = nodes.next();
				ArrayList<Gene> gene_list = node.getValue();
				for (Gene gene : gene_list) {
					out &= gene.checkPos();
				}
			}
		}
		return out;
	}
	
	public static int getIPSup(IntervalTree<MappingStat> mapping_tree, int start, int end) {
		if (InParam.getParams().getIp_sup() < 0) {
			double out = 0.0;
			for (Iterator<Node<MappingStat>> start_nodes = mapping_tree.overlappers(start - getDev(), start + getDev()); 
					start_nodes.hasNext();) {
				Node<MappingStat> start_node = start_nodes.next();
				if (start_node.getStart() < start - getDev()) {
					continue;
				}
				out += start_node.getValue().getReadsIP();
			}
			for (Iterator<Node<MappingStat>> end_nodes = mapping_tree.overlappers(end - getDev(), end + getDev()); 
					end_nodes.hasNext();) {
				Node<MappingStat> end_node = end_nodes.next();
				if (end_node.getEnd() > end + getDev()) {
					continue;
				}
				out += end_node.getValue().getReadsIP();
			}
			return (int) Math.max(out * 0.01, InParam.getParams().getSup_reads());
		}
		return InParam.getParams().getIp_sup();
	}
	
	public static int nearestOne(int target, int a, int b) {
		return Math.abs((long) target - (long) a) < Math.abs((long) target - (long) b) ? a : b;
	}
	
	public static boolean withIn(int target, int lower, int upper) {
		if (lower > upper) {
			lower ^= upper;
			upper ^= lower;
			lower ^= upper;
		}
		return target >= lower && target <= upper;
	}
	
	/**
	 * print some tips with time
	 * @param prefix the tip of time prefix;
	 */
	public static void printNow(String prefix) {
		Date time = new Date();
		System.out.printf("%s %tF %tT\n", prefix, time, time);
	}
	
	public static Method getInstance() {
		return m;
	}
	
	public static RemoverRNA rRNA() {
		return m.r;
	}
	
	public static String getChrSeq(String chr, int start, int end) {
		if (m.chr_seqs.get(chr) == null || end < start) {
			return null;
		}
		return m.chr_seqs.get(chr).substring(start, end);
	}
	
	public static int getChrLength(String chr) {
		return m.chr_lengths.getOrDefault(chr, 0);
	}
	
	public static int getDev() {
		return InParam.getParams().getRead_dev();
	}
	
	public static boolean isIP() {
		return m.ip_flag;
	}
	
	
}
