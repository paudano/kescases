package edu.gatech.kesmlst;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import edu.gatech.kanalyze.comp.reader.FileSequenceSource;
import edu.gatech.kanalyze.comp.reader.SequenceSource;
import edu.gatech.kanalyze.condition.StreamConditionListener;
import edu.gatech.kanalyze.module.count.CountModule;
import edu.gatech.kanalyze.util.kmer.KmerUtil;
import edu.gatech.kestrel.KestrelConstants;
import edu.gatech.kestrel.activeregion.ActiveRegionContainer;
import edu.gatech.kestrel.activeregion.ActiveRegionDetector;
import edu.gatech.kestrel.activeregion.Haplotype;
import edu.gatech.kestrel.activeregion.RegionHaplotype;
import edu.gatech.kestrel.counter.CountMap;
import edu.gatech.kestrel.counter.IkcCountMap;
import edu.gatech.kestrel.counter.MemoryCountMap;
import edu.gatech.kestrel.io.InputSample;
import edu.gatech.kestrel.refreader.ReferenceReader;
import edu.gatech.kestrel.refreader.ReferenceRegion;
import edu.gatech.kestrel.refreader.ReferenceRegionContainer;
import edu.gatech.kestrel.refreader.ReferenceSequence;
//import edu.gatech.kestrel.runner.KestrelRunner.KestrelRunnerException;
import edu.gatech.kestrel.varfilter.VariantFilter;
import edu.gatech.kestrel.varfilter.VariantFilterRunner;
import edu.gatech.kestrel.variant.VariantCall;
import edu.gatech.kestrel.variant.VariantCaller;
import edu.gatech.kestrel.writer.VariantWriter;
import edu.gatech.kestrel.writer.VariantWriterInitException;

public class KesMLST {

	public static void main(String[] args) {
		
		// TODO Auto-generated method stub
		KmerUtil kUtil;  // k-mer utility
		
		ReferenceRegionContainer refRegionContainer;            // Reference sequences regions
		Iterator<ReferenceSequence> refRegionSequenceIterator;  // Iterator for references in refRegionContainer
		
		ReferenceSequence refSequence;     // A reference sequence from refRegionSequenceIterator
		ReferenceRegion[] refRegionArray;  // Reference regions associated with refSequence
		
		CountMap counter;             // K-mer counter
		
		ActiveRegionDetector arDetector;        // Detects active regions
		ActiveRegionContainer regionContainer;  // A container of active regions
		
		VariantCaller varCaller;  // Calls variants from alignments
		VariantFilterRunner variantFilterRunner;
		
		VariantWriter variantWriter;  // Writes variants to output
		
		SampleCall samCall; // Sample call 
		
		
		int flankLength;  // Length of flanks to add to regions extracted from reference sequences
		
		// Init
		//logger.trace("exec(): Started");
		//errThrowable = null;
		
		int kSize = 31;
		int kMinSize = 15;
		int kMinMask = 0; 
		boolean removeIkc = false;
		String tempDirName = ".";
		ArrayList<InputSample> sampleList = new ArrayList<>();
		
		for (int index = 1; index < args.length ; index++)
			sampleList.add(new InputSample(null, new SequenceSource[] {new FileSequenceSource(new File(args[index]), "auto", null, 1, null)}));
		
		//Add input fastq files to sampleList
		//uncomment when accepting arguements
		try {
			
			// Set flank length
			//flankLength = this.flankLength;
			
			//if (flankLength < 0) {
			//flankLength =  77; // (int) (kSize * DEFAULT_FLANK_LENGTH_MULTIPLIER);
			
			//	if (flankLength < 0)
			//		flankLength = Integer.MAX_VALUE;
			//}
			
			// Normalize temporary directory
			/*if (tempDirName == null)
				tempDirName = "";
			
			tempDirName = tempDirName.trim();
			
			if (tempDirName.isEmpty())
				tempDirName = ".";*/
			
			// Setup counter and kUtil
			/*if (kmerCountInMemory) {
				kUtil = KmerUtil.get(kSize);
				
				logger.info("Counting k-mers in memory");
				counter = new MemoryCountMap(kUtil, getCountModule());
				
			} else {
				kUtil = KmerUtil.get(kSize, kMinSize, kMinMask);
				
				logger.info("Counting k-mers from file");
				counter = new IkcCountMap(kUtil, getCountModule(), new File(tempDirName), removeIkc);
			}*/
			
			samCall = new SampleCall();
			kUtil = KmerUtil.get(kSize, kMinSize, kMinMask);
			counter = new IkcCountMap(kUtil, getCountModule(), new File(tempDirName), removeIkc);
			
			
			// Create filter runner
			variantFilterRunner = new VariantFilterRunner();
			variantFilterRunner.addFilter(VariantFilter.getFilter("coverage:0.5,0", null));
			
			// Create aligner and variant caller
			varCaller = new VariantCaller();
			//varCaller.setCallAmbiguousVariant(callAmbiguousVariant);
			
//			if (variantCallByRegion)
//				varCaller.setVariantCallByRegion();
//			else
//				varCaller.setVariantCallByReference();
//			
			// Read reference sequences
			//logger.info("Reading references");
			File referenceFile = new File(args[0]);
			ReferenceReader referenceReader = new ReferenceReader(kUtil);
			//referenceReader.setFlankLength(flankLength);
			//referenceReader.setRemoveDescription(removeReferenceSequenceDescription);
			//referenceReader.setRevComplementNegStrand(reverseComplementNegativeStrand);
			ArrayList<SequenceSource> referenceList = new ArrayList<>();
			referenceList.add(new FileSequenceSource(referenceFile, "auto", null, 1, null));
			try {
//				refRegionContainer = referenceReader.read(
//						referenceList.toArray(new SequenceSource[0]),
//						((intervalContainer.isEmpty()) ? null : intervalContainer.getMap())
//				);
				refRegionContainer = referenceReader.read(
						referenceList.toArray(new SequenceSource[0]),
						null
				);
				
			} catch (IOException ex) {
				err("Error reading reference sequence(s)", ex);
				
				return;
			}
			
			// Check reference sequences
			if (refRegionContainer.isEmpty()) {
				err("No reference sequences (see -r option)", null);
				
				return;
			}
			
			// Get k-mer counter and active region detector
			arDetector = new ActiveRegionDetector(counter);
			
			//arDetector.setAlignmentWeight(alignmentWeight);
//			arDetector.setMinimumDifference(minimumDifference);
//			arDetector.setDifferenceQuantile(differenceQuantile);
//			arDetector.setAnchorBothEnds(anchorBothEnds);
//			arDetector.setCountReverseKmers(countReverseKmers);
//			arDetector.setPeakScanLength(peakScanLength);
//			arDetector.setEndScanLimitFactor(endScanLimitFactor);
//			arDetector.setCallAmbiguousRegions(callAmbiguousRegions);
//			arDetector.setDecayMinimum(expDecayMin);
//			arDetector.setDecayAlpha(expDecayAlpha);
//			arDetector.setMaxAlignerState(maxAlignerState);
//			
			// Open output
//			try {
//				variantWriter = VariantWriter.getWriter(outputFormat, outputFile, variantCallByRegion, loader);
//				
//			} catch (IllegalArgumentException ex) {
//				err("Error opening variant writer: " + ex.getMessage(), ex);
//				return;
//				
//			} catch (FileNotFoundException ex) {
//				err("File not found while variant writer: " + ex.getMessage(), ex);
//				return;
//				
//			} catch (IOException ex) {
//				err("IO error while variant writer: " + ex.getMessage(), ex);
//				return;
//				
//			} catch (VariantWriterInitException ex) {
//				err("Error loading writer: " + ex.getMessage(), ex);
//				return;
//			}
//			
			// Process samples
			for (InputSample sample : sampleList) {
				//logger.info("Processing sample: {}", sample.name);
				
				samCall.reset();
				//variantWriter.setSampleName(sample.name);
				
				// Get counts
				try {
					counter.set(sample);
					
				} catch (FileNotFoundException ex) {
					err(null, ex);
					return;
				
				} catch (IOException ex) {
					err(null, ex);
					return;
				}
				
				// Find active regions
				//logger.trace("Getting active regions: {}", sample.name);
				
				// Iterate over reference sequences
				refRegionSequenceIterator = refRegionContainer.refSequenceIterator();
				
				while (refRegionSequenceIterator.hasNext()) {
					
					refSequence = refRegionSequenceIterator.next();
					
					refRegionArray = refRegionContainer.get(refSequence);
					
					// Iterate over reference regions
					for (ReferenceRegion refRegion : refRegionArray) {
						
						int varCount = 0;
						
						//	logger.trace("Searching for active regions in {}", refRegion);
							
						//variantWriter.setReferenceRegion(refRegion);
						
						regionContainer = arDetector.getActiveRegions(refRegion);
						
						// Find variants in active regions
						//for (int index=0 ; index <= regionContainer.haplotypes.length; index++){
						
						for (RegionHaplotype thisRegionHaplotype : regionContainer.haplotypes) {
							
							
							//logger.trace("Found: {} ({} haplotypes)", thisRegionHaplotype.toString(), thisRegionHaplotype.haplotype.length);
							varCaller.init(thisRegionHaplotype.activeRegion);
							// Find variants in each haplotype
								for (Haplotype haplotype : thisRegionHaplotype.haplotype) {
									varCaller.add(haplotype);
							
							
							// Output variant
								for (VariantCall var : varCaller.getVariants()) {
									var = variantFilterRunner.filter(var);
								
									if (var != null){
										varCount += 1;
										
								}
							}
						}
					}
						
					samCall.addRef(regionContainer, varCount);	// Update sample call
				}
			}
				
				//System.out.println(sample.name);
				samCall.print();
				System.out.println();
			}
			//logger.trace("exec(): Complete");
			
		} catch (Exception ex) {
		//	logger.error("Unexpected exception in run(): {} ({})", ex.getMessage(), ex.getClass().getSimpleName());
		//	errThrowable = ex;
			
			ex.printStackTrace();
		}
		
		return;
	}

	private static CountModule getCountModule() {
		CountModule countModule = new CountModule();
		
		countModule.configure(null);
		int kSize = 31;
		int minKmerCount = 0;
		// Add libraries to KAnalyze
//		for (URL libUrl : libraryUrlList)
//			countModule.addLibraryURL(libUrl);
		
		// Set k-mer size
		countModule.setKSize(kSize);
		
		// Set options
//		countModule.setSequenceBatchSize(sequenceBatchSize);
//		countModule.setSequenceQueueSize(sequenceQueueSize);
//		countModule.setKmerBatchSize(kmerBatchSize);
//		countModule.setKmerQueueSize(kmerQueueSize);
//		countModule.setSplitThreadCount(splitThreadCount);
//		countModule.setKmerThreadCount(kmerThreadCount);
//		countModule.setKmerCacheSize(batchCacheSize);
//		countModule.setSequenceCacheSize(batchCacheSize);
//		countModule.setSegmentSize(segmentSize);
//		countModule.setTempDirName(tempDirName);
//		
		if (minKmerCount > 0)
			countModule.addPostCountFilterDefinition("kmercount:" + minKmerCount);
		
		//countModule.setFreeSegment(freeResources);
		
		countModule.addListener(new StreamConditionListener("Kestrel", System.err, null, true));  // TODO: Replace with a Kestrel condition listener
		
		return countModule;
	}
	private static void err(String message, Throwable ex) {
		
		// Log
		System.err.println("Error: " + message + ((ex != null) ? ": " + ex.getMessage() : ""));
		System.exit(KestrelConstants.ERR_SYSTEM);
		return;
	}
}

class SampleCall {
	
	private HashMap<String, GeneCall> callMap;
	
	public SampleCall() {
		callMap = new HashMap<>();
		
		return;
	}
	
	public void addRef(ActiveRegionContainer region, int varCount) {
		String[] tok = region.refRegion.name.split("_", 2);
		
		GeneCall call = callMap.get(tok[0]);
		
		if (call == null) {
			call = new GeneCall(0, Integer.MAX_VALUE);
			
			callMap.put(tok[0], call);
		}
		
		if (region.stats.min == 0 ) 
			return;
		
		if (varCount < call.distance)
			callMap.put(tok[0], new GeneCall(Integer.parseInt(tok[1]), varCount));
		
		return;
	}
	
	public void print() {
		System.out.println("Gene\tAllele\tDistance");

		
		for (String geneName : callMap.keySet()) {
			GeneCall call = callMap.get(geneName);
			
			if( call.allele != 0)
				System.out.printf("%s\t%d\t%d\n", geneName, call.allele, call.distance);
			else
				System.out.printf("%s\tNA\tNA\n", geneName);
		}
		
		return;
	}
	
	public void reset() {
		callMap.clear();
	}
}

class GeneCall {
	
	public final int allele;
	
	public final int distance;
	
	public GeneCall(int allele, int distance) {
		
		this.allele = allele;
		this.distance = distance;

		return;
	}
	
}
