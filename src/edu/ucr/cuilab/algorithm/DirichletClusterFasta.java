package edu.ucr.cuilab.algorithm;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

public class DirichletClusterFasta {

	private static List<String> getPermutations(char[] charList, int depth) {
		List<String> stringList = new ArrayList<String>();
		if (1 == depth) {
			for (int i = 0; i < charList.length; i++) {
				stringList.add(String.valueOf(charList[i]));
			}
		} else {
			List<String> subStringList = getPermutations(charList, depth - 1);
			for (int i = 0; i < charList.length; i++) {
				for (int j = 0; j < subStringList.size(); j++) {
					stringList.add(String.valueOf(charList[i])
							+ subStringList.get(j));
				}
			}
		}
		return stringList;
	}

	private static double parse(double defaultValue, String input,
			String descript) {
		double result = defaultValue;
		try {
			result = Double.parseDouble(input);
		} catch (NullPointerException npe) {
			System.out.println(descript + " error! Use default value "
					+ defaultValue);
			result = defaultValue;
		} catch (NumberFormatException nfe) {
			System.out.println(descript + " error! Use default value "
					+ defaultValue);
			result = defaultValue;
		}
		return result;
	}

	private static int parse(int defaultValue, String input, String descript) {
		int result = defaultValue;
		try {
			result = Integer.parseInt(input);
		} catch (NullPointerException npe) {
			System.out.println(descript + " error! Use default value "
					+ defaultValue);
			result = defaultValue;
		} catch (NumberFormatException nfe) {
			System.out.println(descript + " error! Use default value "
					+ defaultValue);
			result = defaultValue;
		}
		return result;
	}

	public static void frontEnd(String[] args) throws Exception {

		int transOrder = DefaultConstants.TRANSORDER;
		int particles = DefaultConstants.PARTICLES;
		int neighbor = DefaultConstants.NEIGHBOR;
		int seqs = 0;
		double alpha = DefaultConstants.ALPHALOW;
		double alphaLow = DefaultConstants.ALPHALOW;
		double alphaHigh = DefaultConstants.ALPHAHIGH;
		double threshold = DefaultConstants.THRESHOLD;
		double zero = DefaultConstants.ZERO;
		double majority = DefaultConstants.MAJORITY;
		int clusterNumForTest = DefaultConstants.TESTNUM;
		double coverage = DefaultConstants.COVERAGE;
		int neighborX = 1;
		// String input = "/home/xinping/Desktop/6008/5_test.txt";
		// String input = "/home/xinping/Desktop/6008/5_2345Borr.txt";
		// String input = "/home/xinping/Desktop/6008/Data/5_3910AcaryEhrl.txt";
//		String input = "/home/xinping/Desktop/6008/Data/AcaryEhrlNew.txt";
		// String input = "/home/xinping/Desktop/6008/Data/AcaryEhrlNew_5W.txt";
		// String input =
		// "/home/xinping/Desktop/6008/Data/Acary_Ehrl_9104_pair_5W.txt";
//		String input = "/home/xinping/Desktop/6008/Data/678910_1W.fna";
		String input = "/home/xinping/options/5.fna";
		String output = null;

		Date date = new Date();

		Params params = new Params(particles, neighbor, seqs, transOrder,
				majority, threshold, alpha, alphaLow, alphaHigh, zero,
				clusterNumForTest, coverage, neighborX);

		if (0 == args.length) {
			printHelp();
//			params.setAlpha(0.0000000001528);
//			params.setAlphaLow(0.0000000001528);
			return;
		}

		int pos = 0;
		while (pos < args.length - 1) {
			if (args[pos].charAt(0) == '-') {
				switch (args[pos].charAt(1)) {
				case 'a':
				case 'A':
					params.setAlpha(parse(DefaultConstants.ALPHA,
							args[pos + 1], "Parameter alpha"));
					pos += 2;
					break;
				case 'c':
				case 'C':
					params.setCoverage(parse(DefaultConstants.COVERAGE,
							args[pos + 1], "Parameter coverage"));
					pos += 2;
					break;
				case 'i':
				case 'I':
					input = args[pos + 1];
					pos += 2;
					break;
				case 'l':
				case 'L':
					params.setAlphaLow(parse(DefaultConstants.ALPHALOW,
							args[pos + 1], "Parameter alpha_low"));
					pos += 2;
					break;
				case 'm':
				case 'M':
					params.setMajority(parse(DefaultConstants.MAJORITY,
							args[pos + 1], "Parameter majority"));
					pos += 2;
					break;
				case 'n':
				case 'N':
					params.setNeighbor(parse(DefaultConstants.NEIGHBOR,
							args[pos + 1], "Parameter neighbor_number"));
					pos += 2;
					break;
				case 'o':
				case 'O':
					params.setTransOrder(parse(DefaultConstants.TRANSORDER,
							args[pos + 1], "Parameter transorder"));
					pos += 2;
					break;
				case 'p':
				case 'P':
					params.setParticles(parse(DefaultConstants.PARTICLES,
							args[pos + 1], "Parameter particle_amount"));
					pos += 2;
					break;
				case 't':
				case 'T':
					params.setThreshold(parse(DefaultConstants.THRESHOLD,
							args[pos + 1], "Parameter threshold"));
					pos += 2;
					break;
				case 'u':
				case 'U':
					output = args[pos + 1];
					pos += 2;
					break;
				case 'x':
				case 'X':
					params.setNeighborX(parse(DefaultConstants.NEIGHBORX,
							args[pos + 1], "Parameter neighborX"));
					pos += 2;
					break;
				case 'z':
				case 'Z':
					params.setZero(parse(DefaultConstants.ZERO, args[pos + 1],
							"Parameter zero"));
					pos += 2;
					break;
				default:
					pos++;
					break;
				}
			} else {
				pos++;
			}
		}

		if (null == output) {
			output = input + "." + date.toString() + ".output";
		}

		mainJobFasta(input, output, params);
	}

	private static List<String[]> readFasta(String inputFile) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(new File(
				inputFile)));
		List<String[]> readList = new ArrayList<String[]>();
		String line = null;
		int counter = -1;
		while (null != (line = br.readLine())) {
			if (line.startsWith(">")) {
				readList.add(new String[] { line, "" });
				counter++;
			} else {
				readList.get(counter)[1] += line;
			}
		}
		br.close();
		return readList;
	}

	public static void printHelp() {
		System.out.println("Usage: ");
		System.out.println("\tjava -jar " + DefaultConstants.PACKAGENAME
				+ " -i abundance_species_equal.fna");
		System.out.println("The FASTA file contains pair-end sequence\n");
		System.out.println("More Usage: ");
		System.out.println("\tjava -jar " + DefaultConstants.PACKAGENAME
				+ " [OPTION] abundance_species_equal.fna\n");
		System.out.println("Options:");
		System.out
				.println("\t-a:\tThe value followed will be the alpha parameter");
		System.out.println("\t   \tDefault value " + DefaultConstants.ALPHA);
		System.out
				.println("\t-c:\tThe value followed will be the coverage parameter");
		System.out.println("\t   \tDefault value " + DefaultConstants.COVERAGE);
		System.out
				.println("\t-l:\tThe value followed will be the alpha_low parameter");
		System.out.println("\t   \tDefault value " + DefaultConstants.ALPHALOW);
		System.out
				.println("\t-m:\tThe value followed will be the majority parameter");
		System.out.println("\t   \tDefault value " + DefaultConstants.MAJORITY);
		System.out
				.println("\t-n:\tThe value followed will be the neighbor_num parameter");
		System.out.println("\t   \tDefault value " + DefaultConstants.NEIGHBOR);
		System.out
				.println("\t-o:\tThe value followed will be the transform_order parameter");
		System.out.println("\t   \tDefault value "
				+ DefaultConstants.TRANSORDER);
		System.out
				.println("\t-p:\tThe value followed will be particle_number parameter");
		System.out
				.println("\t   \tDefault value " + DefaultConstants.PARTICLES);
		System.out
				.println("\t-t:\tThe value followed will be the threshold parameter");
		System.out
				.println("\t   \tDefault value " + DefaultConstants.THRESHOLD);
		System.out
				.println("\t-u:\tThe string followed will be the output file");
		System.out
				.println("\t   \tDefault value is the input file name with additional time");
		System.out
				.println("\t-x:\tThe value followed will be neighborX parameter");
		System.out
				.println("\t   \tDefault value " + DefaultConstants.NEIGHBORX);
		System.out.println("\t-z:\tThe value followed will be zero parameter");
		System.out.println("\t   \tDefault value " + DefaultConstants.ZERO);
	}

	private static List<Set<Integer>> updateOverlap(
			List<DoubleRead> doubleReadList, List<Set<Integer>> overlapList) {
		Map<Integer, Integer> reverseMap = new HashMap<Integer, Integer>();
		for (int i = 0; i < doubleReadList.size(); i++) {
			reverseMap.put(doubleReadList.get(i).getId(), i);
		}
		List<Set<Integer>> resultList = new ArrayList<Set<Integer>>();
		for (int i = 0; i < overlapList.size(); i++) {
			Set<Integer> temp = new HashSet<Integer>();
			for (Integer id : overlapList.get(i)) {
				if (reverseMap.containsKey(id)) {
					temp.add(reverseMap.get(id));
				}
			}
			if (temp.size() > 0) {
				resultList.add(temp);
			}
		}
		Collections.sort(resultList, new ListCompareLength());
		return resultList;
	}

	private static int max(int[][] data) {
		int result = 0;
		for (int i = 0; i < data.length; i++) {
			for (int d : data[i]) {
				if (d > result) {
					result = d;
				}
			}
		}
		return result;
	}

	public static int[] dirichletClusterSingle(List<DoubleRead> doubleReadList,
			List<Set<Integer>> overlapList, Params params, int zModeLower) {

		int[][] z = new int[params.getSeqs()][params.getParticles()];
		double[] w = new double[params.getParticles()];
		for (int i = 0; i < params.getSeqs(); i++) {
			for (int j = 0; j < params.getParticles(); j++) {
				z[i][j] = 0;
			}
		}
		for (int i = 0; i < params.getParticles(); i++) {
			w[i] = 1.0 / (0.0 + params.getParticles());
		}

		Map<Integer, String> idTagMap = new HashMap<Integer, String>();
		for (int i = 0; i < params.getSeqs(); i++) {
			idTagMap.put(i, "");
		}

		List<Map<Integer, List<Integer>>> tempAccumCountList = new ArrayList<Map<Integer, List<Integer>>>(
				params.getParticles());
		List<Map<Integer, Integer>> mapCountList = new ArrayList<Map<Integer, Integer>>(
				params.getParticles());

		for (int i = 0; i < params.getParticles(); i++) {
			tempAccumCountList.add(new TreeMap<Integer, List<Integer>>());
			mapCountList.add(new TreeMap<Integer, Integer>());
		}
		for (Integer elem : overlapList.get(0)) {
			for (int j = 0; j < params.getParticles(); j++) {
				z[elem][j] = 1;
			}
			idTagMap.put(elem, idTagMap.get(elem) + "Init");
		}
		tempAccumCountList = DirichletClusterSingle.updateCountLists(
				doubleReadList, z, overlapList.get(0), tempAccumCountList);
		mapCountList = DirichletClusterSingle.updateGroupMapList(mapCountList,
				z, overlapList.get(0));
		// DirichletClusterSingle.newGroup(doubleReadList,
		// params.getNeighbor());
		DirichletClusterSingle.newGroup(doubleReadList, params.getSeqs()
				/ (max(z) + 1) / params.getNeighborX());

		List<Set<Integer>> rmList = new ArrayList<Set<Integer>>();
		int accumSeqCount = overlapList.get(0).size();
		int overlapInTotal = 0;
		Map<Integer, Set<Integer>> overlapResult = new HashMap<Integer, Set<Integer>>();
		Map<Integer, Set<Integer>> overallResult = new HashMap<Integer, Set<Integer>>();

		int group = max(z);
		for (int i = 1; i < overlapList.size(); i++) {
			if (overlapList.get(i).size() > 1) {

				overlapInTotal++;
				// update
				int[][] tempZLower = DirichletClusterSingle.clusterOverlapSeqs(
						mapCountList, tempAccumCountList, doubleReadList,
						overlapList.get(i), z, w, params.getParticles(),
						params.getAlphaLow(), accumSeqCount, group);
				int[][] tempZUpper = DirichletClusterSingle.clusterOverlapSeqs(
						mapCountList, tempAccumCountList, doubleReadList,
						overlapList.get(i), z, w, params.getParticles(),
						params.getAlphaHigh(), accumSeqCount, group);
				boolean isMarjorityLower = DirichletClusterSingle
						.checkMajorityVote(tempZLower, params.getMajority());
				boolean isMarjorityUpper = DirichletClusterSingle
						.checkMajorityVote(tempZUpper, params.getMajority());
				int[] voteLower = DirichletClusterSingle
						.getMode(DirichletClusterSingle.getMode(tempZLower));
				int[] voteUpper = DirichletClusterSingle
						.getMode(DirichletClusterSingle.getMode(tempZUpper));
				int majorityVoteLower = voteLower[0];
				int majorityVoteUpper = voteUpper[0];
				if (isMarjorityLower && isMarjorityUpper
						&& (majorityVoteLower == majorityVoteUpper)) {

					for (Integer j : overlapList.get(i)) {
						for (int k = 0; k < params.getParticles(); k++) {
							z[j][k] = majorityVoteLower;
						}
						idTagMap.put(j, idTagMap.get(j) + "Majority");
					}

					tempAccumCountList = DirichletClusterSingle
							.updateCountLists(doubleReadList, z,
									overlapList.get(i), tempAccumCountList);
					mapCountList = DirichletClusterSingle.updateGroupMapList(
							mapCountList, z, overlapList.get(i));

					accumSeqCount += overlapList.get(i).size();
					if (overlapResult.containsKey(majorityVoteLower)) {

						overlapResult.get(majorityVoteLower).addAll(
								overlapList.get(i));

					} else {
						overlapResult.put(majorityVoteLower,
								new HashSet<Integer>(overlapList.get(i)));

					}
					if (max(z) > group) {
						DirichletClusterSingle.newGroup(
								doubleReadList,
								params.getSeqs() / (max(z) + 1)
										/ params.getNeighborX());
					}
					group = max(z);

				} else {
					rmList.add(overlapList.get(i));
				}
			} else {

				for (Integer elem : overlapList.get(i)) {
					idTagMap.put(elem, idTagMap.get(elem) + "Single");
					// update
					z[elem] = DirichletClusterSingle.clusterOverlapSeqs(
							mapCountList, tempAccumCountList, doubleReadList,
							overlapList.get(i), z, w, params.getParticles(),
							params.getAlpha(), accumSeqCount, group)[0];
					w = DirichletClusterSingle.updateWeights(
							tempAccumCountList, doubleReadList, elem, z, w,
							params.getParticles(), group);

					double eff = 0.0;
					for (double d : w) {
						eff += d * d;
					}
					if (eff > 1 / (params.getThreshold() * params
							.getParticles())) {
						int[] resampleIndex = DirichletClusterSingle.sampleInt(
								w, params.getParticles(), 0);
						int[] temp = new int[params.getParticles()];
						for (int j = 0; j < params.getParticles(); j++) {
							w[j] = 1.0 / (0.0 + params.getParticles());
							temp[j] = z[elem][resampleIndex[j]];
						}
						z[elem] = temp;
					}

					tempAccumCountList = DirichletClusterSingle
							.updateCountLists(doubleReadList, z,
									overlapList.get(i), tempAccumCountList);
					mapCountList = DirichletClusterSingle.updateGroupMapList(
							mapCountList, z, overlapList.get(i));
				}
				if (max(z) > group) {
					DirichletClusterSingle.newGroup(
							doubleReadList,
							params.getSeqs() / (max(z) + 1)
									/ params.getNeighborX());
				}
				group = max(z);
				accumSeqCount += overlapList.get(i).size();
			}

		}

		for (Set<Integer> seqList : rmList) {
			// update
			int[][] tempZ = DirichletClusterSingle.clusterOverlapSeqs(
					mapCountList, tempAccumCountList, doubleReadList, seqList,
					z, w, params.getParticles(), params.getAlpha(),
					accumSeqCount, group);
			int majorityVote = DirichletClusterSingle
					.getMode(DirichletClusterSingle.getMode(tempZ))[0];
			for (Integer seqId : seqList) {
				idTagMap.put(seqId, idTagMap.get(seqId) + "Remove");
				for (int k = 0; k < params.getParticles(); k++) {
					z[seqId][k] = majorityVote;
				}
			}
			accumSeqCount += seqList.size();

			tempAccumCountList = DirichletClusterSingle.updateCountLists(
					doubleReadList, z, seqList, tempAccumCountList);
			mapCountList = DirichletClusterSingle.updateGroupMapList(
					mapCountList, z, seqList);
			if (max(z) > group) {
				DirichletClusterSingle
						.newGroup(doubleReadList, params.getSeqs()
								/ (max(z) + 1) / params.getNeighborX());
			}
			group = max(z);
		}

		int[] result = DirichletClusterSingle
				.compressResult(DirichletClusterSingle.getMode(z));
		for (int i = 0; i < result.length; i++) {
			if (overallResult.containsKey(result[i])) {
				overallResult.get(result[i]).add(i);
			} else {
				Set<Integer> temp = new TreeSet<Integer>();
				temp.add(i);
				overallResult.put(result[i], temp);
			}
		}
		for (int i = 0; i < result.length; i++) {
			result[i] = result[i] + zModeLower;
		}

		System.out.println();
		System.out.print(params.getAlphaHigh());
		System.out.print("\t");
		System.out.print(doubleReadList.size());
		System.out.print("\t");
		System.out.print(overlapInTotal);
		System.out.print("\t");
		System.out.print(overlapInTotal - rmList.size());
		System.out.print("\t");
		System.out.print(1.0 - (rmList.size() + 0.0) / overlapInTotal);
		System.out.print("\t");
		System.out.println(String.valueOf((1.0 - (rmList.size() + 0.0)
				/ overlapInTotal) > params.getCoverage()));
		System.out.println("Overlap");
		for (Integer key : overlapResult.keySet()) {
			System.out.print(String.valueOf(key) + ":");
			Map<String, Integer> tempMap = new HashMap<String, Integer>();
			for (Integer temp : overlapResult.get(key)) {
				
				String strKey = doubleReadList.get(temp).getGI();
				if (tempMap.containsKey(strKey)) {
					tempMap.put(strKey, tempMap.get(strKey) + 1);
				} else {
					tempMap.put(strKey, 1);
				}
			}
			for (String tempKey : tempMap.keySet()) {
				System.out.print("\t");
				System.out.print(tempKey);
				System.out.print("\t");
				System.out.print(tempMap.get(tempKey));
				System.out.print("\t");
				System.out.print(tempMap.get(tempKey)
						/ (overlapResult.get(key).size() + 0.0));
			}
			System.out.println();
		}
		System.out.println("Overall");
		for (Integer key : overallResult.keySet()) {
			System.out.print(String.valueOf(key) + ":");
			Map<String, Integer> tempMap = new HashMap<String, Integer>();
			for (Integer temp : overallResult.get(key)) {

				String strKey = doubleReadList.get(temp).getGI();
				if (tempMap.containsKey(strKey)) {
					tempMap.put(strKey, tempMap.get(strKey) + 1);
				} else {
					tempMap.put(strKey, 1);
				}
			}
			for (String tempKey : tempMap.keySet()) {
				System.out.print("\t");
				System.out.print(tempKey);
				System.out.print("\t");
				System.out.print(tempMap.get(tempKey));
				System.out.print("\t");
				System.out.print(tempMap.get(tempKey)
						/ (overallResult.get(key).size() + 0.0));
			}
			System.out.println();
			System.out.println();
		}
		if ((1.0 - (rmList.size() + 0.0) / overlapInTotal) > params
				.getCoverage()) {
			return result;
		} else {
			return null;
		}
	}

	private static int[] dirichletCluster(List<DoubleRead> doubleReadList,
			List<Set<Integer>> overlapList, int[] zMode, Params params)
			throws CloneNotSupportedException {
		int[] updateZMode = zMode.clone();
		int tempZModeLower = 0;
		Set<Integer> zModeSet = new TreeSet<Integer>();
		for (int i = 0; i < updateZMode.length; i++) {
			zModeSet.add(updateZMode[i]);
		}
		double[] alphaHighArr = new double[] { 10, 5, 2, 1, 0.5, 0.1, 0.01,
				0.001 };
		System.out
				.println((new Date()).toString()
						+ "\tDirichlet Cluster start with alpha = "
						+ params.getAlpha());

		int groupNum = zModeSet.size();
		int oldGroupNum = zModeSet.size();

		for (Integer elem : zModeSet) {

			System.out.println((new Date()).toString() + "\tzMode = " + elem);

			List<DoubleRead> subDoubleReadList = new ArrayList<DoubleRead>();
			for (int i = 0; i < zMode.length; i++) {
				if (elem.equals(zMode[i])) {
					subDoubleReadList.add(doubleReadList.get(i));
				}
			}

			System.out.println((new Date()).toString()
					+ "\tUpdating overlap infomation");

			List<Set<Integer>> subOverlapList = updateOverlap(
					subDoubleReadList, overlapList);
			Params tempParams = (Params) params.clone();
			tempParams.setSeqs(subDoubleReadList.size());
			// tempParams.setNeighbor(subDoubleReadList.size()/10);
			// tempParams.setNeighbor(decideNeighbor(subDoubleReadList.size()));
			if (1 < groupNum) {
				tempParams.setNeighbor(subDoubleReadList.size()
						/ (groupNum * params.getNeighborX()));
			}
			System.out.println((new Date()).toString()
					+ "\tDirichlet single round begin");

			int[] subResult = new int[subDoubleReadList.size()];
			for (int i = 0; i < subResult.length; i++) {
				subResult[i] = tempZModeLower + 1;
			}

			boolean flag = false;

			for (double alphaHigh : alphaHighArr) {
				tempParams.setAlphaHigh(alphaHigh);
				int[] tempSubResult = dirichletClusterSingle(subDoubleReadList,
						subOverlapList, tempParams, tempZModeLower);
				if (!flag && (null != tempSubResult)) {
					subResult = tempSubResult;
					flag = true;
					break;
				}
			}
			for (int i = 0; i < subResult.length; i++) {
				updateZMode[subDoubleReadList.get(i).getId()] = subResult[i];
				if (subResult[i] > tempZModeLower) {
					tempZModeLower = subResult[i];
				}

			}

			oldGroupNum--;
			groupNum = oldGroupNum + tempZModeLower;

			System.out.println((new Date()).toString()
					+ "\tDirichlet single round complete");

		}

		System.out.println((new Date()).toString()
				+ "\tDirichlet Cluster complete");

		return updateZMode;
	}

	private static String getGI(String description) {
		int start = description.indexOf("GI=") + 3;
		int end = description.indexOf(',');
		return description.substring(start, end);
	}

	public static void mainJobFasta(String inputFile, String outputFile,
			Params params) throws Exception {

		System.out.println((new Date()).toString()
				+ "\tProgram start... \nReading dataset...");

		List<String[]> readFastaList = readFasta(inputFile);
		List<String> readList = new ArrayList<String>(readFastaList.size());
		List<String> descList = new ArrayList<String>(readFastaList.size());
		params.setSeqs(readFastaList.size() / 2);

		for (int i = 0; i < readFastaList.size(); i++) {
			descList.add(readFastaList.get(i)[0]);
			readList.add(readFastaList.get(i)[1]);
		}

		char[] charList = { 'A', 'C', 'G', 'T' };
		List<String> permutationList = getPermutations(charList, 2);

		System.out.println((new Date()).toString()
				+ "\tBegin to calculate overlapping infos");

		List<Set<Integer>> overlapList = EncodeBinTask.encodeMainJob(readList,
				32);

		System.out.println((new Date()).toString() + "\tComplete");

		List<DoubleRead> doubleReadList = new ArrayList<DoubleRead>();
		for (int i = 0; i < params.getSeqs(); i++) {
			DoubleRead dr = new DoubleRead(i, descList.get(i * 2 + 0),
					descList.get(i * 2 + 1), readList.get(i * 2 + 0),
					readList.get(i * 2 + 1), 2, permutationList);
			dr.setGI(getGI(descList.get(i * 2 + 0)));
			doubleReadList.add(dr);
		}
		
		System.out
		.println((new Date()).toString() + "\tMain processor started");
		
		int[] zMode = new int[params.getSeqs()];
		for (int i = 0; i < params.getSeqs(); i++) {
			zMode[i] = 1;
		}
		zMode = dirichletCluster(doubleReadList, overlapList, zMode, params);

		System.out.println((new Date()).toString()
				+ "\tComplete\nWriting result...");

		PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFile)));
		for (int i = 0; i < zMode.length; i++) {
			pw.println(zMode[i]);
		}
		pw.close();
	}

	public static void main(String[] args) throws Exception {
		frontEnd(args);
	}
}
