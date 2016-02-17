package edu.ucr.cuilab.algorithm;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

public class DirichletClusterSingle {

	private static void newGroup(List<DoubleRead> doubleReadList, int neighbor) {
		Collections.sort(doubleReadList, new DoubleReadCompareGC());
		for (int i = 0; i < doubleReadList.size(); i++) {
			int nflag = 1;
			int pflag = 1;
			List<Integer> accumCountList = new ArrayList<Integer>(doubleReadList.get(i).getCountList());
//			System.out.println(accumCountList.toString());

			while ((nflag <= neighbor / 2) && (i - nflag >= 0)) {
				List<Integer> temp = doubleReadList.get(i - nflag)
						.getCountList();
//				System.out.println("COUNT");
//				System.out.println(temp.toString());
				for (int j = 0; j < accumCountList.size(); j++) {
					accumCountList.set(j, accumCountList.get(j) + temp.get(j));
				}
//				System.out.println("ACCUM");
//				System.out.println(accumCountList.toString());
				nflag++;
			}

			while ((pflag <= neighbor / 2)
					&& (i + pflag < doubleReadList.size())) {
				List<Integer> temp = doubleReadList.get(i + pflag)
						.getCountList();
//				System.out.println("COUNT");
//				System.out.println(temp.toString());
				for (int j = 0; j < accumCountList.size(); j++) {
					accumCountList.set(j, accumCountList.get(j) + temp.get(j));
				}
//				System.out.println("ACCUM");
//				System.out.println(accumCountList.toString());
				pflag++;
			}

			doubleReadList.get(i).setAccumulateCountList(accumCountList);
			doubleReadList.get(i).setNewGroupStartp(
					startPCounts(accumCountList));
			doubleReadList.get(i).setNewGroupTransp(
					transPCounts(accumCountList));
//			doubleReadList.get(i).outputStartp();
		}
		Collections.sort(doubleReadList, new DoubleReadCompareID());
//		doubleReadList.get(0).outputStartp();
//		doubleReadList.get(0).outputAccum();
//		return doubleReadList;
	}

	// return {mode, length of mode}
	private static int[] getMode(int[] datas) {
		Map<Integer, Integer> tempMap = new HashMap<Integer, Integer>();
		for (int i : datas) {
			if (tempMap.containsKey(i)) {
				tempMap.put(i, tempMap.get(i) + 1);
			} else {
				tempMap.put(i, 1);
			}
		}

		int[] result = new int[] { 0, 0 };

		for (Integer key : tempMap.keySet()) {
			if (tempMap.get(key) > result[1]) {
				result[0] = key;
				result[1] = tempMap.get(key);
			}
		}

		return result;
	}

	private static int[] getMode(int[][] datas) {
		int[] result = new int[datas.length];

		for (int i = 0; i < datas.length; i++) {
			result[i] = getMode(datas[i])[0];
		}

		return result;
	}

	private static boolean checkMajorityVote(int[][] datas, double majority) {
		int[] mode = getMode(getMode(datas));
		boolean result = false;
		if (majority * datas.length <= mode[1]) {
			result = true;
		}
		return result;
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

	private static Map<Integer, Integer> groupCol(int[][] data, int k) {
		Map<Integer, Integer> dataSet = new HashMap<Integer, Integer>();
		for (int i = 0; i < data.length; i++) {
			if (data[i][k] != 0) {
				if (dataSet.containsKey(data[i][k])) {
					dataSet.put(data[i][k], dataSet.get(data[i][k]) + 1);
				} else {
					dataSet.put(data[i][k], 1);
				}
			}
		}
		return dataSet;
	}

	private static double[] priorP(int seqIndex, int[][] z, double[] w,
			int particles, double alpha) {
		int group = max(z);
		double[][] pList = new double[particles][group + 1];
		double[] colSum = new double[group + 1];
		for (int j = 0; j < group + 1; j++) {
			colSum[j] = 0.0;
		}
		for (int i = 0; i < particles; i++) {
			Map<Integer, Integer> groupMap = groupCol(z, i);
			int[] groupCol = new int[groupMap.keySet().size()];
			int temp = 0;
			for (Integer id:groupMap.keySet()) {
				groupCol[temp] = id;
				temp++;
			}
			//Integer[] groupCol = (Integer[]) groupMap.keySet().toArray();
			for (int j = 0; j < group + 1; j++) {
				pList[i][j] = w[i] * alpha / (alpha + seqIndex - 1)
						/ (group + 1 - groupCol.length);
			}
			for (int j : groupCol) {
				pList[i][j] = w[i] * groupMap.get(j) / (seqIndex - 1 + alpha);
			}
			for (int j = 0; j < group + 1; j++) {
				colSum[j] += pList[i][j];
			}
		}
		return colSum;
	}

	private static double[] startPCounts(int[] counts) {
		double[] probs = new double[counts.length / 4];
		double sum = 0.0;
		for (int i = 0; i < probs.length; i++) {
			probs[i] = counts[i * 4] + counts[i * 4 + 1] + counts[i * 4 + 2]
					+ counts[i * 4 + 3];
			sum += probs[i];
		}
		for (int i = 0; i < probs.length; i++) {
			probs[i] /= sum;
		}
		return probs;
	}

	private static double[] startPCounts(List<Integer> counts) {
		double[] probs = new double[counts.size() / 4];
		double sum = 0.0;
		for (int i = 0; i < probs.length; i++) {
			probs[i] = counts.get(i * 4) + counts.get(i * 4 + 1)
					+ counts.get(i * 4 + 2) + counts.get(i * 4 + 3);
			sum += probs[i];
		}
		for (int i = 0; i < probs.length; i++) {
			probs[i] /= sum;
		}
		return probs;
	}

	private static double[] transPCounts(int[] counts) {
		double[] probs = new double[counts.length];
		for (int i = 0; i < probs.length; i += 4) {
			double sum = counts[i] + counts[i + 1] + counts[i + 2]
					+ counts[i + 3];
			probs[i] = counts[i] / sum;
			probs[i + 1] = counts[i + 1] / sum;
			probs[i + 2] = counts[i + 2] / sum;
			probs[i + 3] = counts[i + 3] / sum;
		}
		return probs;
	}

	private static double[] transPCounts(List<Integer> counts) {
		double[] probs = new double[counts.size()];
		for (int i = 0; i < probs.length; i += 4) {
			double sum = counts.get(i) + counts.get(i + 1) + counts.get(i + 2)
					+ counts.get(i + 3);
			probs[i] = counts.get(i) / sum;
			probs[i + 1] = counts.get(i + 1) / sum;
			probs[i + 2] = counts.get(i + 2) / sum;
			probs[i + 3] = counts.get(i + 3) / sum;
		}
		return probs;
	}

	private static double[][][] startP(List<DoubleRead> doubleReadList,
			int seqIndex, int[][] z, int particles) {
		int group = max(z);
		double[][][] p = new double[particles][group + 1][];
		double[] newStartP = doubleReadList.get(seqIndex).getNewGroupStartp();
		for (int k = 0; k < particles; k++) {
			for (int j = 0; j < group; j++) {
				List<Integer> clusterIndex = new ArrayList<Integer>();
				for (int i = 0; i < z.length; i++) {
					if (z[i][k] == j) {
						clusterIndex.add(i);
					}
				}
				if (clusterIndex.size() > 0) {
					List<Integer> temp = doubleReadList
							.get(clusterIndex.get(0)).getCountList();
					for (int i = 1; i < clusterIndex.size(); i++) {
						for (int t = 0; t < temp.size(); t++) {
							temp.set(
									t,
									temp.get(t)
											+ doubleReadList
													.get(clusterIndex.get(i))
													.getCountList().get(t));
						}
					}
					p[k][j] = startPCounts(temp);
				} else {
					p[k][j] = new double[newStartP.length];
				}
			}
			p[k][group] = newStartP;
		}
		return p;
	}

	private static double[][][] transP(List<DoubleRead> doubleReadList,
			int seqIndex, int[][] z, int particles) {

		int group = max(z);
		double[] newTransP = doubleReadList.get(seqIndex).getNewGroupTransp();
		double[][][] p = new double[particles][group + 1][];
		for (int k = 0; k < particles; k++) {
			for (int j = 0; j < group; j++) {
				List<Integer> clusterIndex = new ArrayList<Integer>();
				for (int i = 0; i < z.length; i++) {
					if (z[i][k] == j) {
						clusterIndex.add(i);
					}
				}
				if (clusterIndex.size() > 0) {
					List<Integer> temp = doubleReadList
							.get(clusterIndex.get(0)).getCountList();
					for (int i = 1; i < clusterIndex.size(); i++) {
						for (int t = 0; t < temp.size(); t++) {
							temp.set(
									t,
									temp.get(t)
											+ doubleReadList
													.get(clusterIndex.get(i))
													.getCountList().get(t));
						}
					}
					p[k][j] = transPCounts(temp);
				} else {
					p[k][j] = new double[newTransP.length];
				}

			}
			p[k][group] = newTransP;
//			for (double tempP:p[k][group]) {
//				System.out.print(tempP);
//				System.out.print('\t');
//			}
//			System.out.println();
		}

		return p;
	}

	private static double[] postPTemp(List<DoubleRead> doubleReadList,
			int seqIndex, int[][] z, double[] w, int particles) {

		int group = max(z);
		double[][][] tempTransP = transP(doubleReadList, seqIndex, z, particles);
		double[][][] tempStartP = startP(doubleReadList, seqIndex, z, particles);

		double[][] p = new double[particles][group + 1];
		double[] result = new double[group + 1];

		int[] startWhere = doubleReadList.get(seqIndex).getStartWhere();

		for (int k = 0; k < particles; k++) {
			for (int j = 0; j < group + 1; j++) {
				double sumTempStartP = 0.0;
				double sumTempTransP = 0.0;
				for (int i = 0; i < tempTransP[k][j].length; i++) {
					sumTempTransP += tempTransP[k][j][i];
				}
				for (int i = 0; i < tempStartP[k][j].length; i++) {
					sumTempStartP += tempStartP[k][j][i];
				}
				if ((sumTempTransP > DefaultConstants.ZERO)
						&& (sumTempStartP > DefaultConstants.ZERO)) {
					double[] temp = tempTransP[k][j];
					double tempSum = 0.0;
					for (int i = 0; i < temp.length; i++) {
						if (temp[i] < DefaultConstants.ZERO) {
							tempSum += temp[i] * DefaultConstants.ZERO;
						} else {
							tempSum += temp[i] * Math.log(temp[i]);
						}
					}
					p[k][j] = w[k] * Math.exp(tempSum);
					for (int i = 0; i < startWhere.length; i++) {
						p[k][j] *= tempStartP[k][j][startWhere[i]];
					}
				} else {
					p[k][j] = 0.0;
				}
//				System.out.println(sumTempStartP);
//				System.out.println(sumTempTransP);
			}
		}

		for (int i = 0; i < particles; i++) {
			for (int j = 0; j < group + 1; j++) {
				result[j] += p[i][j];
			}
		}

		return result;
	}

	private static double[] postP(List<DoubleRead> doubleReadList,
			int seqIndex, int[][] z, double[] w, int particles, double alpha) {
		int group = max(z);
		double[] p = new double[group + 1];
		double[] tempPostP = postPTemp(doubleReadList, seqIndex, z, w,
				particles);
		double[] tempPriorP = priorP(seqIndex, z, w, particles, alpha);

		double sum = 0.0;
		for (int i = 0; i < p.length; i++) {
			p[i] = tempPostP[i] * tempPriorP[i];
			sum += p[i];
		}
		if (sum < DefaultConstants.ZERO) {
			sum = DefaultConstants.ZERO;
		}
		for (int i = 0; i < p.length; i++) {
			p[i] /= sum;
		}

//		System.out.println(p.toString());
		return p;
	}

	private static int[] sampleInt(double[] probList, int particles, int start) {
		int[] samples = new int[particles];
		double[] accumList = new double[probList.length];
		accumList[0] = probList[0];
		for (int i = 1; i < probList.length; i++) {
			accumList[i] = probList[i] + accumList[i - 1];
		}
		for (int i = 0; i < probList.length; i++) {
			accumList[i] = accumList[i] / accumList[probList.length - 1];
		}
		Random rand = new Random();
		for (int i = 0; i < particles; i++) {
			double r = rand.nextDouble();
			int j = 0;
			while ((j < accumList.length) && (accumList[j] < r)) {
				j++;
			}
			samples[i] = j + start;
		}
		return samples;
	}

	private static int[] clusterOneSeq(List<DoubleRead> doubleReadList,
			int seqIndex, int[][] z, double[] w, int particles, double alpha) {
		double[] posterior = postP(doubleReadList, seqIndex, z, w, particles,
				alpha);
		int[] sample = sampleInt(posterior, particles, 1);

//		System.out.println(posterior.toString());

		return sample;
	}

	private static int[][] clusterOverlapSeqs(List<DoubleRead> doubleReadList,
			Set<Integer> seqList, int[][] z, double[] w, int particles,
			double alpha) {
		int[][] tempZ = new int[seqList.size()][particles];
		int flag = 0;
		for (Integer elem : seqList) {
			tempZ[flag] = clusterOneSeq(doubleReadList, elem, z, w, particles,
					alpha);
			flag++;
		}
		return tempZ;
	}

	private static double[] updateWeights(List<DoubleRead> doubleReadList,
			int seqIndex, int[][] z, double[] w, int particles) {
		double[][][] tempTransP = transP(doubleReadList, seqIndex, z, particles);
		double[][][] tempStartP = startP(doubleReadList, seqIndex, z, particles);

		double[] p = new double[particles];

		int[] startWhere = doubleReadList.get(seqIndex).getStartWhere();

		for (int k = 0; k < particles; k++) {
			double[] temp = tempTransP[k][z[seqIndex][k]];
			double tempSum = 0.0;
			for (int i = 0; i < temp.length; i++) {
				if (temp[i] < DefaultConstants.ZERO) {
					tempSum += temp[i] * DefaultConstants.ZERO;
				} else {
					tempSum += temp[i] * Math.log(temp[i]);
				}
			}
			p[k] = w[k] * Math.exp(tempSum);
			for (int i = 0; i < startWhere.length; i++) {
				p[k] *= tempStartP[k][z[seqIndex][k]][startWhere[i]];
			}
		}
		return p;
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

		for (Integer elem : overlapList.get(0)) {
			for (int j = 0; j < params.getParticles(); j++) {
				z[elem][j] = 1;
			}
		}

		System.out.println("DirichletClusterSingle: before newGroup");
		
		newGroup(doubleReadList, params.getNeighbor());
		
		System.out.println("DirichletClusterSingle: after newGroup");

		List<Set<Integer>> rmList = new ArrayList<Set<Integer>>();

		for (int i = 1; i < overlapList.size(); i++) {
			System.out.println("Overlap: " + i);
			if (overlapList.get(i).size() > 1) {
				int[][] tempZLower = clusterOverlapSeqs(doubleReadList,
						overlapList.get(i), z, w, params.getParticles(),
						params.getAlphaLow());
				int[][] tempZUpper = clusterOverlapSeqs(doubleReadList,
						overlapList.get(i), z, w, params.getParticles(),
						params.getAlphaHigh());
				boolean isMarjorityLower = checkMajorityVote(tempZLower,
						params.getMajority());
				boolean isMarjorityUpper = checkMajorityVote(tempZUpper,
						params.getMajority());
				int[] voteLower = getMode(getMode(tempZLower));
				int[] voteUpper = getMode(getMode(tempZUpper)); 
				int majorityVoteLower = voteLower[0];
				int majorityVoteUpper = voteUpper[0];
//				System.out.print(voteLower[0]);
//				System.out.print("\t");
//				System.out.print(voteLower[1]);
//				System.out.print("\t");
//				System.out.print(voteUpper[0]);
//				System.out.print("\t");
//				System.out.println(voteUpper[1]);
				if (isMarjorityLower && isMarjorityUpper
						&& (majorityVoteLower == majorityVoteUpper)) {
					for (Integer j:overlapList.get(i)) {
						for (int k = 0; k < params.getParticles(); k++) {
							z[j][k] = majorityVoteLower;
						}
					}
//					System.out.println("Succeed");
				} else {
					rmList.add(overlapList.get(i));
//					System.out.println("Remove");
				}
			} else {

				for (Integer elem : overlapList.get(i)) {
					z[elem] = clusterOneSeq(doubleReadList, elem, z, w,
							params.getParticles(), params.getAlpha());
					w = updateWeights(doubleReadList, elem, z, w,
							params.getParticles());
					double eff = 0.0;
					for (double d : w) {
						eff += d * d;
					}
					if (eff > 1 / (params.getThreshold() * params
							.getParticles())) {
						int[] resampleIndex = sampleInt(w,
								params.getParticles(), 0);
						int[] temp = new int[params.getParticles()];
						for (int j = 0; j < params.getParticles(); j++) {
							w[j] = 1.0 / (0.0 + params.getParticles());
							temp[j] = z[elem][resampleIndex[j]];
						}
						z[elem] = temp;
					}
				}
				
//				System.out.println("Single");
			}
		}

		System.out.println("Begin to process removed seqs");
		
		for (Set<Integer> seqList : rmList) {
			int[][] tempZ = clusterOverlapSeqs(doubleReadList, seqList, z, w,
					params.getParticles(), params.getAlpha());
			int majorityVote = getMode(getMode(tempZ))[0];
			for (int j = 0; j < seqList.size(); j++) {
				for (int k = 0; k < params.getParticles(); k++) {
					z[j][k] = majorityVote;
				}
			}
		}
		
		System.out.println("Finish");

		int[] result = getMode(z);
		for (int i = 0; i < result.length; i++) {
			result[i] = result[i] + zModeLower;
		}
		return result;
	}
}
