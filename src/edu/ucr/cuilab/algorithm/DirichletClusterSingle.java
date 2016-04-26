package edu.ucr.cuilab.algorithm;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;

public class DirichletClusterSingle {

	public static void newGroup(List<DoubleRead> doubleReadList, int neighbor) {
		Collections.sort(doubleReadList, new DoubleReadCompareGC());
		for (int i = 0; i < doubleReadList.size(); i++) {
			int nflag = 1;
			int pflag = 1;
			List<Integer> accumCountList = new ArrayList<Integer>(
					doubleReadList.get(i).getCountList());

			while ((nflag <= neighbor / 2) && (i - nflag >= 0)) {
				List<Integer> temp = doubleReadList.get(i - nflag)
						.getCountList();
				for (int j = 0; j < accumCountList.size(); j++) {
					accumCountList.set(j, accumCountList.get(j) + temp.get(j));
				}
				nflag++;
			}

			while ((pflag <= neighbor / 2)
					&& (i + pflag < doubleReadList.size())) {
				List<Integer> temp = doubleReadList.get(i + pflag)
						.getCountList();
				for (int j = 0; j < accumCountList.size(); j++) {
					accumCountList.set(j, accumCountList.get(j) + temp.get(j));
				}
				pflag++;
			}

			doubleReadList.get(i).setAccumulateCountList(accumCountList);
			doubleReadList.get(i).setNewGroupStartp(
					startPCounts(accumCountList));
			doubleReadList.get(i).setNewGroupTransp(
					transPCounts(accumCountList));
		}
		Collections.sort(doubleReadList, new DoubleReadCompareID());
	}

	public static int[] compressResult(int[] datas) {
		Map<Integer, Integer> keyMap = new TreeMap<Integer, Integer>();
		int[] result = new int[datas.length];
		int key = 1;
		for (int i : datas) {
			if (!keyMap.containsKey(i)) {
				keyMap.put(i, key);
				key++;
			}
		}
		for (int i = 0; i < datas.length; i++) {
			result[i] = keyMap.get(datas[i]);
		}
		return result;
	}

	// return {mode, length of mode}
	public static int[] getMode(int[] datas) {
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
			if (tempMap.get(key) >= result[1]) {
				result[0] = key;
				result[1] = tempMap.get(key);
			}
		}
		return result;
	}

	public static int[] getMode(int[][] datas) {
		int[] result = new int[datas.length];

		for (int i = 0; i < datas.length; i++) {
			result[i] = getMode(datas[i])[0];
		}
		return result;
	}

	public static boolean checkMajorityVote(int[][] datas, double majority) {
		return (majority * datas.length < getMode(getMode(datas))[1]);
	}

	public static List<Map<Integer, Integer>> updateGroupMapList(
			List<Map<Integer, Integer>> groupMapList, int[][] z,
			Set<Integer> seqList) {
		for (int i = 0; i < groupMapList.size(); i++) {
			Map<Integer, Integer> temp = groupMapList.get(i);
			for (Integer j : seqList) {
				if (temp.containsKey(z[j][i])) {
					temp.put(z[j][i], temp.get(z[j][i]) + 1);
				} else {
					temp.put(z[j][i], 1);
				}
			}
		}
		return groupMapList;
	}

	private static double[] priorP(List<Map<Integer, Integer>> groupMapList,
			int seqIndex, int[][] z, double[] w, int particles, double alpha,
			int completeSeqCount, int group) {
		double[][] pList = new double[particles][group + 1];
		double[] colSum = new double[group + 1];
		for (int j = 0; j < group + 1; j++) {
			colSum[j] = 0.0;
		}
		for (int i = 0; i < particles; i++) {
			Map<Integer, Integer> groupMap = groupMapList.get(i);
			int[] groupCol = new int[groupMap.keySet().size()];
			int temp = 0;
			for (Integer id : groupMap.keySet()) {
				groupCol[temp] = id;
				temp++;
			}
			for (int j = 0; j < group + 1; j++) {
				pList[i][j] = w[i] * alpha / (alpha + completeSeqCount)
						/ (group + 1 - groupCol.length);
			}
			for (int j : groupCol) {
				pList[i][j - 1] = w[i] * groupMap.get(j)
						/ (completeSeqCount + alpha);
			}
			for (int j = 0; j < group + 1; j++) {
				colSum[j] += pList[i][j];
			}
		}
		return colSum;
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

	public static List<Map<Integer, List<Integer>>> updateCountLists(
			List<DoubleRead> doubleReadList, int[][] z, Set<Integer> seqList,
			List<Map<Integer, List<Integer>>> originCountLists) {

		for (int i = 0; i < originCountLists.size(); i++) {
			Map<Integer, List<Integer>> listMap = originCountLists.get(i);
			for (Integer j : seqList) {
				if (listMap.containsKey(z[j][i])) {
					List<Integer> temp = new ArrayList<Integer>(
							listMap.get(z[j][i]));
					for (int k = 0; k < temp.size(); k++) {
						temp.set(k, temp.get(k)
								+ doubleReadList.get(j).getCountList().get(k));
					}
					listMap.put(z[j][i], temp);
				} else {
					listMap.put(z[j][i], doubleReadList.get(j).getCountList());
				}
			}
			originCountLists.set(i, listMap);
		}

		return originCountLists;
	}

	private static double[][][] startP(
			List<Map<Integer, List<Integer>>> countLists,
			List<DoubleRead> doubleReadList, int seqIndex, int[][] z,
			int particles, int group) {
		double[][][] p = new double[particles][group + 1][];
		double[] newStartP = doubleReadList.get(seqIndex).getNewGroupStartp();
		for (int k = 0; k < particles; k++) {
			for (int j = 0; j < group; j++) {
				p[k][j] = new double[newStartP.length];
			}
			for (Integer j : countLists.get(k).keySet()) {
				p[k][j - 1] = startPCounts(countLists.get(k).get(j));
			}
			p[k][group] = newStartP;
		}
		return p;
	}

	private static double[][][] transP(
			List<Map<Integer, List<Integer>>> countLists,
			List<DoubleRead> doubleReadList, int seqIndex, int[][] z,
			int particles, int group) {
		double[] newTransP = doubleReadList.get(seqIndex).getNewGroupTransp();
		double[][][] p = new double[particles][group + 1][];
		for (int k = 0; k < particles; k++) {
			for (int j = 0; j < group; j++) {
				p[k][j] = new double[newTransP.length];
			}
			for (Integer j : countLists.get(k).keySet()) {
				p[k][j - 1] = transPCounts(countLists.get(k).get(j));
			}
			p[k][group] = newTransP;
		}
		return p;
	}

	private static double[] postPTemp(
			List<Map<Integer, List<Integer>>> countLists,
			List<DoubleRead> doubleReadList, int seqIndex, int[][] z,
			double[] w, int particles, int group) {

		double[][][] tempTransP = transP(countLists, doubleReadList, seqIndex,
				z, particles, group);
		double[][][] tempStartP = startP(countLists, doubleReadList, seqIndex,
				z, particles, group);

		double[][] p = new double[particles][group + 1];
		double[] result = new double[group + 1];

		int[] startWhere = doubleReadList.get(seqIndex).getStartWhere();
		List<Integer> countList = doubleReadList.get(seqIndex).getCountList();
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
							tempSum += countList.get(i) * DefaultConstants.ZERO;
						} else {
							tempSum += countList.get(i) * Math.log(temp[i]);
						}
					}
					p[k][j] = w[k] * Math.exp(tempSum);

					for (int i = 0; i < startWhere.length; i++) {
						p[k][j] *= tempStartP[k][j][startWhere[i]] + Double.MIN_NORMAL;
					}
				} else {
					p[k][j] = Double.MIN_NORMAL;
				}
			}
		}
		for (int i = 0; i < particles; i++) {
			for (int j = 0; j < group + 1; j++) {
				result[j] += p[i][j];
			}
		}
		return result;
	}

	private static double[] postP(List<Map<Integer, Integer>> groupCount,
			List<Map<Integer, List<Integer>>> countLists,
			List<DoubleRead> doubleReadList, int seqIndex, int[][] z,
			double[] w, int particles, double alpha, int completeSeqCount, int group) {
		double[] p = new double[group + 1];
		double[] tempPostP = postPTemp(countLists, doubleReadList, seqIndex, z,
				w, particles, group);

		double[] tempPriorP = priorP(groupCount, seqIndex, z, w, particles,
				alpha, completeSeqCount, group);

		double sum = 0.0;
		for (int i = 0; i < p.length; i++) {
			p[i] = tempPostP[i] * tempPriorP[i] + Double.MIN_NORMAL;
			sum += p[i];
		}
		for (int i = 0; i < p.length; i++) {
			p[i] /= sum;
		}
		return p;
	}

	public static int[] sampleInt(double[] probList, int particles, int start) {
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

	private static int[] clusterOneSeq(List<Map<Integer, Integer>> groupCount,
			List<Map<Integer, List<Integer>>> countLists,
			List<DoubleRead> doubleReadList, int seqIndex, int[][] z,
			double[] w, int particles, double alpha, int completeSeqCount, int group) {
		double[] posterior = postP(groupCount, countLists, doubleReadList,
				seqIndex, z, w, particles, alpha, completeSeqCount, group);
		int[] sample = sampleInt(posterior, particles, 1);
		return sample;
	}

	public static int[][] clusterOverlapSeqs(
			List<Map<Integer, Integer>> groupCount,
			List<Map<Integer, List<Integer>>> countLists,
			List<DoubleRead> doubleReadList, Set<Integer> seqList, int[][] z,
			double[] w, int particles, double alpha, int completeSeqCount, int group) {
		int[][] tempZ = new int[seqList.size()][particles];
		int flag = 0;
		for (Integer elem : seqList) {
			tempZ[flag] = clusterOneSeq(groupCount, countLists, doubleReadList,
					elem, z, w, particles, alpha, completeSeqCount, group);
			flag++;
		}
		return tempZ;
	}

	public static double[] updateWeights(
			List<Map<Integer, List<Integer>>> countLists,
			List<DoubleRead> doubleReadList, int seqIndex, int[][] z,
			double[] w, int particles, int group) {
		double[][][] tempTransP = transP(countLists, doubleReadList, seqIndex,
				z, particles, group);
		double[][][] tempStartP = startP(countLists, doubleReadList, seqIndex,
				z, particles, group);

		double[] p = new double[particles];

		int[] startWhere = doubleReadList.get(seqIndex).getStartWhere();

		List<Integer> countList = doubleReadList.get(seqIndex).getCountList();
		for (int k = 0; k < particles; k++) {
			double[] temp = tempTransP[k][z[seqIndex][k] - 1];
			double tempSum = 0.0;
			for (int i = 0; i < temp.length; i++) {
				if (temp[i] < DefaultConstants.ZERO) {
					tempSum += countList.get(i) * DefaultConstants.ZERO;
				} else {
					tempSum += countList.get(i) * Math.log(temp[i]);
				}
			}
			p[k] = w[k] * Math.exp(tempSum);
			for (int i = 0; i < startWhere.length; i++) {
				p[k] *= tempStartP[k][z[seqIndex][k] - 1][startWhere[i]];
			}
			p[k] += Double.MIN_NORMAL;
		}

		double sum = 0.0;
		for (double elem : p) {
			sum += elem;
		}
		for (int i = 0; i < p.length; i++) {
			p[i] /= sum;
		}
		return p;
	}
}
