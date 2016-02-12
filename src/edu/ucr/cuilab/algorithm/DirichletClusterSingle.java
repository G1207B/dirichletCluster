package edu.ucr.cuilab.algorithm;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

public class DirichletClusterSingle {

	private static void newGroup(List<DoubleRead> doubleReadList, Params params) {

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

	private static boolean checkMajorityVote(int[][] datas, Params params) {
		int[] mode = getMode(getMode(datas));
		boolean result = false;
		if (params.getMajority() * datas.length <= mode[1]) {
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

	// private static double[][] postPTemp(int seqIndex, Params params, int[][]
	// z) {
	// int group = max(z);
	// double[][] pList = new double[params.getParticles()][group + 1];
	// for (int i = 0; i < params.getParticles(); i++) {
	// Map<Integer, Integer> groupMap = groupCol(z, i);
	// Integer[] groupCol = (Integer[])groupMap.keySet().toArray();
	// for (int j = 0; j < group + 1; j++) {
	// pList[i][j] = params.getAlpha()
	// / (params.getAlpha() + seqIndex - 1)
	// / (group + 1 - groupCol.length);
	// }
	// for (int j:groupCol) {
	// pList[i][j] = groupMap.get(j) / (seqIndex - 1 + params.getAlpha());
	// }
	// }
	// return pList;
	// }

	private static double[] priorP(int seqIndex, Params params, int[][] z,
			int[] w) {
		int group = max(z);
		double[][] pList = new double[params.getParticles()][group + 1];
		double[] colSum = new double[group + 1];
		for (int j = 0; j < group + 1; j++) {
			colSum[j] = 0.0;
		}
		for (int i = 0; i < params.getParticles(); i++) {
			Map<Integer, Integer> groupMap = groupCol(z, i);
			Integer[] groupCol = (Integer[]) groupMap.keySet().toArray();
			for (int j = 0; j < group + 1; j++) {
				pList[i][j] = w[i] * params.getAlpha()
						/ (params.getAlpha() + seqIndex - 1)
						/ (group + 1 - groupCol.length);
			}
			for (int j : groupCol) {
				pList[i][j] = w[i] * groupMap.get(j)
						/ (seqIndex - 1 + params.getAlpha());
			}
			for (int j = 0; j < group + 1; j++) {
				colSum[j] += pList[i][j];
			}
		}
		return colSum;
	}

	private static double[] startPCounts(int[] counts, Params params) {
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

	private static double[] transPCounts(int[] counts, Params params) {
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

	private static void transP(List<DoubleRead> doubleReadList, Params params,
			int seqIndex, int[][] z) {

	}

	private static void postPTemp() {

	}

	private static void postP() {

	}

	private static int[] sampleInt(double[] probList, Params params) {
		int[] samples = new int[params.getParticles()];
		double[] accumList = new double[probList.length];
		accumList[0] = probList[0];
		for (int i = 1; i < probList.length; i++) {
			accumList[i] = probList[i] + accumList[i - 1];
		}
		Random rand = new Random();
		for (int i = 0; i < params.getParticles(); i++) {
			double r = rand.nextDouble();
			int j = 0;
			while ((accumList[j] < r) && (j < accumList.length)) {
				j++;
			}
			samples[i] = j + 1;
		}
		return samples;
	}

	private static int[] clusterOneSeq(DoubleRead doubleRead, Params params,
			int[][] z, int[] w) {
		int[] sample = new int[params.getParticles()];

		return sample;
	}

	private static int[][] clusterOverlapSeqs(List<DoubleRead> doubleReadList,
			Params params, int[][] z, int[] w) {
		int[][] tempZ = new int[doubleReadList.size()][params.getParticles()];
		for (int i = 0; i < doubleReadList.size(); i++) {
			tempZ[i] = clusterOneSeq(doubleReadList.get(i), params, z, w);
		}
		return tempZ;
	}

	public static void dirichletClusterSingle(List<DoubleRead> doubleReadList,
			List<List<Integer>> overlapList, Params params, int zModeLower) {
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

	}
}
