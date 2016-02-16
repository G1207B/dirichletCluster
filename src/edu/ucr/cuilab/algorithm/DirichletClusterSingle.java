package edu.ucr.cuilab.algorithm;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

public class DirichletClusterSingle {

	private static void newGroup(List<DoubleRead> doubleReadList, Params params) {
		Collections.sort(doubleReadList, new DoubleReadCompareGC());
		for (int i = 0; i < doubleReadList.size(); i++) {
			int nflag = 1;
			int pflag = 1;
			List<Integer> accumCountList = doubleReadList.get(i).getCountList();

			while ((nflag <= params.getNeighbor() / 2) && (i - nflag >= 0)) {
				List<Integer> temp = doubleReadList.get(i - nflag)
						.getCountList();
				for (int j = 0; j < accumCountList.size(); j++) {
					accumCountList.set(j, accumCountList.get(j) + temp.get(j));
				}
				nflag++;
			}

			while ((pflag <= params.getNeighbor() / 2)
					&& (i + nflag < doubleReadList.size())) {
				List<Integer> temp = doubleReadList.get(i + pflag)
						.getCountList();
				for (int j = 0; j < accumCountList.size(); j++) {
					accumCountList.set(j, accumCountList.get(j) + temp.get(j));
				}
				pflag++;
			}

			doubleReadList.get(i).setNewGroupStartp(
					startPCounts(accumCountList, params));
			doubleReadList.get(i).setNewGroupTransp(
					transPCounts(accumCountList, params));
		}
		Collections.sort(doubleReadList, new DoubleReadCompareID());
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

	private static double[] priorP(int seqIndex, Params params, int[][] z,
			double[] w) {
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

	private static double[] startPCounts(List<Integer> counts, Params params) {
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

	private static double[] transPCounts(List<Integer> counts, Params params) {
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
			Params params, int seqIndex, int[][] z) {
		int group = max(z) + 1;
		double[][][] p = new double[params.getParticles()][group + 1][];
		double[] newStartP = doubleReadList.get(seqIndex).getNewGroupStartp();
		for (int k = 0; k < params.getParticles(); k++) {
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
					p[k][j] = startPCounts(temp, params);
				} else {
					p[k][j] = new double[newStartP.length];
				}
			}
			p[k][group] = newStartP;
		}
		return p;
	}

	private static double[][][] transP(List<DoubleRead> doubleReadList,
			Params params, int seqIndex, int[][] z) {

		int group = max(z) + 1;
		double[] newTransP = doubleReadList.get(seqIndex).getNewGroupTransp();
		double[][][] p = new double[params.getParticles()][group + 1][];
		for (int k = 0; k < params.getParticles(); k++) {
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
					p[k][j] = transPCounts(temp, params);
				} else {
					p[k][j] = new double[newTransP.length];
				}

			}
			p[k][group] = newTransP;
		}

		return p;
	}

	private static double[] postPTemp(List<DoubleRead> doubleReadList,
			Params params, int seqIndex, int[][] z, double[] w) {

		int group = max(z) + 1;
		double[][][] tempTransP = transP(doubleReadList, params, seqIndex, z);
		double[][][] tempStartP = startP(doubleReadList, params, seqIndex, z);

		double[][] p = new double[params.getParticles()][group + 1];
		double[] result = new double[group + 1];

		int[] startWhere = doubleReadList.get(seqIndex).getStartWhere();

		for (int k = 0; k < params.getParticles(); k++) {
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
			}
		}

		for (int i = 0; i < params.getParticles(); i++) {
			for (int j = 0; j < group + 1; j++) {
				result[j] += p[i][j];
			}
		}

		return result;
	}

	private static double[] postP(List<DoubleRead> doubleReadList,
			Params params, int seqIndex, int[][] z, double[] w) {
		int group = max(z) + 1;
		double[] p = new double[group + 1];
		double[] tempPostP = postPTemp(doubleReadList, params, seqIndex, z, w);
		double[] tempPriorP = priorP(seqIndex, params, z, w);
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
		return p;
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

	private static int[] clusterOneSeq(List<DoubleRead> doubleReadList,
			Params params, int seqIndex, int[][] z, double[] w) {
		double[] posterior = postP(doubleReadList, params, seqIndex, z, w);
		int[] sample = sampleInt(posterior, params);

		return sample;
	}

	private static int[][] clusterOverlapSeqs(List<DoubleRead> doubleReadList,
			Params params, int[][] z, double[] w) {
		int[][] tempZ = new int[doubleReadList.size()][params.getParticles()];
		for (int i = 0; i < doubleReadList.size(); i++) {
			tempZ[i] = clusterOneSeq(doubleReadList, params, i, z, w);
		}
		return tempZ;
	}

	private static double[] updateWeights(List<DoubleRead> doubleReadList,
			Params params, int seqIndex, int[][] z, double[] w) {
		double[][][] tempTransP = transP(doubleReadList, params, seqIndex, z);
		double[][][] tempStartP = startP(doubleReadList, params, seqIndex, z);

		double[] p = new double[params.getParticles()];

		int[] startWhere = doubleReadList.get(seqIndex).getStartWhere();

		for (int k = 0; k < params.getParticles(); k++) {
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

	private static List<List<Integer>> updateOverlap(
			List<DoubleRead> doubleReadList, List<List<Integer>> overlapList) {
		return new ArrayList<List<Integer>>();
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
